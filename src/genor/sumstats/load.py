"""Load a user-supplied summary statistics file into GeNor's standard table.

Standard columns (GWAS-SSF names): variant_id, chromosome, base_pair_location,
effect_allele, other_allele, beta, standard_error, p_value, z, n, and when present
effect_allele_frequency, info, n_cases, n_controls.

A loaded trait is a directory holding ``sumstats.parquet`` plus ``sumstats.json``
(trait metadata: build, ancestry, N, case/control counts, column mapping, source).
Other commands take that directory as input.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from genor.errors import GenorError, NeedsUserDecision, Warnings
from genor.guards import norm_ancestry, norm_build, refuse_individual_level
from genor.sumstats.columns import detect, normalise

STANDARD = ["variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "p_value", "z", "n", "effect_allele_frequency", "info",
            "n_cases", "n_controls"]

CHROM_MAP = {"X": "23", "Y": "24", "XY": "25", "MT": "26", "M": "26"}


@dataclass
class TraitMeta:
    trait: str
    build: str | None = None
    ancestry: str | None = None
    n: float | None = None
    n_cases: float | None = None
    n_controls: float | None = None
    trait_type: str = "quantitative"  # quantitative | binary
    sample_prevalence: float | None = None
    population_prevalence: float | None = None
    source: str | None = None
    source_sha256: str | None = None
    columns: dict[str, str] = field(default_factory=dict)
    effect_type: str = "beta"
    build_inferred: bool = False
    n_variants: int = 0
    qc: dict[str, Any] | None = None

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, default=str)


def _open_text(path: Path) -> io.TextIOBase:
    name = path.name.lower()
    if name.endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    if name.endswith(".zip"):
        zf = zipfile.ZipFile(path)
        members = [m for m in zf.namelist() if not m.endswith("/") and "__MACOSX" not in m]
        if len(members) != 1:
            raise GenorError(f"{path.name} contains {len(members)} files; expected one.",
                             code="zip_multiple_members", hint="Unzip it and pass the file you mean.")
        return io.TextIOWrapper(zf.open(members[0]), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _sniff(path: Path) -> tuple[str | None, int]:
    """Return (separator, number of leading comment lines). None means whitespace."""
    skip = 0
    with _open_text(path) as f:
        lines = []
        for line in f:
            if line.startswith("##"):
                skip += 1
                continue
            lines.append(line)
            if len(lines) >= 5:
                break
    if not lines:
        raise GenorError(f"{path.name} is empty", code="empty_file")
    head = lines[0]
    if "\t" in head:
        return "\t", skip
    if "," in head and head.count(",") >= 3:
        return ",", skip
    try:
        d = csv.Sniffer().sniff("".join(lines), delimiters=" ;|")
        if d.delimiter != " ":
            return d.delimiter, skip
    except csv.Error:
        pass
    return None, skip


def read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith(".parquet"):
        return pd.read_parquet(path)
    if ".vcf" in name:
        return read_gwas_vcf(path)
    sep, skip = _sniff(path)
    kw: dict[str, Any] = {"sep": sep} if sep else {"sep": r"\s+"}
    with _open_text(path) as f:
        return pd.read_csv(f, skiprows=skip, dtype=str, na_values=["NA", "NaN", "nan", ".", "", "-"],
                           keep_default_na=False, low_memory=False, **kw)


def read_gwas_vcf(path: Path) -> pd.DataFrame:
    """Read a single-trait GWAS-VCF (MRC IEU OpenGWAS format)."""
    rows = []
    opener = gzip.open if path.name.lower().endswith(("gz", "bgz")) else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            fmt = p[8].split(":")
            val = dict(zip(fmt, p[9].split(":"), strict=False))
            alts = p[4].split(",")
            if len(alts) != 1:
                continue
            rows.append({
                "variant_id": val.get("ID", p[2]),
                "chromosome": p[0], "base_pair_location": p[1],
                "effect_allele": alts[0], "other_allele": p[3],
                "beta": val.get("ES"), "standard_error": val.get("SE"),
                "neg_log10_p": val.get("LP"), "effect_allele_frequency": val.get("AF"),
                "n": val.get("SS"),
            })
    return pd.DataFrame(rows).replace({".": None})


def _num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _parse_variant_ids(df: pd.DataFrame, w: Warnings) -> pd.DataFrame:
    """Recover chromosome/position (and alleles) from IDs like 1:12345:A:G or chr1_12345_A_G."""
    if "variant_id" not in df or ("chromosome" in df and "base_pair_location" in df):
        return df
    parts = df["variant_id"].astype(str).str.replace("chr", "", regex=False).str.split(r"[:_]", expand=True)
    if parts.shape[1] >= 2 and pd.to_numeric(parts[1], errors="coerce").notna().mean() > 0.9:
        df["chromosome"] = parts[0]
        df["base_pair_location"] = parts[1]
        if parts.shape[1] >= 4 and "effect_allele" not in df:
            df["other_allele"], df["effect_allele"] = parts[2], parts[3]
            w.add("alleles_from_id", "Alleles were taken from the variant ID as REF:ALT with ALT as "
                  "the effect allele. Check this against the file's documentation.", "serious")
        w.add("position_from_id", "Chromosome and position were parsed from the variant ID.", "info")
    return df


def standardise(raw: pd.DataFrame, overrides: dict[str, str] | None, w: Warnings,
                n: float | None = None, n_cases: float | None = None,
                n_controls: float | None = None, effect_type: str | None = None,
                ) -> tuple[pd.DataFrame, dict[str, str], str]:
    try:
        mapping, notes = detect(list(raw.columns), overrides)
    except KeyError as e:
        raise GenorError(str(e), code="bad_column_override") from e
    for note in notes:
        w.add("column_note", note, "info")
    df = raw[list(mapping)].rename(columns=mapping).copy()
    norm = {normalise(c): c for c in raw.columns}
    if "a1" in norm and "ref" in norm and "alt" in norm and "other_allele" not in df:
        a1 = raw[norm["a1"]].astype(str).str.upper()
        ref = raw[norm["ref"]].astype(str).str.upper()
        alt = raw[norm["alt"]].astype(str).str.upper()
        df["other_allele"] = np.where(a1 == alt, ref, alt)
        if "a1freq" in norm and "effect_allele_frequency" not in df:
            df["effect_allele_frequency"] = raw[norm["a1freq"]]
    df = _parse_variant_ids(df, w)

    missing_core = [c for c in ("effect_allele", "other_allele") if c not in df]
    if missing_core:
        raise GenorError(
            f"Could not find columns for: {', '.join(missing_core)}",
            code="columns_not_detected",
            hint=("Name them explicitly, e.g. --col effect_allele=ALT --col other_allele=REF. "
                  f"File columns: {', '.join(map(str, raw.columns))}"),
            detected=mapping,
        )
    for c in ("effect_allele", "other_allele"):
        df[c] = df[c].astype(str).str.upper()
    if "chromosome" in df:
        ch = df["chromosome"].astype(str).str.upper().str.replace("CHR", "", regex=False)
        df["chromosome"] = ch.replace(CHROM_MAP)
    for c in ("beta", "odds_ratio", "standard_error", "p_value", "neg_log10_p", "z", "n",
              "effect_allele_frequency", "info", "n_cases", "n_controls", "base_pair_location"):
        if c in df:
            df[c] = _num(df[c])

    # Effect sizes ----------------------------------------------------------
    etype = effect_type or "beta"
    if "beta" not in df and "odds_ratio" in df:
        etype = "or"
    if etype == "or":
        src = df["odds_ratio"] if "odds_ratio" in df else df["beta"]
        bad = (src <= 0).sum()
        if bad:
            w.add("nonpositive_or", f"{bad} odds ratios <= 0 set to missing.", "warning")
        df["beta"] = np.log(src.where(src > 0))
        w.add("or_to_logor", "Odds ratios converted to log-odds (beta = ln OR).", "info")
    elif "beta" in df and effect_type is None:
        b = df["beta"].dropna()
        if len(b) > 100 and (b > 0).all() and 0.8 < b.median() < 1.25:
            raise NeedsUserDecision(
                "The 'beta' column is always positive with median near 1, which looks like odds ratios.",
                question="Is the effect column an odds ratio (OR) or a beta / log-odds?",
                options=["or", "beta"],
            )
    df = df.drop(columns=[c for c in ("odds_ratio",) if c in df])

    # P-values ---------------------------------------------------------------
    if "p_value" not in df and "neg_log10_p" in df:
        df["p_value"] = np.power(10.0, -df["neg_log10_p"])
    df = df.drop(columns=[c for c in ("neg_log10_p",) if c in df])

    has_b, has_se, has_p, has_z = ("beta" in df), ("standard_error" in df), ("p_value" in df), ("z" in df)
    if has_b and has_se:
        z = df["beta"] / df["standard_error"]
        if has_z:
            w.add("z_recomputed", "Z recomputed as beta/SE (file Z column ignored).", "info")
        df["z"] = z
    elif has_z and not has_b:
        pass
    elif has_b and has_p:
        df["z"] = np.sign(df["beta"]) * stats.norm.isf(df["p_value"].clip(lower=1e-300) / 2)
        df["standard_error"] = (df["beta"] / df["z"]).abs()
        w.add("se_from_p", "SE derived from beta and P (no SE column).", "warning")
    else:
        raise GenorError("Need beta+SE, beta+P, or Z to derive signed statistics.",
                         code="no_signed_statistic",
                         hint="Check that the file has an effect column (beta/OR/Z).")
    if "p_value" not in df:
        df["p_value"] = 2 * stats.norm.sf(df["z"].abs())
        w.add("p_from_z", "P computed from Z.", "info")
    elif has_b and has_se:
        # consistency check: P from Z should agree with the file's P
        sub = df[["z", "p_value"]].dropna().head(20000)
        if len(sub) > 100:
            p2 = 2 * stats.norm.sf(sub["z"].abs())
            ok = np.isclose(np.log10(p2.clip(1e-300)), np.log10(sub["p_value"].clip(1e-300)), atol=0.5)
            if ok.mean() < 0.8:
                w.add("p_z_inconsistent", f"P and beta/SE disagree for {100 * (1 - ok.mean()):.0f}% of "
                      "variants. The SE may be on another scale, or P may be GC-corrected.", "serious")

    # Sample size -------------------------------------------------------------
    if n_cases is not None:
        df["n_cases"] = n_cases
    if n_controls is not None:
        df["n_controls"] = n_controls
    if n is not None:
        if "n" in df and df["n"].notna().any():
            w.add("n_overridden", f"--n {n:g} replaced the per-variant N column from the file "
                  f"(median {df['n'].median():g}).", "warning")
        df["n"] = n
    elif "n" not in df and "n_cases" in df and "n_controls" in df:
        df["n"] = df["n_cases"] + df["n_controls"]
    cols = [c for c in STANDARD if c in df]
    return df[cols], mapping, etype


# ---------------------------------------------------------------------------
# Build inference


def _markers() -> pd.DataFrame:
    from importlib import resources

    with resources.files("genor.data").joinpath("build_markers.tsv.gz").open("rb") as f:
        return pd.read_csv(f, sep="\t", compression="gzip", dtype={"chromosome": str})


def infer_build(df: pd.DataFrame) -> tuple[str | None, dict[str, float]]:
    if "chromosome" not in df or "base_pair_location" not in df:
        return None, {}
    m = _markers()
    if "variant_id" in df:
        sub = df[["variant_id", "chromosome", "base_pair_location"]].merge(
            m, left_on="variant_id", right_on="rsid", suffixes=("", "_m"))
        if len(sub) >= 50:
            f37 = float((sub["base_pair_location"] == sub["pos37"]).mean())
            f38 = float((sub["base_pair_location"] == sub["pos38"]).mean())
            return _decide(f37, f38, len(sub))
    keys = df["chromosome"].astype(str) + ":" + df["base_pair_location"].astype("Int64").astype(str)
    s = set(keys)
    k37 = m["chromosome"].astype(str) + ":" + m["pos37"].astype(str)
    k38 = m["chromosome"].astype(str) + ":" + m["pos38"].astype(str)
    f37, f38 = float(k37.isin(s).mean()), float(k38.isin(s).mean())
    return _decide(f37, f38, len(m))


def _decide(f37: float, f38: float, n: int) -> tuple[str | None, dict[str, float]]:
    ev = {"match_grch37": round(f37, 3), "match_grch38": round(f38, 3), "markers": n}
    if f37 > 0.5 and f37 > 3 * max(f38, 0.01):
        return "GRCh37", ev
    if f38 > 0.5 and f38 > 3 * max(f37, 0.01):
        return "GRCh38", ev
    if max(f37, f38) > 0.05:  # sparse files: few markers present, but they agree
        if f37 > 5 * max(f38, 0.001):
            return "GRCh37", ev
        if f38 > 5 * max(f37, 0.001):
            return "GRCh38", ev
    return None, ev


# ---------------------------------------------------------------------------


def load(path: Path, trait: str, out_dir: Path, *, build: str | None = None,
         ancestry: str | None = None, n: float | None = None, n_cases: float | None = None,
         n_controls: float | None = None, population_prevalence: float | None = None,
         effect_type: str | None = None, columns: dict[str, str] | None = None,
         w: Warnings | None = None) -> tuple[pd.DataFrame, TraitMeta, Warnings]:
    w = w if w is not None else Warnings()
    path = Path(path)
    if not path.exists():
        raise GenorError(f"{path} does not exist", code="file_not_found")
    refuse_individual_level(path)
    raw = read_table(path)
    if len(raw) == 0:
        raise GenorError(f"{path.name} has no rows", code="empty_file")
    df, mapping, etype = standardise(raw, columns, w, n=n, n_cases=n_cases,
                                     n_controls=n_controls, effect_type=effect_type)

    meta = TraitMeta(trait=trait, source=str(path.resolve()), columns=mapping, effect_type=etype,
                     ancestry=norm_ancestry(ancestry), population_prevalence=population_prevalence)
    b = norm_build(build)
    inferred, evidence = infer_build(df)
    if b is None:
        b = inferred
        meta.build_inferred = b is not None
        if b is None:
            w.add("build_unknown", "Genome build could not be inferred from positions.", "serious",
                  **evidence)
        else:
            w.add("build_inferred", f"Genome build inferred as {b}.", "info", **evidence)
    elif inferred and inferred != b:
        w.add("build_conflict", f"You said {b} but positions match {inferred}.", "serious", **evidence)
    meta.build = b

    if "n" in df and df["n"].notna().any():
        meta.n = float(df["n"].median())
    if "n_cases" in df and df["n_cases"].notna().any():
        meta.n_cases = float(df["n_cases"].median())
        meta.n_controls = float(df["n_controls"].median()) if "n_controls" in df else None
        meta.trait_type = "binary"
        if meta.n_cases and meta.n_controls:
            meta.sample_prevalence = meta.n_cases / (meta.n_cases + meta.n_controls)
    elif etype == "or":
        meta.trait_type = "binary"
    if meta.n is None:
        w.add("n_missing", "No sample size found. LDSC and some MR diagnostics need N.", "serious")
    if meta.ancestry is None:
        w.add("ancestry_unknown", "GWAS ancestry not given. Methods that use an LD panel will "
              "ask for it.", "warning")

    meta.n_variants = len(df)
    out_dir.mkdir(parents=True, exist_ok=True)
    return df, meta, w


def save_trait(df: pd.DataFrame, meta: TraitMeta, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    meta.n_variants = len(df)
    df.reset_index(drop=True).to_parquet(out_dir / "sumstats.parquet", index=False)
    (out_dir / "sumstats.json").write_text(meta.to_json())


def open_trait(path: Path | str) -> tuple[pd.DataFrame, TraitMeta]:
    """Open a trait directory written by ``genor sumstats load`` (or a raw parquet)."""
    p = Path(path)
    if p.is_dir():
        pq, js = p / "sumstats.parquet", p / "sumstats.json"
        if not pq.exists():
            raise GenorError(f"{p} has no sumstats.parquet", code="not_a_trait_dir",
                             hint="Load the file first: genor sumstats load <file> --trait NAME --out DIR")
        meta = TraitMeta(**json.loads(js.read_text())) if js.exists() else TraitMeta(trait=p.name)
        return pd.read_parquet(pq), meta
    raise GenorError(f"{p} is not a trait directory", code="not_a_trait_dir",
                     hint="Pass the --out directory of 'genor sumstats load', not the raw file.")
