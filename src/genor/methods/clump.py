"""LD clumping.

``plink`` mode runs PLINK 2 ``--clump`` per chromosome against a 1000 Genomes panel.
``distance`` mode keeps the strongest variant per window with no LD information; it is
a labelled fallback for when the panel is unavailable, never a silent substitute.

Defaults:
- loci (``--preset loci``): p1 5e-8, r2 0.1, 250 kb  (FUMA-like independent signals)
- MR instruments (``--preset mr``): p1 5e-8, r2 0.001, 10,000 kb (TwoSampleMR defaults)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from genor import backends, refdata
from genor.errors import GenorError, Warnings
from genor.guards import check_ancestry, check_build
from genor.sumstats.load import TraitMeta

PRESETS = {
    "loci": {"p1": 5e-8, "r2": 0.1, "kb": 250},
    "mr": {"p1": 5e-8, "r2": 0.001, "kb": 10_000},
}

PANELS = {"EUR": "g1000-eur", "EAS": "g1000-eas"}


def distance_clump(df: pd.DataFrame, p1: float, kb: int) -> pd.DataFrame:
    d = df[df["p_value"] < p1].sort_values("p_value")
    kept: list[int] = []
    taken: dict[str, list[int]] = {}
    for idx, row in d.iterrows():
        c = str(row["chromosome"])
        pos = int(row["base_pair_location"])
        if any(abs(pos - q) <= kb * 1000 for q in taken.get(c, [])):
            continue
        taken.setdefault(c, []).append(pos)
        kept.append(idx)
    return df.loc[kept].sort_values(["chromosome", "base_pair_location"])


def panel_for(meta: TraitMeta, panel: str | None) -> str:
    if panel:
        return panel
    anc = meta.ancestry
    if anc in PANELS:
        return PANELS[anc]
    check_ancestry(anc, "EUR", "the default LD panel (1000G EUR)")
    return "g1000-eur"


def plink_clump(df: pd.DataFrame, meta: TraitMeta, work: Path, p1: float, r2: float, kb: int,
                panel: str | None, w: Warnings) -> tuple[pd.DataFrame, dict[str, Any]]:
    panel = panel_for(meta, panel)
    e = refdata.entry(panel)
    check_ancestry(meta.ancestry, e["ancestry"], f"the LD panel '{panel}'")
    check_build(meta.build, e["build"], f"the LD panel '{panel}'")
    prefix = refdata.plink_prefix(panel)
    exe = backends.require("plink2")
    cand = df[df["p_value"] < p1].copy()
    if "variant_id" not in cand or not cand["variant_id"].astype(str).str.startswith("rs").any():
        raise GenorError("PLINK clumping against 1000G matches by rsID; this trait has none.",
                         code="clump_needs_rsid", hint="Use --mode distance, or add rsIDs.")
    work.mkdir(parents=True, exist_ok=True)
    kept_ids: list[str] = []
    missing = 0
    for c, g in cand.groupby(cand["chromosome"].astype(str)):
        if not Path(f"{prefix}{c}.bed").exists():
            missing += len(g)
            continue
        f = work / f"cand.{c}.tsv"
        # include sub-threshold SNPs on this chromosome only if needed by p2 (we use p2 = p1)
        g[["variant_id", "p_value"]].rename(columns={"variant_id": "SNP", "p_value": "P"}).to_csv(
            f, sep="\t", index=False)
        out = work / f"clump.{c}"
        backends.run("plink2", ["--bfile", f"{prefix}{c}", "--clump", str(f),
                                "--clump-id-field", "SNP", "--clump-p-field", "P",
                                "--clump-p1", str(p1), "--clump-p2", str(p1),
                                "--clump-r2", str(r2), "--clump-kb", str(kb),
                                "--out", str(out)], log=work / f"plink.{c}.log")
        res = Path(f"{out}.clumps")
        if res.exists():
            cl = pd.read_csv(res, sep=r"\s+")
            kept_ids.extend(cl["ID"].astype(str))
    in_panel = cand["variant_id"].astype(str).isin(kept_ids)
    absent = cand[~in_panel]
    # PLINK silently skips variants absent from the panel: report them
    info = {"panel": panel, "plink2": str(exe), "n_candidates": int(len(cand)),
            "n_index": int(in_panel.sum())}
    if len(absent):
        info["n_candidates_removed_by_clumping_or_absent_from_panel"] = int(len(absent))
    if missing:
        w.add("chromosome_not_in_panel", f"{missing} candidates on chromosomes without panel "
              "files (e.g. X) were skipped.", "info")
    out_df = cand[in_panel].sort_values(["chromosome", "base_pair_location"])
    return out_df, info


def clump(df: pd.DataFrame, meta: TraitMeta, work: Path, *, mode: str = "plink",
          preset: str = "loci", p1: float | None = None, r2: float | None = None,
          kb: int | None = None, panel: str | None = None,
          w: Warnings | None = None) -> tuple[pd.DataFrame, dict[str, Any]]:
    w = w if w is not None else Warnings()
    pr = PRESETS[preset]
    p1, r2, kb = p1 or pr["p1"], r2 if r2 is not None else pr["r2"], kb or pr["kb"]
    params = {"mode": mode, "preset": preset, "p1": p1, "r2": r2, "kb": kb}
    if mode == "distance":
        out = distance_clump(df, p1, kb)
        w.add("distance_clumping", f"Distance-only pruning ({kb} kb window), no LD information. "
              "Index variants may still be in LD; results are approximate.", "serious")
        return out, {**params, "n_index": int(len(out))}
    out, info = plink_clump(df, meta, work, p1, r2, kb, panel, w)
    return out, {**params, **info}
