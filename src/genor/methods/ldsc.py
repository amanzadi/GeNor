# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""LD score regression: SNP heritability and genetic correlation.

Data handling mirrors ``munge_sumstats.py`` + ``ldsc.py --h2/--rg`` from the
reference implementation; the regression itself is the vendored LDSC code in
``genor.methods._ldsc``. Defaults: 200 jackknife blocks, two-step estimator with
cutoff 30 for single-annotation models, HapMap3 SNPs, MHC excluded via the weights.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from genor import refdata
from genor.errors import GenorError, Warnings
from genor.guards import check_ancestry
from genor.methods._ldsc import regressions as reg
from genor.sumstats.load import TraitMeta
from genor.sumstats.qc import is_palindromic, is_snv

DEFAULT_REF = "ldsc-ldscores-eur"
DEFAULT_W = "ldsc-weights-eur"


@dataclass
class LDRef:
    ref_ld: pd.DataFrame  # SNP, CHR, BP, L2
    w_ld: pd.DataFrame  # SNP, w
    M: float
    ancestry: str
    build: str


def _read_prefix(prefix: str, col: str) -> pd.DataFrame:
    parts = []
    for c in range(1, 23):
        f = Path(f"{prefix}{c}.l2.ldscore.gz")
        if not f.exists():
            raise GenorError(f"Missing LD score file {f}", code="refdata_corrupt")
        d = pd.read_csv(f, sep="\t", usecols=["CHR", "SNP", "BP", "L2"])
        parts.append(d.rename(columns={"L2": col}))
    return pd.concat(parts, ignore_index=True)


def _read_m(prefix: str) -> float:
    total = 0.0
    for c in range(1, 23):
        f = Path(f"{prefix}{c}.l2.M_5_50")
        if not f.exists():
            raise GenorError(f"Missing {f} (M_5_50).", code="refdata_no_m",
                             hint="Use an LD score set that ships M_5_50 files as --ref.")
        total += float(f.read_text().split()[0])
    return total


def load_reference(ref: str = DEFAULT_REF, w: str = DEFAULT_W) -> LDRef:
    rp, wp = refdata.ldscore_prefix(ref), refdata.ldscore_prefix(w)
    ref_ld = _read_prefix(rp, "L2")
    w_ld = _read_prefix(wp, "w")[["SNP", "w"]]
    e = refdata.entry(ref)
    return LDRef(ref_ld, w_ld, _read_m(rp), e["ancestry"], e["build"])


def munge(df: pd.DataFrame, meta: TraitMeta, ld: LDRef, w: Warnings, label: str) -> pd.DataFrame:
    """Reduce a GeNor table to LDSC's SNP/A1/A2/Z/N with munge_sumstats.py's filters."""
    d = df.copy()
    if "n" not in d or d["n"].isna().all():
        raise GenorError(f"{label}: no sample size; LDSC needs N.", code="n_missing",
                         hint="Reload with --n (or --n-cases/--n-controls).")
    if "variant_id" in d and d["variant_id"].astype(str).str.startswith("rs").mean() > 0.5:
        d["SNP"] = d["variant_id"].astype(str)
    elif meta.build == ld.build and "chromosome" in d:
        key = ld.ref_ld.assign(k=ld.ref_ld["CHR"].astype(str) + ":" + ld.ref_ld["BP"].astype(str))
        k = d["chromosome"].astype(str) + ":" + d["base_pair_location"].astype("Int64").astype(str)
        d["SNP"] = k.map(dict(zip(key["k"], key["SNP"], strict=True)))
        w.add("rsid_from_position", f"{label}: rsIDs assigned from {ld.build} positions.", "info")
    else:
        raise GenorError(f"{label}: no rsIDs and positions are not {ld.build}.",
                         code="cannot_match_ldscores",
                         hint=f"Lift over to {ld.build} (genor sumstats liftover) or provide rsIDs.")
    a1, a2 = d["effect_allele"], d["other_allele"]
    keep = is_snv(a1, a2) & ~is_palindromic(a1, a2) & d["z"].notna() & d["SNP"].notna() & d["n"].notna()
    if "info" in d and d["info"].notna().any():
        keep &= ~(d["info"] < 0.9)
    if "effect_allele_frequency" in d and d["effect_allele_frequency"].notna().any():
        f = d["effect_allele_frequency"]
        keep &= ~(np.minimum(f, 1 - f) <= 0.01)
    d = d[keep]
    # munge_sumstats: drop SNPs with N < 90th percentile / 1.5
    n90 = d["n"].quantile(0.9)
    low_n = d["n"] < n90 / 1.5
    if low_n.any():
        w.add("ldsc_low_n_removed", f"{label}: {int(low_n.sum())} SNPs with N < N90/1.5 removed "
              "(munge_sumstats default).", "info")
    d = d[~low_n].drop_duplicates("SNP")
    out = pd.DataFrame({"SNP": d["SNP"], "A1": a1.loc[d.index], "A2": a2.loc[d.index],
                        "Z": d["z"], "N": d["n"]})
    return out


def _merge(ss: pd.DataFrame, ld: LDRef) -> pd.DataFrame:
    # reference-panel order, as ldsc.py does: jackknife blocks are contiguous in this order
    m = ld.ref_ld[["SNP", "L2"]].merge(ss, on="SNP").merge(ld.w_ld, on="SNP")
    return m.reset_index(drop=True)


def _hsq_dict(h: reg.Hsq, n_snp: int, P: float | None, K: float | None) -> dict[str, Any]:
    tot, se = float(np.ravel(h.tot)[0]), float(np.ravel(h.tot_se)[0])
    out: dict[str, Any] = {
        "h2_observed": tot, "h2_observed_se": se,
        "h2_z": tot / se if se > 0 else None,
        "intercept": float(np.ravel(h.intercept)[0]),
        "intercept_se": float(np.ravel(h.intercept_se)[0]) if h.intercept_se != "NA" else None,
        "lambda_gc": float(h.lambda_gc), "mean_chi2": float(h.mean_chisq),
        "n_snps": int(n_snp),
    }
    if h.mean_chisq > 1:
        out["ratio"] = float(np.ravel(h.ratio)[0])
        out["ratio_se"] = float(np.ravel(h.ratio_se)[0])
    if P is not None and K is not None:
        c = reg.h2_obs_to_liab(1, P, K)
        out["h2_liability"] = tot * c
        out["h2_liability_se"] = se * c
        out["sample_prevalence"] = P
        out["population_prevalence"] = K
    return out


def _interpret_h2(res: dict[str, Any], label: str, w: Warnings) -> None:
    if res["mean_chi2"] < 1.02:
        w.add("ldsc_low_chi2", f"{label}: mean chi2 = {res['mean_chi2']:.3f} < 1.02. LDSC is "
              "unreliable at this signal level (too little polygenic signal or N too small).", "serious")
    if res.get("h2_z") is not None and res["h2_z"] < 4:
        w.add("ldsc_h2_weak", f"{label}: h2 z-score {res['h2_z']:.1f} < 4. Genetic correlations "
              "with this trait will be imprecise and hard to interpret.", "warning")
    if res["intercept"] > 1.1:
        w.add("ldsc_intercept_high", f"{label}: LDSC intercept {res['intercept']:.3f}. Suggests "
              "confounding (stratification, cryptic relatedness) or misspecified N.", "serious")
    if res["intercept"] < 0.95:
        w.add("ldsc_intercept_low", f"{label}: intercept {res['intercept']:.3f} < 1: sumstats may "
              "be GC-corrected already.", "warning")


def _two_step(n_annot: int, intercept: float | None) -> int | None:
    return 30 if n_annot == 1 and intercept is None else None


def h2(df: pd.DataFrame, meta: TraitMeta, ld: LDRef, w: Warnings, *, n_blocks: int = 200,
       intercept: float | None = None, pop_prev: float | None = None) -> dict[str, Any]:
    check_ancestry(meta.ancestry, ld.ancestry, "the LD score reference")
    ss = _merge(munge(df, meta, ld, w, meta.trait), ld)
    if len(ss) < 200_000:
        w.add("ldsc_few_snps", f"{meta.trait}: only {len(ss)} SNPs overlap the LD scores (LDSC "
              "expects ~1M HapMap3 SNPs); standard errors will be larger.", "warning")
    n = len(ss)
    s = lambda x: np.asarray(x, dtype=float).reshape((n, 1))  # noqa: E731
    h = reg.Hsq(s(ss["Z"] ** 2), s(ss["L2"]), s(ss["w"]), s(ss["N"]), np.array([[ld.M]]),
                n_blocks=min(n, n_blocks), intercept=intercept, twostep=_two_step(1, intercept))
    P = meta.sample_prevalence
    K = pop_prev if pop_prev is not None else meta.population_prevalence
    res = _hsq_dict(h, n, P, K)
    _interpret_h2(res, meta.trait, w)
    return res


def _align(m: pd.DataFrame) -> pd.DataFrame:
    """Align trait 2 Z to trait 1 alleles (LDSC FLIP_ALLELES logic), drop incompatible."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    a1, a2, b1, b2 = m["A1"], m["A2"], m["A1x"], m["A2x"]
    same = (a1 == b1) & (a2 == b2)
    flip = (a1 == b2) & (a2 == b1)
    cs = (a1 == b1.map(comp)) & (a2 == b2.map(comp))
    cf = (a1 == b2.map(comp)) & (a2 == b1.map(comp))
    ok = same | flip | cs | cf
    m = m[ok].copy()
    sign = np.where((flip | cf)[ok], -1.0, 1.0)
    m["Z2"] = m["Z2"] * sign
    return m


def rg(d1: pd.DataFrame, m1: TraitMeta, d2: pd.DataFrame, m2: TraitMeta, ld: LDRef, w: Warnings,
       *, n_blocks: int = 200, intercept_gencov: float | None = None,
       pop_prev: tuple[float | None, float | None] = (None, None)) -> dict[str, Any]:
    for m in (m1, m2):
        check_ancestry(m.ancestry, ld.ancestry, "the LD score reference")
    if m1.ancestry != m2.ancestry:
        w.add("rg_cross_ancestry", "The two GWAS have different ancestries; standard LDSC rg "
              "assumes shared LD. Use a cross-ancestry method (e.g. Popcorn) instead.", "serious")
    s1 = _merge(munge(d1, m1, ld, w, m1.trait), ld)
    s2 = munge(d2, m2, ld, w, m2.trait).rename(columns={"A1": "A1x", "A2": "A2x", "Z": "Z2", "N": "N2"})
    m = s1.rename(columns={"Z": "Z1", "N": "N1"}).merge(s2, on="SNP")
    m = _align(m).reset_index(drop=True)
    n = len(m)
    if n < 1000:
        raise GenorError(f"Only {n} SNPs shared between the traits after alignment.",
                         code="ldsc_too_few_shared")
    s = lambda x: np.asarray(x, dtype=float).reshape((n, 1))  # noqa: E731
    r = reg.RG(s(m["Z1"]), s(m["Z2"]), s(m["L2"]), s(m["w"]), s(m["N1"]), s(m["N2"]),
               np.array([[ld.M]]), intercept_gencov=intercept_gencov,
               n_blocks=min(n, n_blocks), twostep=None)  # ldsc.py --rg never enables two-step (intercept_h2 is a list there)
    P = (m1.sample_prevalence, m2.sample_prevalence)
    h1 = _hsq_dict(r.hsq1, n, P[0], pop_prev[0] or m1.population_prevalence)
    h2_ = _hsq_dict(r.hsq2, n, P[1], pop_prev[1] or m2.population_prevalence)
    _interpret_h2(h1, m1.trait, w)
    _interpret_h2(h2_, m2.trait, w)
    gc = r.gencov
    out: dict[str, Any] = {
        "trait1": m1.trait, "trait2": m2.trait, "n_snps": n,
        "gencov": float(np.ravel(gc.tot)[0]), "gencov_se": float(np.ravel(gc.tot_se)[0]),
        "gencov_intercept": float(np.ravel(gc.intercept)[0]),
        "gencov_intercept_se": (float(np.ravel(gc.intercept_se)[0])
                                if not isinstance(gc.intercept_se, str) else None),
        "h2_trait1": h1, "h2_trait2": h2_,
    }
    if r._negative_hsq:
        out.update(rg=None, rg_se=None, p=None)
        w.add("rg_negative_h2", "One trait has h2 <= 0; rg is undefined.", "serious")
    else:
        out.update(rg=r.rg_ratio, rg_se=r.rg_se, z=float(r.z), p=float(r.p))
        if abs(r.rg_ratio) > 1.2:
            w.add("rg_out_of_bounds", f"rg = {r.rg_ratio:.2f} is outside [-1, 1]; usually low "
                  "h2 in one trait. Do not interpret.", "serious")
    gi, gse = out["gencov_intercept"], out["gencov_intercept_se"]
    if gse and abs(gi / gse) > 3:
        w.add("sample_overlap_signal", f"Cross-trait intercept {gi:.3f} (SE {gse:.3f}) differs "
              "from 0: evidence of sample overlap or shared confounding. Relevant for MR.",
              "warning", gencov_intercept=gi)
    out["gencov_intercept_p"] = float(2 * stats.norm.sf(abs(gi / gse))) if gse else None
    return out
