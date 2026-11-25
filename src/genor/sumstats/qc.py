# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Summary-statistics QC.

Every filter is counted; nothing is silently dropped. Defaults and why:

- MAF >= 0.01 (when frequency is present): rare-variant estimates are noisy and LD
  reference panels of ~500 people cannot tag them.
- INFO >= 0.8 (when present): poorly imputed variants inflate noise.
- P in (0, 1], SE > 0, finite beta: malformed rows.
- Duplicate variant IDs / positions: keep the row with the smallest P.
- Non-SNV alleles are kept but counted; LDSC and clumping use SNVs only.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from genor.errors import Warnings
from genor.guards import MHC

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def is_palindromic(a1: pd.Series, a2: pd.Series) -> pd.Series:
    return a1.map(COMPLEMENT).eq(a2)


def is_snv(a1: pd.Series, a2: pd.Series) -> pd.Series:
    return a1.isin(list("ACGT")) & a2.isin(list("ACGT"))


def lambda_gc(p: pd.Series) -> float:
    p = p.dropna()
    p = p[(p > 0) & (p <= 1)]
    if len(p) == 0:
        return float("nan")
    chi = stats.chi2.isf(np.median(p), 1)
    return float(chi / stats.chi2.ppf(0.5, 1))


def mhc_mask(df: pd.DataFrame, build: str | None) -> pd.Series:
    if "chromosome" not in df or "base_pair_location" not in df:
        return pd.Series(False, index=df.index)
    c, start, end = MHC.get(build or "GRCh37", MHC["GRCh37"])
    return (df["chromosome"].astype(str) == str(c)) & df["base_pair_location"].between(start, end)


def run_qc(df: pd.DataFrame, build: str | None, *, maf: float = 0.01, info: float = 0.8,
           drop_mhc: bool = False, drop_palindromic_ambiguous: bool = False,
           w: Warnings | None = None) -> tuple[pd.DataFrame, dict[str, Any], Warnings]:
    w = w if w is not None else Warnings()
    report: dict[str, Any] = {"n_input": int(len(df)), "removed": {}, "flags": {}}
    keep = pd.Series(True, index=df.index)

    def drop(mask: pd.Series, reason: str) -> None:
        nonlocal keep
        m = mask & keep
        report["removed"][reason] = int(m.sum())
        keep &= ~mask

    drop(df["p_value"].isna() | ~df["p_value"].between(0, 1, inclusive="right"), "bad_p")
    if "beta" in df:
        drop(~np.isfinite(df["beta"]), "bad_beta")
    if "standard_error" in df:
        drop(~(df["standard_error"] > 0), "bad_se")
    drop(df["effect_allele"].isna() | df["other_allele"].isna()
         | (df["effect_allele"] == df["other_allele"]), "bad_alleles")
    if "effect_allele_frequency" in df and df["effect_allele_frequency"].notna().any():
        f = df["effect_allele_frequency"]
        drop(~f.between(0, 1), "bad_eaf")
        drop(np.minimum(f, 1 - f) < maf, f"maf_lt_{maf}")
    else:
        w.add("no_eaf", "No allele frequency column: MAF filter skipped and palindromic SNPs "
              "cannot be strand-resolved (they will be dropped in harmonisation).", "warning")
    if "info" in df and df["info"].notna().any():
        drop(df["info"] < info, f"info_lt_{info}")
    if "chromosome" in df:
        drop(~df["chromosome"].astype(str).isin([str(i) for i in range(1, 27)]), "bad_chromosome")

    d = df[keep].copy()
    # duplicates: keep the smallest P
    d = d.sort_values("p_value")
    before = len(d)
    if "variant_id" in d and d["variant_id"].notna().any():
        d = d[~(d["variant_id"].notna() & d["variant_id"].duplicated())]
    if "chromosome" in d and "base_pair_location" in d:
        a1, a2 = d["effect_allele"].astype(str), d["other_allele"].astype(str)
        lo = a1.where(a1 < a2, a2)
        hi = a2.where(a1 < a2, a1)
        key = (d["chromosome"].astype(str) + ":" + d["base_pair_location"].astype("Int64").astype(str)
               + ":" + lo + ":" + hi)
        d = d[~key.duplicated()]
    report["removed"]["duplicates"] = int(before - len(d))

    snv = is_snv(d["effect_allele"], d["other_allele"])
    pal = is_palindromic(d["effect_allele"], d["other_allele"])
    report["flags"]["non_snv"] = int((~snv).sum())
    report["flags"]["palindromic"] = int(pal.sum())
    if "effect_allele_frequency" in d:
        f = d["effect_allele_frequency"]
        amb = pal & f.between(0.4, 0.6)
        report["flags"]["palindromic_ambiguous_freq"] = int(amb.sum())
        if drop_palindromic_ambiguous:
            d = d[~amb]
            report["removed"]["palindromic_ambiguous"] = int(amb.sum())
    mhc = mhc_mask(d, build)
    report["flags"]["mhc"] = int(mhc.sum())
    if drop_mhc:
        d = d[~mhc]
        report["removed"]["mhc"] = int(mhc.sum())

    d = d.sort_values(["chromosome", "base_pair_location"]) if "chromosome" in d else d
    report["n_output"] = int(len(d))
    lam = lambda_gc(d["p_value"])
    chi2 = (d["z"] ** 2).mean() if "z" in d else float("nan")
    report["lambda_gc"] = round(lam, 4)
    report["mean_chi2"] = round(float(chi2), 4)
    report["n_genome_wide_sig"] = int((d["p_value"] < 5e-8).sum())
    if "n" in d and d["n"].notna().any():
        report["n_median"] = float(d["n"].median())
        if "n_cases" in d and "n_controls" in d and d["n_cases"].notna().any() and d["n_controls"].notna().any():
            # case-control form: 1000 cases and 1000 controls
            nca, nco = float(d["n_cases"].median()), float(d["n_controls"].median())
            report["lambda_1000"] = round(1 + (lam - 1) * (1 / nca + 1 / nco) * 500, 4)
        else:
            report["lambda_1000"] = round(1 + (lam - 1) * 1000 / float(d["n"].median()), 4)
    report["min_p"] = float(d["p_value"].min()) if len(d) else None

    removed = report["n_input"] - report["n_output"]
    if report["n_input"] and removed / report["n_input"] > 0.3:
        w.add("qc_removed_many", f"QC removed {removed / report['n_input']:.0%} of variants; check "
              "column mapping and frequency/INFO scales.", "serious", removed=removed)
    if lam > 1.3:
        w.add("inflation_high", f"lambda_GC = {lam:.2f}. High for most GWAS; can reflect "
              "polygenicity in large samples. LDSC intercept separates the two.", "warning")
    if lam < 0.9:
        w.add("deflation", f"lambda_GC = {lam:.2f} (< 0.9). The file may already be GC-corrected "
              "or filtered to significant hits.", "serious")
    if report["n_output"] < 50_000:
        w.add("few_variants", f"Only {report['n_output']} variants after QC; genome-wide methods "
              "(LDSC) need a genome-wide file, not a hit list.", "warning")
    if report["flags"]["mhc"]:
        w.add("mhc_present", f"{report['flags']['mhc']} variants in the MHC (chr6:25-34 Mb). "
              "LDSC weights exclude it; consider excluding it from MR instruments.", "info")
    return d.reset_index(drop=True), report, w
