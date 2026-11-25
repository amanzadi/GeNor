# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Wakefield ABF fine-mapping (single causal variant) and coloc ABF colocalisation.

Ports of ``coloc::finemap.abf`` and ``coloc::coloc.abf`` (coloc 5.x) for summary
statistics with beta/varbeta (preferred) or p-values + MAF + N. Default priors:
p1 = p2 = 1e-4, p12 = 1e-5 (Giambartolomei 2014; Wallace 2020 discusses p12 = 5e-6).
Prior SD of effect: 0.15 x sdY for quantitative traits, 0.2 for case-control (log OR).

Both assume at most one causal variant per trait in the region. Multiple signals
break that assumption: check with COJO or SuSiE before trusting H3 vs H4.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import logsumexp

from genor.errors import GenorError, Warnings
from genor.sumstats.load import TraitMeta


def parse_region(region: str) -> tuple[str, int, int]:
    try:
        c, rng = region.replace("chr", "").split(":")
        a, b = rng.replace(",", "").split("-")
        return c, int(a), int(b)
    except ValueError as e:
        raise GenorError(f"Bad region '{region}'", code="bad_region", hint="Use chr:start-end") from e


def subset_region(df: pd.DataFrame, region: str) -> pd.DataFrame:
    c, a, b = parse_region(region)
    m = (df["chromosome"].astype(str) == c) & df["base_pair_location"].between(a, b)
    return df[m].copy()


def sdy_est(varbeta: np.ndarray, maf: np.ndarray, n: np.ndarray) -> float:
    """coloc:::sdY.est: regress 2*N*MAF*(1-MAF) on 1/varbeta through the origin."""
    oneover = 1 / varbeta
    nvx = 2 * n * maf * (1 - maf)
    cf = np.sum(oneover * nvx) / np.sum(oneover**2)
    return float(np.sqrt(cf))


def labf(d: pd.DataFrame, meta: TraitMeta, w: Warnings, sdy: float | None = None) -> tuple[np.ndarray, dict[str, Any]]:
    binary = meta.trait_type == "binary"
    info: dict[str, Any] = {"type": "cc" if binary else "quant"}
    if "beta" in d and "standard_error" in d and d["beta"].notna().all():
        vb = d["standard_error"].to_numpy(float) ** 2
        z = d["beta"].to_numpy(float) / np.sqrt(vb)
        if binary:
            sd_prior = 0.2
        else:
            if sdy is None:
                if "effect_allele_frequency" in d and "n" in d and d["n"].notna().all():
                    f = d["effect_allele_frequency"].to_numpy(float)
                    maf = np.minimum(f, 1 - f)
                    sdy = sdy_est(vb, maf, d["n"].to_numpy(float))
                    info["sdY_estimated"] = sdy
                else:
                    sdy = 1.0
                    w.add("sdy_assumed", f"{meta.trait}: sdY unknown and not estimable (needs "
                          "MAF and N); assuming 1 (standardised trait).", "warning")
            sd_prior = 0.15 * sdy
        r = sd_prior**2 / (sd_prior**2 + vb)
        info["sd_prior"] = sd_prior
        return 0.5 * (np.log(1 - r) + r * z**2), info
    # p-value route
    if "effect_allele_frequency" not in d or "n" not in d:
        raise GenorError("coloc needs beta+SE, or P with MAF and N.", code="coloc_inputs")
    f = d["effect_allele_frequency"].to_numpy(float)
    maf = np.minimum(f, 1 - f)
    n = d["n"].to_numpy(float)
    z = stats.norm.isf(d["p_value"].to_numpy(float) / 2)
    if binary:
        s = meta.sample_prevalence
        if s is None:
            raise GenorError("Case-control p-value route needs the case fraction.", code="coloc_inputs",
                             hint="Reload with --n-cases/--n-controls.")
        V = 1 / (2 * n * maf * (1 - maf) * s * (1 - s))
        sd_prior = 0.2
    else:
        V = 1 / (2 * n * maf * (1 - maf))
        sd_prior = 0.15
        w.add("sdy_assumed", f"{meta.trait}: p-value route assumes a standardised trait (sdY = 1).",
              "warning")
    r = sd_prior**2 / (sd_prior**2 + V)
    return 0.5 * (np.log(1 - r) + r * z**2), info


def credible_set(pp: np.ndarray, ids: np.ndarray, level: float = 0.95) -> list[str]:
    o = np.argsort(-pp)
    cs = np.cumsum(pp[o])
    k = int(np.searchsorted(cs, level) + 1)
    return [str(x) for x in ids[o][:k]]


def finemap_abf(d: pd.DataFrame, meta: TraitMeta, w: Warnings, p1: float = 1e-4,
                level: float = 0.95) -> tuple[pd.DataFrame, dict[str, Any]]:
    if len(d) < 2:
        raise GenorError("Region has fewer than 2 variants.", code="region_empty")
    lab, info = labf(d, meta, w)
    nsnp = len(lab)
    if nsnp * p1 >= 1:
        raise GenorError(f"{nsnp} variants x prior {p1:g} >= 1: the null prior would be negative.",
                         code="bad_param", hint="Use a smaller region or a smaller --p1 (e.g. 1/(10*nsnp)).")
    prior = np.full(nsnp, np.log(p1))
    null = np.log(1 - nsnp * p1)
    allv = np.append(lab + prior, 0 + null)
    pp = np.exp(allv - logsumexp(allv))
    res = d[["variant_id", "chromosome", "base_pair_location", "p_value"]].copy()
    res["lABF"] = lab
    res["PP"] = pp[:-1]
    res = res.sort_values("PP", ascending=False)
    ids = res["variant_id"].to_numpy()
    # credible set conditional on one causal variant being present (PPs renormalised without the null)
    pp_c = res["PP"].to_numpy() / res["PP"].sum()
    cs = credible_set(pp_c, ids, level)
    summary = {"n_variants": nsnp, "pp_null": float(pp[-1]), "credible_set_level": level,
               "credible_set_size": len(cs), "credible_set": cs[:50],
               "top_variant": str(ids[0]), "top_pp": float(res["PP"].iloc[0]), **info}
    if summary["pp_null"] > 0.5:
        w.add("finemap_no_signal", "Posterior on 'no causal variant' > 0.5: no clear signal here.",
              "warning")
    return res, summary


def coloc_abf(d1: pd.DataFrame, m1: TraitMeta, d2: pd.DataFrame, m2: TraitMeta, w: Warnings,
              p1: float = 1e-4, p2: float = 1e-4, p12: float = 1e-5) -> tuple[pd.DataFrame, dict[str, Any]]:
    key = "variant_id"
    m = d1.merge(d2, on=key, suffixes=("_1", "_2"))
    if len(m) < 10:
        raise GenorError(f"Only {len(m)} shared variants in the region.", code="coloc_too_few",
                         hint="Check that both traits use the same build/IDs and cover the region.")
    def side(sfx: str) -> pd.DataFrame:
        cols = {c: c[:-2] for c in m.columns if c.endswith(sfx)}
        return m[[*cols]].rename(columns=cols)
    l1, i1 = labf(side("_1"), m1, w)
    l2, i2 = labf(side("_2"), m2, w)
    lsum = l1 + l2
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(lsum)
    lh3 = np.log(p1) + np.log(p2) + (s1 + s2) + np.log1p(-np.exp(s12 - (s1 + s2))) \
        if s1 + s2 > s12 else -np.inf
    labs = np.array([0.0, np.log(p1) + s1, np.log(p2) + s2, lh3, np.log(p12) + s12])
    pp = np.exp(labs - logsumexp(labs))
    snp_pp = np.exp(lsum - s12)
    res = pd.DataFrame({"variant_id": m[key], "SNP.PP.H4": snp_pp}).sort_values("SNP.PP.H4", ascending=False)
    summary = {"nsnps": int(len(m)), "PP.H0": pp[0], "PP.H1": pp[1], "PP.H2": pp[2],
               "PP.H3": pp[3], "PP.H4": pp[4], "priors": {"p1": p1, "p2": p2, "p12": p12},
               "top_shared_variant": str(res["variant_id"].iloc[0]),
               "prior_trait1": i1, "prior_trait2": i2}
    if pp[4] > 0.8:
        summary["interpretation"] = "Strong support for a shared causal variant (PP.H4 > 0.8)."
    elif pp[3] > 0.8:
        summary["interpretation"] = "Strong support for distinct causal variants (PP.H3 > 0.8)."
    elif pp[0] + pp[1] + pp[2] > 0.5:
        summary["interpretation"] = "At least one trait lacks a clear signal in this region (H0-H2)."
    else:
        summary["interpretation"] = "Inconclusive."
    return res, {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in summary.items()}
