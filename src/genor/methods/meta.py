# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Variant-level meta-analysis of K GWAS (native replacement for METAL, which cannot be shipped).

Schemes:
- ``ivw``: fixed-effects inverse-variance weighting (METAL SCHEME STDERR), plus
  DerSimonian-Laird random effects, Cochran's Q and I^2 (METAL HETEROGENEITY ON).
- ``samplesize``: weighted Z (METAL SCHEME SAMPLESIZE), weights sqrt(N).

All studies are aligned to the first study's alleles before combining.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from genor.errors import GenorError, Warnings
from genor.sumstats.load import TraitMeta

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _key(d: pd.DataFrame) -> pd.Series:
    if "variant_id" in d and d["variant_id"].astype(str).str.startswith("rs").mean() > 0.5:
        return d["variant_id"].astype(str)
    return d["chromosome"].astype(str) + ":" + d["base_pair_location"].astype("Int64").astype(str)


def _align_all(tables: list[pd.DataFrame], w: Warnings) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Stack all studies, aligned to the alleles of the first study that reports each variant.

    Like METAL: matching alleles keep their sign, swapped alleles flip it, strand-complemented
    alleles are accepted, anything else is dropped and counted. Palindromic SNPs are matched
    as written (METAL does the same), so files must be on the same strand.
    """
    cols = ["beta", "standard_error", "z", "n", "effect_allele_frequency"]
    parts = []
    for i, t in enumerate(tables):
        d = t.reindex(columns=["effect_allele", "other_allele", *cols]).copy()
        d["_k"] = _key(t).to_numpy()
        d["_study"] = i
        parts.append(d.drop_duplicates("_k"))
    lt = pd.concat(parts, ignore_index=True)
    anchor = lt.drop_duplicates("_k").set_index("_k")[["effect_allele", "other_allele"]]
    a1 = anchor["effect_allele"].reindex(lt["_k"]).to_numpy()
    a2 = anchor["other_allele"].reindex(lt["_k"]).to_numpy()
    e, o = lt["effect_allele"].to_numpy(), lt["other_allele"].to_numpy()
    ce = pd.Series(e).map(COMP).to_numpy()
    co = pd.Series(o).map(COMP).to_numpy()
    # direct matches win over strand complements (palindromes match both ways)
    same_d = (e == a1) & (o == a2)
    swap_d = (e == a2) & (o == a1)
    same = same_d | ((ce == a1) & (co == a2) & ~swap_d)
    swap = (swap_d | ((ce == a2) & (co == a1) & ~same_d)) & ~same
    for c in ("beta", "z"):
        lt.loc[swap, c] = -lt.loc[swap, c]
    lt.loc[swap, "effect_allele_frequency"] = 1 - lt.loc[swap, "effect_allele_frequency"]
    bad = ~(same | swap)
    reports = []
    for i in range(len(tables)):
        m = lt["_study"].to_numpy() == i
        reports.append({"study": i, "n": int(m.sum()), "flipped": int((swap & m).sum()),
                        "allele_mismatch_dropped": int((bad & m).sum())})
    if bad.sum():
        w.add("meta_allele_mismatch", f"{int(bad.sum())} study-variant rows had alleles "
              "incompatible with the first study and were dropped.", "warning")
    lt = lt[~bad].copy()
    lt["effect_allele"] = a1[~bad]
    lt["other_allele"] = a2[~bad]
    return lt, reports


def meta_analyse(tables: list[pd.DataFrame], metas: list[TraitMeta], scheme: str = "ivw",
                 w: Warnings | None = None) -> tuple[pd.DataFrame, dict[str, Any]]:
    w = w if w is not None else Warnings()
    if len(tables) < 2:
        raise GenorError("Meta-analysis needs at least two studies.", code="meta_too_few")
    builds = {m.build for m in metas}
    if len(builds) > 1:
        raise GenorError(f"Studies are on different builds: {sorted(map(str, builds))}.",
                         code="build_mismatch", hint="Lift all studies to one build first.")
    ancs = {m.ancestry for m in metas}
    if len(ancs) > 1:
        w.add("cross_ancestry_meta", f"Studies span ancestries {sorted(map(str, ancs))}. Fixed-effect "
              "IVW assumes one true effect; inspect I^2 and random-effects results. Trans-ancestry "
              "methods (MR-MEGA) model ancestry-related heterogeneity better.", "warning")
    lt, reports = _align_all(tables, w)
    g = lt.groupby("_k", sort=False)
    k = g["_study"].nunique()
    out = pd.DataFrame({"variant_id": k.index, "n_studies": k.to_numpy()})
    first = lt.drop_duplicates("_k").set_index("_k")
    out["effect_allele"] = first.loc[out["variant_id"], "effect_allele"].to_numpy()
    out["other_allele"] = first.loc[out["variant_id"], "other_allele"].to_numpy()
    if scheme == "ivw":
        lt["w"] = 1 / lt["standard_error"] ** 2
        lt["wb"] = lt["w"] * lt["beta"]
        lt["wb2"] = lt["w"] * lt["beta"] ** 2
        lt["w2"] = lt["w"] ** 2
        s = lt.groupby("_k", sort=False)[["w", "wb", "wb2", "w2"]].sum()
        beta = s["wb"] / s["w"]
        se = np.sqrt(1 / s["w"])
        Q = (s["wb2"] - s["wb"] ** 2 / s["w"]).clip(lower=0)
        df_q = (k - 1).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            i2 = np.where(Q > 0, np.clip((Q - df_q) / Q, 0, 1), 0.0)
            tau2 = np.clip((Q - df_q) / (s["w"] - s["w2"] / s["w"]), 0, None).fillna(0)
        lt = lt.join(tau2.rename("tau2"), on="_k")
        lt["wr"] = 1 / (lt["standard_error"] ** 2 + lt["tau2"])
        lt["wrb"] = lt["wr"] * lt["beta"]
        r = lt.groupby("_k", sort=False)[["wr", "wrb"]].sum()
        out["beta"] = beta.to_numpy()
        out["standard_error"] = se.to_numpy()
        out["z"] = out["beta"] / out["standard_error"]
        out["p_value"] = 2 * stats.norm.sf(np.abs(out["z"]))
        out["het_Q"] = Q.to_numpy()
        out["het_df"] = df_q
        out["het_p"] = np.where(df_q > 0, stats.chi2.sf(Q.to_numpy(), np.maximum(df_q, 1)), np.nan)
        out["het_I2"] = i2
        out["tau2"] = tau2.to_numpy()
        out["beta_random"] = (r["wrb"] / r["wr"]).to_numpy()
        out["se_random"] = np.sqrt(1 / r["wr"]).to_numpy()
        out["p_random"] = 2 * stats.norm.sf(np.abs(out["beta_random"] / out["se_random"]))
    elif scheme == "samplesize":
        if lt["n"].isna().any():
            raise GenorError("Sample-size scheme needs N for every study.", code="n_missing")
        lt["wz"] = np.sqrt(lt["n"]) * lt["z"]
        s = lt.groupby("_k", sort=False)[["wz", "n"]].sum()
        z = s["wz"] / np.sqrt(s["n"])
        out["z"] = z.to_numpy()
        out["p_value"] = 2 * stats.norm.sf(np.abs(out["z"]))
    else:
        raise GenorError(f"Unknown scheme '{scheme}'", code="bad_param", hint="ivw | samplesize")
    out["n"] = lt.groupby("_k", sort=False)["n"].sum(min_count=1).to_numpy()
    lt["fw"] = lt["effect_allele_frequency"] * lt["n"]
    fsum = lt.groupby("_k", sort=False)[["fw", "n"]].sum(min_count=1)
    out["effect_allele_frequency"] = (fsum["fw"] / fsum["n"]).to_numpy()
    # carry positions from whichever study had them
    with_pos = [t.assign(_k=_key(t).to_numpy()) for t in tables if "chromosome" in t]
    if with_pos:
        pos = pd.concat(with_pos)
        pos = pos.drop_duplicates("_k").set_index("_k")
        out["chromosome"] = pos["chromosome"].reindex(out["variant_id"]).to_numpy()
        out["base_pair_location"] = pos["base_pair_location"].reindex(out["variant_id"]).to_numpy()
    summary = {"scheme": scheme, "n_studies": len(tables), "n_variants": int(len(out)),
               "n_in_all_studies": int((out["n_studies"] == len(tables)).sum()),
               "n_genome_wide_sig": int((out["p_value"] < 5e-8).sum()), "alignment": reports}
    if scheme == "ivw":
        sig = out[out["p_value"] < 5e-8]
        summary["n_sig_with_high_I2"] = int((sig["het_I2"] > 0.5).sum())
        if summary["n_sig_with_high_I2"]:
            w.add("heterogeneous_hits", f"{summary['n_sig_with_high_I2']} genome-wide significant "
                  "variants have I^2 > 50%: check random-effects results.", "warning")
    return out, summary
