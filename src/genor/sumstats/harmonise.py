# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Allele harmonisation between two sets of summary statistics.

Follows TwoSampleMR's ``action = 2`` semantics: align the second dataset to the
effect allele of the first, fix strand flips, and resolve palindromic SNPs with
allele frequencies (dropped when MAF > 0.42 or frequencies are missing).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from genor.errors import GenorError, Warnings
from genor.sumstats.qc import COMPLEMENT

PALINDROME_MAF = 0.42


def _comp(s: pd.Series) -> pd.Series:
    return s.map(COMPLEMENT)


def join_key(a: pd.DataFrame, b: pd.DataFrame) -> str:
    if "variant_id" in a and "variant_id" in b:
        ia = a["variant_id"].dropna().astype(str)
        if ia.str.startswith("rs").mean() > 0.5 and len(set(ia) & set(b["variant_id"].dropna())) > 0:
            return "variant_id"
    if all(c in x for x in (a, b) for c in ("chromosome", "base_pair_location")):
        return "chrpos"
    raise GenorError("Cannot match variants: no shared rsIDs and no chromosome/position in both.",
                     code="cannot_match_variants")


def _with_key(d: pd.DataFrame, key: str) -> pd.DataFrame:
    d = d.copy()
    if key == "chrpos":
        d["_key"] = d["chromosome"].astype(str) + ":" + d["base_pair_location"].astype("Int64").astype(str)
    else:
        d["_key"] = d["variant_id"].astype(str)
    return d.drop_duplicates("_key")


def harmonise(a: pd.DataFrame, b: pd.DataFrame, suffixes: tuple[str, str] = ("_x", "_y"),
              w: Warnings | None = None) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return a merged table with dataset ``b`` aligned to ``a``'s effect allele."""
    w = w if w is not None else Warnings()
    key = join_key(a, b)
    sa, sb = suffixes
    m = _with_key(a, key).merge(_with_key(b, key), on="_key", suffixes=(sa, sb))
    ea, oa = m[f"effect_allele{sa}"], m[f"other_allele{sa}"]
    eb, ob = m[f"effect_allele{sb}"], m[f"other_allele{sb}"]

    same = (ea == eb) & (oa == ob)
    swap = (ea == ob) & (oa == eb)
    strand = (ea == _comp(eb)) & (oa == _comp(ob)) & ~same & ~swap
    strand_swap = (ea == _comp(ob)) & (oa == _comp(eb)) & ~same & ~swap & ~strand
    pal = _comp(ea).eq(oa)

    flip = swap | strand_swap
    fb = f"effect_allele_frequency{sb}"
    fa = f"effect_allele_frequency{sa}"
    for c in (f"beta{sb}", f"z{sb}"):
        if c in m:
            m.loc[flip, c] = -m.loc[flip, c]
    if fb in m:
        m.loc[flip, fb] = 1 - m.loc[flip, fb]

    # palindromic: alleles always look "same" or "swap"; use frequencies to decide
    drop_pal = pd.Series(False, index=m.index)
    if pal.any():
        if fa in m and fb in m:
            fa_, fb_ = m[fa], m[fb]
            maf_a = np.minimum(fa_, 1 - fa_)
            ambiguous = pal & ((maf_a > PALINDROME_MAF) | fa_.isna() | fb_.isna())
            disagree = pal & ~ambiguous & ((fa_ < 0.5) != (fb_ < 0.5))
            for c in (f"beta{sb}", f"z{sb}"):
                if c in m:
                    m.loc[disagree, c] = -m.loc[disagree, c]
            m.loc[disagree, fb] = 1 - m.loc[disagree, fb]
            drop_pal = ambiguous
        else:
            drop_pal = pal
            w.add("palindromic_dropped_no_freq", f"{int(pal.sum())} palindromic SNPs dropped "
                  "(no allele frequencies to resolve strand).", "warning")

    ok = (same | swap | strand | strand_swap) & ~drop_pal
    report = {
        "matched_on": key,
        "n_a": int(len(a)), "n_b": int(len(b)), "n_shared": int(len(m)),
        "aligned": int(same.sum()), "swapped": int(swap.sum()),
        "strand_flipped": int((strand | strand_swap).sum()),
        "palindromic": int(pal.sum()), "palindromic_dropped": int(drop_pal.sum()),
        "allele_mismatch_dropped": int((~(same | swap | strand | strand_swap)).sum()),
    }
    report["n_kept"] = int(ok.sum())
    if report["n_shared"] and report["allele_mismatch_dropped"] / report["n_shared"] > 0.05:
        w.add("allele_mismatch", f"{report['allele_mismatch_dropped']} of {report['n_shared']} "
              "shared variants have incompatible alleles: possible build or ID mismatch.", "serious")
    out = m[ok].copy()
    # canonical alleles = dataset a's
    out["effect_allele"] = out[f"effect_allele{sa}"]
    out["other_allele"] = out[f"other_allele{sa}"]
    return out.drop(columns=["_key"]), report
