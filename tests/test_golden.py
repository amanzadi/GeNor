# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Native ports vs the reference R implementations.

Reference values in tests/golden/expected_r.json come from scripts/golden/make_golden.R
(TwoSampleMR 0.7.14, coloc 5.2.3, meta 9.0.0). Deterministic quantities must agree to
numerical precision; bootstrap standard errors use different RNGs and are compared
loosely.
"""

import numpy as np
import pandas as pd
import pytest

from genor.errors import Warnings
from genor.methods import mr
from genor.methods.coloc import coloc_abf, finemap_abf
from genor.methods.meta import meta_analyse
from genor.sumstats.load import TraitMeta

from .conftest import GOLDEN

TIGHT = dict(rel=1e-8, abs=1e-12)


@pytest.fixture(scope="module")
def harm():
    h = pd.read_csv(GOLDEN / "mr_harmonised.tsv", sep="\t")
    return (h["beta_exp"].to_numpy(), h["beta_out"].to_numpy(), h["se_exp"].to_numpy(),
            h["se_out"].to_numpy())


def test_ivw(harm, expected_r):
    e = expected_r["mr"]["ivw"]
    r = mr.ivw(*harm)
    for k in ("b", "se", "pval", "Q", "Q_pval"):
        assert r[k] == pytest.approx(e[k], **TIGHT), k
    assert r["Q_df"] == e["Q_df"]
    fe = mr.ivw(*harm, fixed=True)
    assert fe["se"] == pytest.approx(expected_r["mr"]["ivw_fe"]["se"], **TIGHT)


def test_egger(harm, expected_r):
    e = expected_r["mr"]["egger"]
    r = mr.egger(*harm)
    pairs = {"b": "b", "se": "se", "pval": "pval", "intercept": "b_i", "intercept_se": "se_i",
             "intercept_pval": "pval_i", "Q": "Q", "Q_pval": "Q_pval"}
    for ours, theirs in pairs.items():
        assert r[ours] == pytest.approx(e[theirs], **TIGHT), ours


def test_wald_and_loo(harm, expected_r):
    bx, by, sx, sy = harm
    w = mr.wald_ratio(bx[0], by[0], sx[0], sy[0])
    for k in ("b", "se", "pval"):
        assert w[k] == pytest.approx(expected_r["mr"]["wald_first"][k], **TIGHT)
    loo = mr.ivw(bx[1:], by[1:], sx[1:], sy[1:])
    assert loo["b"] == pytest.approx(expected_r["mr"]["loo_drop_first"]["b"], **TIGHT)


def test_weighted_median(harm, expected_r):
    e = expected_r["mr"]["weighted_median"]
    r = mr.weighted_median_est(*harm, nboot=1000, seed=1)
    assert r["b"] == pytest.approx(e["b"], rel=1e-10)
    assert r["se"] == pytest.approx(e["se"], rel=0.15)  # bootstrap


@pytest.mark.parametrize("weighted,key", [(True, "weighted_mode"), (False, "simple_mode")])
def test_mode(harm, expected_r, weighted, key):
    e = expected_r["mr"][key]
    r = mr.mode_est(*harm, weighted=weighted, nboot=300, seed=1)
    # point estimate comes from a port of R's density(); grid-level agreement expected
    assert r["b"] == pytest.approx(e["b"], rel=1e-6)
    assert r["se"] == pytest.approx(e["se"], rel=0.3)  # bootstrap, fewer replicates


@pytest.fixture(scope="module")
def region():
    r = pd.read_csv(GOLDEN / "coloc_region.tsv", sep="\t", dtype={"chromosome": str})
    base = r[["variant_id", "chromosome", "base_pair_location"]]
    d1 = base.assign(beta=r["beta1"], standard_error=r["se1"], effect_allele_frequency=r["maf"],
                     n=r["n1"], p_value=0.5)
    d2 = base.assign(beta=r["beta2"], standard_error=r["se2"], effect_allele_frequency=r["maf"],
                     n=r["n2"], p_value=0.5)
    m1 = TraitMeta(trait="q", trait_type="quantitative")
    m2 = TraitMeta(trait="cc", trait_type="binary", sample_prevalence=0.4)
    return d1, m1, d2, m2


def test_coloc_estimated_sdy(region, expected_r):
    d1, m1, d2, m2 = region
    _, s = coloc_abf(d1, m1, d2, m2, Warnings())
    e = expected_r["coloc"]["summary_sdy_estimated"]
    for h in range(5):
        assert s[f"PP.H{h}"] == pytest.approx(e[f"PP.H{h}.abf"], rel=1e-6, abs=1e-14)


def test_coloc_sdy1(region, expected_r, monkeypatch):
    import genor.methods.coloc as C

    d1, m1, d2, m2 = region
    monkeypatch.setattr(C, "sdy_est", lambda *a: 1.0)
    res, s = coloc_abf(d1.drop(columns=["n"]).assign(n=np.nan).fillna({"n": 1}), m1, d2, m2, Warnings())
    e = expected_r["coloc"]["summary"]
    for h in range(5):
        assert s[f"PP.H{h}"] == pytest.approx(e[f"PP.H{h}.abf"], rel=1e-6, abs=1e-14)
    top = res["SNP.PP.H4"].to_numpy()[:5]
    np.testing.assert_allclose(top, expected_r["coloc"]["snp_pp_h4_top"], rtol=1e-6, atol=1e-14)


def test_finemap(region, expected_r, monkeypatch):
    import genor.methods.coloc as C

    d1, m1, _, _ = region
    monkeypatch.setattr(C, "sdy_est", lambda *a: 1.0)
    res, s = finemap_abf(d1, m1, Warnings(), p1=1e-4)
    ours = np.sort(np.append(res["PP"].to_numpy(), s["pp_null"]))[::-1][:5]  # R lists 'null' as a row
    np.testing.assert_allclose(ours, expected_r["coloc"]["finemap_top_pp"], rtol=1e-6)
    assert s["pp_null"] == pytest.approx(expected_r["coloc"]["finemap_null_pp"], rel=1e-6, abs=1e-14)


def test_meta_ivw_random(expected_r):
    m = pd.read_csv(GOLDEN / "meta_studies.tsv", sep="\t")
    tabs = [g.assign(effect_allele="A", other_allele="G", standard_error=g["se"],
                     z=g["beta"] / g["se"], p_value=0.5).drop(columns=["study", "se"])
            for _, g in m.groupby("study")]
    metas = [TraitMeta(trait=f"s{i}", build="GRCh37", ancestry="EUR") for i in range(3)]
    res, _ = meta_analyse(tabs, metas, "ivw", Warnings())
    exp = pd.DataFrame(expected_r["meta"]).set_index("variant_id")
    res = res.set_index("variant_id").loc[exp.index]
    np.testing.assert_allclose(res["beta"], exp["beta"], rtol=1e-10)
    np.testing.assert_allclose(res["standard_error"], exp["se"], rtol=1e-10)
    np.testing.assert_allclose(res["het_Q"], exp["Q"], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(res["het_I2"], exp["I2"], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(res["tau2"], exp["tau2"], rtol=1e-8, atol=1e-14)
    np.testing.assert_allclose(res["beta_random"], exp["beta_random"], rtol=1e-8)
    np.testing.assert_allclose(res["se_random"], exp["se_random"], rtol=1e-8)
