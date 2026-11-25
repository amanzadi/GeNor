# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""GeNor's LDSC wrapper vs the original ldsc.py (CBIIT/ldsc py3 fork) on identical munged input.

Needs the LDSC reference files (~25 MB download), so it is marked slow:
    uv run pytest -m slow tests/test_ldsc_golden.py
"""

import json

import pandas as pd
import pytest

from genor import refdata
from genor.errors import Warnings
from genor.methods import ldsc
from genor.sumstats.load import TraitMeta

from .conftest import GOLDEN

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def ld(tmp_path_factory):
    import os

    os.environ["GENOR_HOME"] = str(tmp_path_factory.mktemp("ldschome"))
    for name in ("ldsc-ldscores-eur", "ldsc-weights-eur"):
        refdata.fetch(name)
    return ldsc.load_reference()


def _trait(name):
    d = pd.read_csv(GOLDEN / f"ldsc_{name}.sumstats.gz", sep="\t")
    df = pd.DataFrame({"variant_id": d["SNP"], "effect_allele": d["A1"], "other_allele": d["A2"],
                       "z": d["Z"], "n": d["N"]})
    return df, TraitMeta(trait=name, build="GRCh37", ancestry="EUR")


def test_h2_matches_ldsc(ld):
    e = json.loads((GOLDEN / "ldsc_expected.json").read_text())["h2"]
    df, m = _trait("bmi")
    r = ldsc.h2(df, m, ld, Warnings())
    assert r["h2_observed"] == pytest.approx(e["h2"], rel=1e-6)
    assert r["h2_observed_se"] == pytest.approx(e["se"], rel=1e-6)
    assert r["intercept"] == pytest.approx(e["intercept"], rel=1e-6)
    assert r["intercept_se"] == pytest.approx(e["intercept_se"], rel=1e-6)


def test_rg_matches_ldsc(ld):
    e = json.loads((GOLDEN / "ldsc_expected.json").read_text())["rg"]
    d1, m1 = _trait("bmi")
    d2, m2 = _trait("mdd")
    r = ldsc.rg(d1, m1, d2, m2, ld, Warnings())
    assert r["rg"] == pytest.approx(e["rg"], rel=1e-6)
    assert r["rg_se"] == pytest.approx(e["se"], rel=1e-6)
    assert r["h2_trait1"]["h2_observed"] == pytest.approx(e["h2_1"], rel=1e-6)
    assert r["h2_trait2"]["h2_observed_se"] == pytest.approx(e["h2_2_se"], rel=1e-6)
    assert r["gencov"] == pytest.approx(e["gencov"], rel=1e-6)
    assert r["gencov_intercept"] == pytest.approx(e["gencov_intercept"], rel=1e-6)
