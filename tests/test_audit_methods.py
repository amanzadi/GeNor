# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Regression tests for audited bugs in sumstats loading, liftover, meta, pleio and LDSC."""

import bz2
import lzma

import numpy as np
import pandas as pd
import pytest

from genor.errors import GenorError, Warnings
from genor.guards import norm_build
from genor.sumstats.liftover import lift, load_chain
from genor.sumstats.load import TraitMeta, load, standardise


def _tsv(n=120, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"SNP": [f"rs{i}" for i in range(n)], "CHR": 1, "BP": np.arange(n) * 1000 + 10_000,
                         "A1": "A", "A2": "G", "BETA": rng.normal(0, 0.02, n), "SE": 0.01,
                         "P": rng.uniform(0.01, 1, n), "N": 10_000})


# A
def test_norm_build_rejects_unknown():
    assert norm_build(None) is None
    assert norm_build("hg19") == "GRCh37" and norm_build("b38") == "GRCh38"
    with pytest.raises(GenorError) as e:
        norm_build("hg18")
    assert e.value.code == "bad_build"


# B
def test_p_underflow_kept():
    raw = pd.DataFrame({"SNP": ["rs1", "rs2", "rs3"], "A1": "A", "A2": "G",
                        "BETA": ["2.0", "0.1", "0"], "SE": ["0.04", "0.05", "0.05"],
                        "P": ["1e-400", "0.05", "1e-400"]})
    w = Warnings()
    df, _, _ = standardise(raw, None, w)
    assert df["p_value"].iloc[0] == 1e-300
    assert df["p_value"].iloc[2] == 0.0  # z = 0: no usable z, left for QC to drop
    assert any(x["code"] == "p_underflow" for x in w.to_list())


# C
def test_liftover_chrx_and_mt(tmp_path):
    chain = tmp_path / "c.chain"
    chain.write_text("chain 1000 chrX 1000 + 100 600 chrX 1000 + 200 700 1\n500\n\n"
                     "chain 1000 chrM 1000 + 100 600 chrM 1000 + 200 700 2\n500\n\n")
    c, p = lift(load_chain(chain), pd.Series(["23", "26"]), pd.Series([150, 150]))
    assert c.tolist() == ["23", "26"] and p.tolist() == [250, 250]


# D
def _study(ids, chrom=True, seed=0):
    rng = np.random.default_rng(seed)
    n = len(ids)
    d = pd.DataFrame({"variant_id": ids, "effect_allele": "A", "other_allele": "G",
                      "beta": rng.normal(0, 0.1, n), "standard_error": 0.05, "n": 1000.0})
    d["z"] = d["beta"] / d["standard_error"]
    if chrom:
        d["chromosome"] = "1"
        d["base_pair_location"] = np.arange(n) + 100
    return d


def test_meta_mixed_keys():
    from genor.methods.meta import meta_analyse

    metas = [TraitMeta(trait="a", build="GRCh37"), TraitMeta(trait="b", build="GRCh37")]
    a = _study([f"rs{i}" for i in range(10)])
    b = _study([f"1:{i + 100}" for i in range(10)], seed=1)
    _, s = meta_analyse([a, b], metas)
    assert s["n_in_all_studies"] == 10
    with pytest.raises(GenorError) as e:
        meta_analyse([a, _study([f"rs{i}" for i in range(10)], chrom=False)
                      .assign(variant_id=[f"1:{i}" for i in range(10)])], metas)
    assert e.value.code == "meta_needs_position"
    with pytest.raises(GenorError) as e:
        meta_analyse([a, _study([f"rs{i}" for i in range(100, 110)])], metas)
    assert e.value.code == "meta_no_overlap"


# E
def test_pleio_name_collision(tmp_path):
    from genor.methods import pleio

    m1 = TraitMeta(trait="MDD-2019", ancestry="EUR")
    m2 = TraitMeta(trait="MDD_2019", ancestry="EUR")
    with pytest.raises(GenorError) as e:
        pleio.run(pd.DataFrame(), m1, pd.DataFrame(), m2, tmp_path, Warnings(), prepare_only=True,
                  template=tmp_path / "missing.txt", ref_mat=tmp_path / "missing.mat")
    assert e.value.code == "pleio_same_trait_name"


# F
def test_ldsc_munge_drops_missing_n():
    from genor.methods.ldsc import LDRef, munge

    n = 20
    d = pd.DataFrame({"variant_id": [f"rs{i}" for i in range(n)], "effect_allele": "A",
                      "other_allele": "G", "z": 1.0, "n": [1000.0] * (n - 2) + [np.nan, np.nan]})
    ld = LDRef(pd.DataFrame(), pd.DataFrame(), 1.0, "EUR", "GRCh37")
    out = munge(d, TraitMeta(trait="t"), ld, Warnings(), "t")
    assert len(out) == n - 2 and out["N"].notna().all()


# G
def test_sites_only_vcf(tmp_path):
    f = tmp_path / "x.vcf"
    f.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "1\t100\trs1\tA\tG\t.\tPASS\t.\n")
    with pytest.raises(GenorError) as e:
        load(f, "t", tmp_path / "o")
    assert e.value.code == "vcf_not_gwas_vcf"


# H
@pytest.mark.parametrize("ext,opener", [(".bz2", bz2.open), (".xz", lzma.open)])
def test_bz2_xz(tmp_path, ext, opener):
    f = tmp_path / f"x.tsv{ext}"
    with opener(f, "wt") as fh:
        _tsv().to_csv(fh, sep="\t", index=False)
    df, _, _ = load(f, "t", tmp_path / "o", build="GRCh37")
    assert len(df) == 120


def test_zst_refused(tmp_path):
    f = tmp_path / "x.tsv.zst"
    f.write_bytes(b"\x28\xb5\x2f\xfd")
    with pytest.raises(GenorError) as e:
        load(f, "t", tmp_path / "o")
    assert e.value.code == "unsupported_compression"
