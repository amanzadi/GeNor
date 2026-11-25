# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

import gzip

import numpy as np
import pandas as pd
import pytest

from genor.errors import GenorError, IndividualLevelDataError, NeedsUserDecision, Warnings
from genor.guards import check_ancestry, check_build, refuse_individual_level
from genor.sumstats.columns import detect
from genor.sumstats.harmonise import harmonise
from genor.sumstats.liftover import lift, load_chain
from genor.sumstats.load import load, open_trait, save_trait
from genor.sumstats.qc import lambda_gc, run_qc


def _write(path, df, sep="\t"):
    df.to_csv(path, sep=sep, index=False)
    return path


def _base(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "rsid": [f"rs{i}" for i in range(1, n + 1)], "chr": 1, "pos": np.arange(n) * 1000 + 10_000,
        "a1": rng.choice(list("AC"), n), "a2": "G", "eaf": rng.uniform(0.05, 0.95, n),
        "beta": rng.normal(0, 0.02, n), "se": 0.01, "N": 10_000})


@pytest.mark.parametrize("header,expected", [
    (["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"], {"SNP": "variant_id", "A1": "effect_allele"}),
    (["MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value"], {"P-value": "p_value", "StdErr": "standard_error"}),
    (["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"],
     {"ALLELE1": "effect_allele", "ALLELE0": "other_allele", "LOG10P": "neg_log10_p", "A1FREQ": "effect_allele_frequency"}),
    (["#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"],
     {"alt": "effect_allele", "ref": "other_allele", "sebeta": "standard_error", "af_alt": "effect_allele_frequency"}),
    (["variant_id", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value"],
     {"hm_beta": "beta", "hm_effect_allele": "effect_allele"}),
])
def test_column_detection(header, expected):
    mapping, _ = detect(header)
    for col, canon in expected.items():
        assert mapping[col] == canon


def test_load_basic_and_derivations(tmp_path):
    d = _base()
    f = _write(tmp_path / "s.tsv", d.drop(columns=["se"]).assign(P=0.01))
    w = Warnings()
    df, meta, _ = load(f, "T", tmp_path / "o", ancestry="eur", w=w)
    assert meta.ancestry == "EUR" and meta.n == 10_000
    assert {"z", "standard_error", "p_value"} <= set(df.columns)
    assert any(x.code == "se_from_p" for x in w.items)


def test_odds_ratio_converted(tmp_path):
    d = _base().assign(OR=lambda x: np.exp(x["beta"])).drop(columns=["beta"])
    df, meta, _ = load(_write(tmp_path / "s.tsv", d), "T", tmp_path / "o")
    assert meta.effect_type == "or"
    np.testing.assert_allclose(df["beta"], np.log(d["OR"]))


def test_beta_that_looks_like_or_asks(tmp_path):
    d = _base().assign(beta=lambda x: np.exp(x["beta"]))
    with pytest.raises(NeedsUserDecision):
        load(_write(tmp_path / "s.tsv", d), "T", tmp_path / "o")


def test_chr_pos_ids_and_whitespace(tmp_path):
    d = pd.DataFrame({"MarkerName": ["1:1000:A:G", "2:2000:C:T"] * 60, "Effect": 0.1, "StdErr": 0.05})
    f = tmp_path / "s.txt"
    f.write_text(d.to_string(index=False))
    df, _, _ = load(f, "T", tmp_path / "o")
    assert set(df["chromosome"]) == {"1", "2"}
    assert df["effect_allele"].iloc[0] == "G"


def test_gwas_vcf(tmp_path):
    f = tmp_path / "t.vcf.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tieu-a-2\n")
        for i in range(60):
            fh.write(f"1\t{1000 + i}\trs{i}\tA\tG\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.1:0.02:5.2:0.3:5000:rs{i}\n")
    df, meta, _ = load(f, "T", tmp_path / "o")
    assert len(df) == 60 and meta.n == 5000
    assert df["p_value"].iloc[0] == pytest.approx(10 ** -5.2)


@pytest.mark.parametrize("name", ["data.bed", "x.pgen", "y.bgen", "z.fam"])
def test_refuse_individual_level_ext(tmp_path, name):
    p = tmp_path / name
    p.write_text("x")
    with pytest.raises(IndividualLevelDataError) as e:
        refuse_individual_level(p)
    assert "summary statistics only" in e.value.hint


def test_refuse_genotype_vcf(tmp_path):
    p = tmp_path / "g.vcf"
    p.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                 "1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\n")
    with pytest.raises(IndividualLevelDataError):
        refuse_individual_level(p)


def test_build_inference(example_dir, tmp_path):
    _, meta, w = load(example_dir / "bmi_sim.tsv.gz", "B", tmp_path / "o")
    assert meta.build == "GRCh37" and meta.build_inferred


def test_qc_filters_and_lambda():
    rng = np.random.default_rng(1)
    n = 2000
    d = pd.DataFrame({"variant_id": [f"rs{i}" for i in range(n)], "chromosome": "1",
                      "base_pair_location": np.arange(n) + 1, "effect_allele": "A", "other_allele": "G",
                      "beta": rng.normal(size=n), "standard_error": 1.0,
                      "effect_allele_frequency": rng.uniform(0.001, 0.999, n)})
    d["z"] = d["beta"]
    from scipy import stats
    d["p_value"] = 2 * stats.norm.sf(d["z"].abs())
    d.loc[0, "p_value"] = np.nan
    d.loc[1, "standard_error"] = 0
    d = pd.concat([d, d.iloc[[5]]])
    out, rep, _ = run_qc(d, "GRCh37")
    assert rep["removed"]["bad_p"] == 1 and rep["removed"]["bad_se"] == 1
    assert rep["removed"]["duplicates"] == 1
    assert rep["removed"]["maf_lt_0.01"] > 0
    assert 0.9 < rep["lambda_gc"] < 1.15
    assert lambda_gc(pd.Series(np.full(100, 0.5))) == pytest.approx(1.0, rel=1e-6)


def test_harmonise_cases():
    a = pd.DataFrame({"variant_id": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
                      "effect_allele": ["A", "A", "A", "A", "A", "C"],
                      "other_allele": ["G", "G", "G", "T", "T", "G"],
                      "beta": [1.0] * 6, "effect_allele_frequency": [0.2, 0.2, 0.2, 0.2, 0.45, 0.2]})
    b = pd.DataFrame({"variant_id": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
                      "effect_allele": ["A", "G", "T", "T", "A", "C"],
                      "other_allele": ["G", "A", "C", "A", "T", "A"],
                      "beta": [0.5] * 6, "effect_allele_frequency": [0.2, 0.8, 0.2, 0.2, 0.45, 0.2]})
    h, rep = harmonise(a, b, suffixes=("_a", "_b"))
    got = dict(zip(h["variant_id_a"], h["beta_b"], strict=True))
    assert got["rs1"] == 0.5          # aligned
    assert got["rs2"] == -0.5         # swapped
    assert got["rs3"] == 0.5          # strand flip (T/C vs A/G)
    assert got["rs4"] == 0.5          # palindromic A/T: outcome T at 0.2 = reverse strand, same allele
    assert "rs5" not in got           # palindromic with MAF > 0.42 dropped
    assert "rs6" not in got           # incompatible alleles
    assert rep["palindromic_dropped"] == 1 and rep["allele_mismatch_dropped"] == 1


def test_liftover_point(tmp_path):
    chain = tmp_path / "c.chain"
    chain.write_text("chain 1000 chr1 1000 + 100 600 chr1 1000 + 200 700 1\n200\t50\t50\n250\n\n")
    blocks = load_chain(chain)
    c, p = lift(blocks, pd.Series(["1", "1", "1", "2"]), pd.Series([150, 320, 360, 150]))
    assert p.tolist()[0] == 250      # in first block: +100
    assert pd.isna(p.iloc[1])        # in the gap
    assert p.tolist()[2] == 460      # second block
    assert pd.isna(p.iloc[3])        # chromosome not in chain


def test_guards():
    with pytest.raises(NeedsUserDecision):
        check_build(None, "GRCh37", "panel")
    with pytest.raises(GenorError) as e:
        check_build("GRCh38", "GRCh37", "panel")
    assert e.value.code == "build_mismatch"
    with pytest.raises(GenorError) as e:
        check_ancestry("EAS", "EUR", "panel")
    assert e.value.code == "ancestry_mismatch"
    with pytest.raises(NeedsUserDecision):
        check_ancestry(None, "EUR", "panel")


def test_trait_roundtrip(tmp_path):
    from genor.sumstats.load import TraitMeta

    d = pd.DataFrame({"variant_id": ["rs1"], "p_value": [0.1]})
    save_trait(d, TraitMeta(trait="x", build="GRCh37"), tmp_path / "t")
    df, m = open_trait(tmp_path / "t")
    assert m.build == "GRCh37" and len(df) == 1
    with pytest.raises(GenorError):
        open_trait(tmp_path / "missing")


def test_plink2_glm_effect_allele(tmp_path):
    # A1 is the tested allele; it equals REF on half the rows
    n = 120
    d = pd.DataFrame({"#CHROM": 1, "POS": np.arange(n) + 100, "ID": [f"rs{i}" for i in range(n)],
                      "REF": "A", "ALT": "G", "A1": ["A", "G"] * (n // 2), "A1_FREQ": 0.3,
                      "OBS_CT": 1000, "BETA": 0.1, "SE": 0.05, "P": 0.04})
    w = Warnings()
    df, _, _ = load(_write(tmp_path / "x.glm.linear", d), "T", tmp_path / "o", w=w)
    assert (df["effect_allele"] == d["A1"]).all()
    assert (df["other_allele"] == np.where(d["A1"] == "A", "G", "A")).all()
    assert (df["effect_allele_frequency"] == 0.3).all()


def test_saige_effect_allele_is_allele2(tmp_path):
    n = 120
    d = pd.DataFrame({"CHR": 1, "POS": np.arange(n) + 100, "MarkerID": [f"rs{i}" for i in range(n)],
                      "Allele1": "A", "Allele2": "G", "AC_Allele2": 50, "AF_Allele2": 0.25,
                      "N": 1000, "BETA": 0.1, "SE": 0.05, "p.value": 0.04})
    df, _, _ = load(_write(tmp_path / "saige.txt", d), "T", tmp_path / "o")
    assert (df["effect_allele"] == "G").all() and (df["other_allele"] == "A").all()
    assert (df["effect_allele_frequency"] == 0.25).all()


def test_pleio_template_alignment(tmp_path):
    from genor.methods.pleio import align_to_template, read_template

    (tmp_path / "t.ref").write_text("CHR\tSNP\tGP\tBP\tA1\tA2\n1\trs1\t0\t10\tA\tG\n1\trs2\t0\t20\tA\tT\n"
                                    "1\trs3\t0\t30\tC\tT\n1\trs4\t0\t40\tC\tG\n")
    tmpl = read_template(tmp_path / "t.ref")
    d = pd.DataFrame({"variant_id": ["rs1", "rs2", "rs3"], "effect_allele": ["G", "T", "A"],
                      "other_allele": ["A", "A", "C"], "z": [2.0, 3.0, 1.0], "p_value": [0.05, 0.003, 0.3]})
    out = align_to_template(d, tmpl, Warnings(), "x")
    z = out["zvec"].ravel()
    assert z[0] == -2.0                # swapped
    assert z[1] == -3.0                # palindromic A/T reported as T/A: swap, not strand flip
    assert np.isnan(z[2])              # incompatible alleles
    assert np.isnan(z[3])              # not in the GWAS
    assert out["logpvec"].ravel()[1] == pytest.approx(-np.log10(0.003))
