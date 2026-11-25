# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""End-to-end on the bundled example, without downloads (distance clumping, no LDSC)."""

import json
import subprocess
import sys

import pandas as pd
import pytest

from genor.commands import coloc_run, finemap_abf, meta_run, mr_run, report_build


def test_mr_recovers_planted_effect(loaded_pair, tmp_path, example_dir):
    bmi, mdd = loaded_pair
    r = mr_run(bmi, mdd, tmp_path / "mr", clump_mode="distance", overlap="none", nboot=200)
    assert r["status"] == "ok"
    s = r["summary"]
    truth = json.loads((example_dir / "truth.json").read_text())["causal_effect"]
    ivw = next(e for e in s["estimates"] if e["method"] == "IVW")
    assert abs(ivw["b"] - truth) < 3 * ivw["se"]
    assert s["steiger"]["correct_causal_direction"]
    assert s["strobe_mr"]["filled"] + s["strobe_mr"]["partial"] + s["strobe_mr"]["human"] == 20
    codes = {w["code"] for w in r["warnings"]}
    assert "distance_clumping" in codes
    manifest = json.loads((tmp_path / "mr" / "manifest.json").read_text())
    assert manifest["inputs"][0]["sha256"] and any(c["key"] == "strobe-mr" for c in manifest["citations"])
    assert (tmp_path / "mr" / "strobe_mr.md").exists()


def test_mr_asks_about_overlap(loaded_pair, tmp_path):
    bmi, mdd = loaded_pair
    r = mr_run(bmi, mdd, tmp_path / "mr", clump_mode="distance", nboot=50)
    assert any("overlap" in a["question"] for a in r["ask_user"])


def test_plink_clump_without_panel_is_actionable(loaded_pair, tmp_path):
    from genor.errors import GenorError

    bmi, mdd = loaded_pair
    with pytest.raises(GenorError) as e:
        mr_run(bmi, mdd, tmp_path / "mr", clump_mode="plink", overlap="none")
    assert e.value.code == "refdata_missing"
    assert "g1000-eur" in e.value.hint
    m = json.loads((tmp_path / "mr" / "manifest.json").read_text())
    assert m["status"] == "error" and m["error"]["code"] == "refdata_missing"


def test_finemap_and_coloc(loaded_pair, tmp_path):
    bmi, mdd = loaded_pair
    lead = pd.read_parquet(bmi / "sumstats.parquet").nsmallest(1, "p_value")["variant_id"].iloc[0]
    f = finemap_abf(bmi, tmp_path / "fm", lead=lead, window_kb=1000)
    assert f["summary"]["top_variant"] == lead
    c = coloc_run(bmi, mdd, tmp_path / "co", lead=lead, window_kb=1000)
    assert abs(sum(c["summary"][f"PP.H{i}"] for i in range(5)) - 1) < 1e-9


def test_meta_of_two_halves(loaded_pair, tmp_path):
    bmi, _ = loaded_pair
    r = meta_run([bmi, bmi], tmp_path / "meta")
    out = pd.read_csv(tmp_path / "meta" / "meta_results.tsv.gz", sep="\t")
    src = pd.read_parquet(bmi / "sumstats.parquet").set_index("variant_id")
    o = out.set_index("variant_id").loc[src.index[:100]]
    # same study twice: same beta, SE / sqrt(2), I2 = 0
    pd.testing.assert_series_equal(o["beta"], src["beta"].iloc[:100], check_names=False, rtol=1e-9)
    assert (o["het_I2"] == 0).all()
    assert r["summary"]["n_studies"] == 2


def test_report(loaded_pair, tmp_path):
    bmi, mdd = loaded_pair
    mr_run(bmi, mdd, tmp_path / "mr", clump_mode="distance", overlap="none", nboot=50)
    r = report_build([bmi, tmp_path / "mr"], tmp_path / "report.md")
    md = (tmp_path / "report.md").read_text()
    assert "STROBE-MR" in md and "References" in md and "not a medical device" in md
    assert (tmp_path / "report.html").exists() and r["n_runs"] == 2


def _cli(*args, cwd):
    p = subprocess.run([sys.executable, "-m", "genor", *args, "--json"], capture_output=True,
                       text=True, cwd=cwd)
    return p.returncode, json.loads(p.stdout) if p.stdout.strip() else None


def test_cli_json_and_exit_codes(tmp_path):
    code, out = _cli("doctor", cwd=tmp_path)
    assert code == 0 and out["status"] == "ok" and "blocked_until_setup" in out
    code, out = _cli("example", "--out", "ex", cwd=tmp_path)
    assert code == 0
    (tmp_path / "geno.bed").write_text("x")
    code, out = _cli("sumstats", "load", "geno.bed", "--trait", "x", cwd=tmp_path)
    assert code == 2 and out["error"]["code"] == "individual_level_refused"
    code, out = _cli("sumstats", "load", "ex/bmi_sim.tsv.gz", "--trait", "B", "--no-plots",
                     "--out", "b", cwd=tmp_path)
    assert code == 0 and out["summary"]["build"] == "GRCh37"
    code, out = _cli("clump", "run", "b", "--mode", "plink", cwd=tmp_path)
    assert code == 3 and out["status"] == "needs_user_decision"  # ancestry unknown -> ask


def test_jobs_roundtrip(tmp_path):
    from genor import jobs

    j = jobs.submit(["doctor"], cwd=tmp_path)
    r = jobs.wait(j["job_id"], 60)
    assert r["status"] == "ok" and r["job_id"] == j["job_id"]
    assert jobs.status(j["job_id"])["status"] == "done"
