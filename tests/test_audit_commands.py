# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Regression tests for audited bugs in commands, MCP, jobs, runs, CLI, backends and lookups."""

import asyncio
import json
import os
import signal
import subprocess

import httpx
import numpy as np
import pandas as pd
import pytest

from genor.errors import GenorError
from genor.sumstats.load import TraitMeta, open_trait, save_trait


def _strict(text: str):
    def bad(c):
        raise ValueError(f"non-standard JSON constant {c}")
    return json.loads(text, parse_constant=bad)


# ---------------------------------------------------------------- A: MCP report_build


def test_mcp_report_build_confined_to_workspace(loaded_pair, tmp_path, monkeypatch):
    from fastmcp import Client

    from genor.mcp_server import mcp

    ws = tmp_path / "ws"
    monkeypatch.setenv("GENOR_WORKDIR", str(ws))
    bmi, _ = loaded_pair
    outside = tmp_path / "outside.md"

    async def go():
        async with Client(mcp) as c:
            bad = (await c.call_tool("report_build", {"run_dirs": [str(bmi)], "out": str(outside)})).data
            ok = (await c.call_tool("report_build", {"run_dirs": [str(bmi)], "out": "r.md",
                                                     "title": "T\nexport EVIL=1\r"})).data
            return bad, ok

    bad, ok = asyncio.run(go())
    assert bad["status"] == "error" and bad["error"]["code"] == "bad_param"
    assert not outside.exists()
    assert ok["status"] == "ok"
    assert "export EVIL=1" not in (ws / "r.md").read_text().splitlines()


# ---------------------------------------------------------------- B, K: jobs


def test_job_id_traversal_rejected(tmp_path):
    from genor import jobs
    from genor.config import jobs_dir

    target = jobs_dir().parent / "secret"
    target.mkdir()
    (target / "stdout.log").write_text("secret")
    (target / "stderr.log").write_text("secret")
    with pytest.raises(GenorError) as e:
        jobs.logs("../secret", 5)
    assert e.value.code == "job_not_found"


def test_cancel_does_not_signal_finished_job(tmp_path):
    from genor import jobs
    from genor.config import jobs_dir

    other = subprocess.Popen(["sleep", "30"], start_new_session=True)  # stands in for a reused PID
    try:
        d = jobs_dir() / "20260101-000000-abcdef"
        d.mkdir()
        (d / "stdout.log").write_text(json.dumps({"status": "ok", "run_dir": "x"}))
        (d / "stderr.log").write_text("")
        (d / "job.json").write_text(json.dumps({"id": d.name, "status": "running", "pid": other.pid,
                                                "label": "doctor", "submitted": "t"}))
        r = jobs.cancel(d.name)
        assert r["status"] == "done"
        assert other.poll() is None
    finally:
        os.killpg(other.pid, signal.SIGKILL)
        other.wait()


# ---------------------------------------------------------------- C: backend_add


def test_backend_add_requires_registry_executable_name(tmp_path):
    from genor import backends

    marker = tmp_path / "ran"
    exe = tmp_path / "evil"
    exe.write_text(f"#!/bin/sh\ntouch {marker}\n")
    exe.chmod(0o755)
    with pytest.raises(GenorError) as e:
        backends.add("magma", exe)
    assert e.value.code == "backend_not_found"
    assert not marker.exists() and "magma" not in backends.registrations()


# ---------------------------------------------------------------- D: meta keeps binary


def test_meta_of_binary_traits_stays_binary(loaded_pair, tmp_path):
    from genor.commands import meta_run

    _, mdd = loaded_pair
    meta_run([mdd, mdd], tmp_path / "meta")
    assert open_trait(tmp_path / "meta")[1].trait_type == "binary"


# ---------------------------------------------------------------- E: coloc keeps palindromes


def _region_trait(path, trait, seed, with_ids=True):
    rng = np.random.default_rng(seed)
    k = 30
    d = pd.DataFrame({
        "chromosome": "1", "base_pair_location": np.arange(1, k + 1) * 1000,
        "effect_allele": ["A", "G"] * (k // 2), "other_allele": ["T", "A"] * (k // 2),
        "beta": rng.normal(0, 0.02, k), "standard_error": 0.01, "p_value": 0.5, "n": 10_000.0,
        "effect_allele_frequency": 0.45,
    })
    if with_ids:
        d.insert(0, "variant_id", [f"rs{i}" for i in range(k)])
    save_trait(d, TraitMeta(trait=trait, build="GRCh37", n=10_000), path)
    return path


def test_coloc_keeps_palindromic_snps(tmp_path):
    from genor.commands import coloc_run

    a = _region_trait(tmp_path / "a", "A", 1)
    b = _region_trait(tmp_path / "b", "B", 2)
    r = coloc_run(a, b, tmp_path / "co", region="1:1-100000")
    assert r["summary"]["nsnps"] == 30


def test_coloc_without_variant_id(tmp_path):
    from genor.commands import coloc_run

    a = _region_trait(tmp_path / "a", "A", 1, with_ids=False)
    b = _region_trait(tmp_path / "b", "B", 2, with_ids=False)
    r = coloc_run(a, b, tmp_path / "co", region="1:1-100000")
    assert r["summary"]["nsnps"] == 30
    with pytest.raises(GenorError) as e:
        coloc_run(a, b, tmp_path / "co2", lead="rs1")
    assert e.value.code == "bad_param"


# ---------------------------------------------------------------- F, G: MR


def test_mr_needs_beta_se(loaded_pair, tmp_path):
    from genor.commands import mr_run

    bmi, mdd = loaded_pair
    d, m = open_trait(bmi)
    save_trait(d.drop(columns=["beta", "standard_error"]), m, tmp_path / "z_only")
    with pytest.raises(GenorError) as e:
        mr_run(tmp_path / "z_only", mdd, tmp_path / "mr", clump_mode="distance", overlap="none")
    assert e.value.code == "mr_needs_beta_se"


def test_mr_records_ld_panel(loaded_pair, tmp_path, monkeypatch):
    import genor.methods.clump as C
    from genor.commands import mr_run

    bmi, mdd = loaded_pair
    real = C.clump

    def fake(df, meta, work, **kw):
        res, info = real(df, meta, work, **{**kw, "mode": "distance"})
        return res, {**info, "mode": "plink", "panel": "g1000-eur"}

    monkeypatch.setattr(C, "clump", fake)
    mr_run(bmi, mdd, tmp_path / "mr", clump_mode="plink", overlap="none", nboot=20)
    m = json.loads((tmp_path / "mr" / "manifest.json").read_text())
    assert [r["name"] for r in m["refdata"]] == ["g1000-eur"]


# ---------------------------------------------------------------- H: list_runs root


def test_list_runs_searches_root_only(tmp_path, monkeypatch):
    from genor.commands import list_runs

    for d in ("a/run1", "b/run2"):
        (tmp_path / d).mkdir(parents=True)
        (tmp_path / d / "manifest.json").write_text(json.dumps({"command": "x", "status": "ok"}))
    monkeypatch.chdir(tmp_path / "b")
    runs = list_runs(tmp_path / "a")["runs"]
    assert [r["run_dir"] for r in runs] == [str(tmp_path / "a" / "run1")]


# ---------------------------------------------------------------- I: no NaN in JSON


def test_no_nan_tokens_in_json(tmp_path, monkeypatch, capsys):
    from genor import cli
    from genor.runs import Run, _jsonable

    val = {"a": float("nan"), "b": [np.float64("nan"), {"c": float("inf")}], "d": np.float32(1.5)}
    assert _jsonable(val) == {"a": None, "b": [None, {"c": None}], "d": 1.5}
    with Run("x", tmp_path / "run") as run:
        run.summary = val
        run.meta = val
    _strict((tmp_path / "run" / "result.json").read_text())
    _strict((tmp_path / "run" / "manifest.json").read_text())
    assert run.result()["summary"]["a"] is None
    monkeypatch.setattr(cli, "JSON_MODE", True)
    cli.emit(val)
    assert _strict(capsys.readouterr().out)["b"][0] is None


# ---------------------------------------------------------------- J: lookups


def _vep_ok(req):
    return httpx.Response(200, json=[{"id": "rs1", "most_severe_consequence": "x"}])


def test_lookup_build_alias_and_http_errors(monkeypatch):
    from genor.lookup import http as H
    from genor.lookup import services as S

    seen = []
    monkeypatch.setattr(H, "TRANSPORT", httpx.MockTransport(lambda r: (seen.append(r.url.host), _vep_ok(r))[1]))
    assert S.vep(["rs1"], "hg38")[0]["id"] == "rs1" and seen == ["rest.ensembl.org"]
    with pytest.raises(GenorError):
        S.genes_in_region("1", 1, 1000, "hg99")

    monkeypatch.setattr(H, "TRANSPORT", httpx.MockTransport(
        lambda r: httpx.Response(503, headers={"Retry-After": "Wed, 21 Oct 2015 07:28:00 GMT"})))
    with pytest.raises(GenorError) as e:
        S.vep(["rs1"], "GRCh38")
    assert e.value.code == "remote_unavailable"

    monkeypatch.setattr(H, "TRANSPORT", httpx.MockTransport(lambda r: httpx.Response(200, text="<html>")))
    with pytest.raises(GenorError) as e:
        S.vep(["rs1"], "GRCh38")
    assert e.value.code == "remote_unavailable"
