import asyncio
import json

import httpx
import pytest

from genor.errors import GenorError
from genor.lookup import http as H
from genor.lookup import services as S
from genor.mcp_server import METHODS, mcp, to_args

from .conftest import GOLDEN

FIX = GOLDEN.parent / "fixtures" / "http"


def test_to_args():
    a = to_args("mr", {"exposure": "/a", "outcome": "/b", "overlap": "none", "exclude_mhc": True,
                       "nboot": 10})
    assert a[:2] == ["mr", "run"] and "--exclude-mhc" in a and a[a.index("--nboot") + 1] == "10"
    a = to_args("meta", {"traits": ["/x", "/y"], "scheme": "ivw"})
    assert a[:4] == ["meta", "run", "/x", "/y"]
    a = to_args("sumstats.load", {"file": "f", "trait": "T", "columns": {"beta": "B"}})
    assert ["--col", "beta=B"] == a[a.index("--col"):a.index("--col") + 2]
    with pytest.raises(GenorError):
        to_args("nope", {})
    with pytest.raises(GenorError):
        to_args("ldsc.rg", {"trait1": "a"})


def test_methods_map_to_real_cli_commands():
    import typer

    from genor.cli import app

    cli = typer.main.get_command(app)
    for method, (path, _, _) in METHODS.items():
        group = cli.commands[path[0]]
        assert path[1] in group.commands, method


def test_tool_list_and_guidance():
    from fastmcp import Client

    async def go():
        async with Client(mcp) as c:
            names = {t.name for t in await c.list_tools()}
            g = (await c.call_tool("get_guidance", {"topic": "genor"})).data
            return names, g

    names, g = asyncio.run(go())
    assert {"doctor", "run_analysis", "job_status", "get_guidance", "sumstats_load"} <= names
    assert "GeNor" in g["text"]


# ---------------------------------------------------------------- recorded HTTP fixtures


def _transport():
    def handler(req: httpx.Request) -> httpx.Response:
        host = req.url.host
        if "opentargets" in host:
            q = json.loads(req.content)["query"]
            key = "ot_search" if "search(" in q else "ot_target" if "target(" in q else "ot_other"
        elif "ensembl" in host:
            key = "ensembl_overlap" if "overlap" in req.url.path else "ensembl_vep"
        elif "string-db" in host:
            key = "string_ppi" if "ppi_enrichment" in req.url.path else "string_network"
        elif "biit" in host:
            key = "gprofiler"
        else:
            return httpx.Response(404)
        return httpx.Response(200, content=(FIX / f"{key}.json").read_bytes(),
                              headers={"content-type": "application/json"})
    return httpx.MockTransport(handler)


@pytest.fixture
def mocked(monkeypatch):
    monkeypatch.setattr(H, "TRANSPORT", _transport())


def test_open_targets(mocked):
    hits = S.ot_search("PCSK9")
    assert hits and "id" in hits[0]
    t = S.ot_target("ENSG00000169174")
    assert t["approvedSymbol"] == "PCSK9"


def test_ensembl(mocked):
    g = S.genes_in_region("1", 55_000_000, 55_100_000, "GRCh38")
    assert any(x["symbol"] == "PCSK9" for x in g)
    v = S.vep(["rs11591147"])
    assert v[0]["most_severe_consequence"]


def test_string_and_gprofiler(mocked):
    s = S.string_network(["PCSK9", "LDLR", "APOB"])
    assert s["n_edges"] >= 1 and s["string_version"] == "12.0"
    e = S.gprofiler(["PCSK9", "LDLR", "APOB"])
    assert e["n_terms"] >= 1


def test_http_errors(monkeypatch):
    monkeypatch.setattr(H, "TRANSPORT", httpx.MockTransport(lambda r: httpx.Response(400, text="bad")))
    with pytest.raises(GenorError) as e:
        S.ot_search("x")
    assert e.value.code == "remote_bad_request"


@pytest.mark.live
def test_live_apis():
    """Nightly: the real services still answer in the shape GeNor expects."""
    assert S.ot_search("PCSK9", "target")
    assert S.genes_in_region("1", 55_000_000, 55_100_000, "GRCh38")
    assert S.string_network(["PCSK9", "LDLR"])["n_edges"] >= 1
    assert S.gprofiler(["PCSK9", "LDLR", "APOB"])["n_terms"] >= 1
