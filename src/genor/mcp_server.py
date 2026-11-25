# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""GeNor MCP server (stdio).

Tools mirror the CLI. Analyses run as background jobs: a tool waits up to
``wait_seconds`` (default 40, below Codex's 60 s tool timeout) and then returns either
the result or a ``job_id`` to poll with ``job_status`` / ``job_result``.

Paths: give absolute paths. Relative paths resolve against the workspace, which is
``$GENOR_WORKDIR`` if set, else the server's working directory (or ``~/genor-workspace``
when that is ``/``).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from fastmcp import FastMCP

from genor import RESEARCH_USE_NOTICE, __version__, jobs
from genor.errors import GenorError

INSTRUCTIONS = f"""GeNor runs post-GWAS analyses on GWAS summary statistics.
Start every session with `doctor`, then read `get_guidance` (topic 'genor') unless the
GeNor skills are already loaded. Policy: summary statistics only; the user supplies their
own GWAS files (use `sumstats_sources` to tell them where to look); never guess build,
ancestry, sample overlap, or exposure/outcome definitions - ask the human. Results with
`ask_user` entries need the human's answer before you interpret them. Long analyses return
a job_id: poll `job_status`, then `job_result`. {RESEARCH_USE_NOTICE}"""

mcp = FastMCP("genor", instructions=INSTRUCTIONS, version=__version__)

METHODS: dict[str, tuple[list[str], list[str], str]] = {
    "sumstats.load": (["sumstats", "load"], ["file"],
                      "params: file, trait, out, build, ancestry, n, n_cases, n_controls, "
                      "population_prevalence, effect_type, columns{canonical:filecol}, maf, info, drop_mhc"),
    "sumstats.liftover": (["sumstats", "liftover"], ["trait"], "params: trait, to (GRCh37|GRCh38), out"),
    "meta": (["meta", "run"], ["traits"], "params: traits[list of trait dirs], scheme (ivw|samplesize), name, out"),
    "ldsc.h2": (["ldsc", "h2"], ["trait"], "params: trait, pop_prev, intercept, ref, w, out"),
    "ldsc.rg": (["ldsc", "rg"], ["trait1", "trait2"], "params: trait1, trait2, pop_prev1, pop_prev2, out"),
    "clump": (["clump", "run"], ["trait"], "params: trait, mode (plink|distance), preset (loci|mr), p1, r2, kb, panel, out"),
    "mr": (["mr", "run"], [], "params: exposure, outcome, p_threshold, clump_mode (plink|distance), panel, "
           "instruments, exclude_mhc, overlap (none|partial|substantial|unknown), nboot, seed, out"),
    "pleio": (["pleio", "run"], ["trait1", "trait2"], "params: trait1, trait2, stat (condfdr|conjfdr), fdr, "
              "randprune_n, seed, keep_mhc, prepare_only, ref_mat, template, out"),
    "finemap.abf": (["finemap", "abf"], ["trait"], "params: trait, region (chr:start-end) or lead, window_kb, p1, level, out"),
    "coloc.abf": (["coloc", "abf"], ["trait1", "trait2"], "params: trait1, trait2, region or lead, window_kb, p1, p2, p12, out"),
    "magma.gene": (["magma", "gene"], ["trait"], "params: trait, gene_loc, ref_bfile, up_kb, down_kb, out"),
    "magma.sets": (["magma", "sets"], ["genes_raw"], "params: genes_raw, set_annot or covar, out"),
}


def workspace() -> Path:
    if w := os.environ.get("GENOR_WORKDIR"):
        p = Path(w).expanduser()
    else:
        p = Path.cwd()
        if str(p) == "/" or not os.access(p, os.W_OK):
            p = Path.home() / "genor-workspace"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _err(e: GenorError) -> dict[str, Any]:
    from genor.errors import NeedsUserDecision

    return {"status": "needs_user_decision" if isinstance(e, NeedsUserDecision) else "error",
            "error": e.to_dict()}


def to_args(method: str, params: dict[str, Any]) -> list[str]:
    if method not in METHODS:
        raise GenorError(f"Unknown method '{method}'", code="unknown_method",
                         hint="Methods: " + ", ".join(METHODS))
    prefix, positional, _ = METHODS[method]
    params = dict(params)
    args = list(prefix)
    for p in positional:
        if p not in params:
            raise GenorError(f"Missing parameter '{p}' for {method}", code="bad_param",
                             hint=METHODS[method][2])
        v = params.pop(p)
        args += [str(x) for x in v] if isinstance(v, list) else [str(v)]
    for k, v in params.items():
        if v is None:
            continue
        flag = "--" + k.replace("_", "-")
        if k == "columns" and isinstance(v, dict):
            for ck, cv in v.items():
                args += ["--col", f"{ck}={cv}"]
        elif isinstance(v, bool):
            args.append(flag if v else "--no-" + k.replace("_", "-"))
        elif isinstance(v, list):
            for x in v:
                args += [flag, str(x)]
        else:
            args += [flag, str(v)]
    return args


def _submit(args: list[str], wait_seconds: float) -> dict[str, Any]:
    j = jobs.submit(args, cwd=workspace())
    return jobs.wait(j["job_id"], max(0.0, min(wait_seconds, 50.0)))


# ------------------------------------------------------------------ tools


@mcp.tool
def doctor() -> dict[str, Any]:
    """Platform, available backends, downloaded reference data, and what each missing piece blocks. Call first."""
    from genor.commands import doctor as d

    out = d()
    out["workspace"] = str(workspace())
    return out


@mcp.tool
def get_guidance(topic: str | None = None) -> dict[str, Any]:
    """GeNor skill text for clients without skills support. topic: 'genor' (start here), a workflow such as
    'mendelian-randomisation', or '<skill>/references/<file>.md'. No topic lists all topics."""
    from genor import skills

    try:
        return skills.guidance(topic)
    except GenorError as e:
        return _err(e)


@mcp.tool
def sumstats_sources(name: str | None = None) -> dict[str, Any]:
    """Where users can download GWAS summary statistics (GWAS Catalog, PGC, OpenGWAS, FinnGen, Pan-UKB, ...).
    GeNor never downloads GWAS itself; the user downloads and gives a local path."""
    from genor.commands import sources

    try:
        return sources(name)
    except GenorError as e:
        return _err(e)


@mcp.tool
def example_data(out_dir: str = "genor_example") -> dict[str, Any]:
    """Copy the bundled simulated GWAS pair (BMI_sim exposure, MDD_sim outcome) into the workspace."""
    from genor.commands import example

    return example(workspace() / out_dir)


@mcp.tool
def sumstats_load(file: str, trait: str, out: str | None = None, build: str | None = None,
                  ancestry: str | None = None, n: float | None = None, n_cases: float | None = None,
                  n_controls: float | None = None, population_prevalence: float | None = None,
                  effect_type: str | None = None, columns: dict[str, str] | None = None,
                  maf: float = 0.01, info: float = 0.8, drop_mhc: bool = False,
                  wait_seconds: float = 40) -> dict[str, Any]:
    """Load a user's summary statistics file: detect columns, standardise, infer build, QC, plots.
    Returns the trait directory (pass it to other tools). Do not guess ancestry or N: ask the user."""
    params = {"file": file, "trait": trait, "out": out or f"genor_runs/{trait}", "build": build,
              "ancestry": ancestry, "n": n, "n_cases": n_cases, "n_controls": n_controls,
              "population_prevalence": population_prevalence, "effect_type": effect_type,
              "columns": columns, "maf": maf, "info": info, "drop_mhc": drop_mhc}
    try:
        return _submit(to_args("sumstats.load", params), wait_seconds)
    except GenorError as e:
        return _err(e)


@mcp.tool
def run_analysis(method: str, params: dict[str, Any], wait_seconds: float = 40) -> dict[str, Any]:
    """Run an analysis as a background job. method is one of:
    sumstats.load, sumstats.liftover, meta, ldsc.h2, ldsc.rg, clump, mr, pleio, finemap.abf, coloc.abf, magma.gene, magma.sets.
    params use the CLI option names with underscores (see list_methods). Returns the result if done within
    wait_seconds, else a job_id to poll."""
    try:
        return _submit(to_args(method, params), wait_seconds)
    except GenorError as e:
        return _err(e)


@mcp.tool
def list_methods() -> dict[str, Any]:
    """Methods accepted by run_analysis and their parameters."""
    return {m: {"cli": " ".join(v[0]), "params": v[2]} for m, v in METHODS.items()}


@mcp.tool
def job_status(job_id: str) -> dict[str, Any]:
    """State of a background job: running, done, failed or cancelled."""
    try:
        return jobs.status(job_id)
    except GenorError as e:
        return _err(e)


@mcp.tool
def job_result(job_id: str, wait_seconds: float = 0) -> dict[str, Any]:
    """Result JSON of a finished job (optionally wait up to wait_seconds)."""
    try:
        return jobs.wait(job_id, min(wait_seconds, 50)) if wait_seconds else jobs.result(job_id)
    except GenorError as e:
        return _err(e)


@mcp.tool
def job_logs(job_id: str, n: int = 40) -> dict[str, Any]:
    """Last log lines of a job."""
    try:
        return jobs.logs(job_id, n)
    except GenorError as e:
        return _err(e)


@mcp.tool
def job_cancel(job_id: str) -> dict[str, Any]:
    """Cancel a running job."""
    try:
        return jobs.cancel(job_id)
    except GenorError as e:
        return _err(e)


@mcp.tool
def list_runs() -> dict[str, Any]:
    """Run directories (with manifest.json) in the workspace."""
    from genor.commands import list_runs as lr

    cwd = os.getcwd()
    try:
        os.chdir(workspace())
        out = lr()
        out["workspace"] = str(workspace())
        return out
    finally:
        os.chdir(cwd)


@mcp.tool
def read_run(run_dir: str) -> dict[str, Any]:
    """Short result of a run directory (result.json) plus the list of files it wrote."""
    import json

    d = Path(run_dir)
    if not d.is_absolute():
        d = workspace() / d
    r = d / "result.json"
    if not r.exists():
        return {"status": "error", "error": {"code": "not_a_run_dir", "message": f"No result.json in {d}"}}
    return {**json.loads(r.read_text()), "all_files": sorted(p.name for p in d.iterdir())}


@mcp.tool
def refdata_list() -> dict[str, Any]:
    """Reference datasets (LD scores, 1000G panels, chain files), sizes, and whether downloaded."""
    from genor import refdata

    return {"datasets": refdata.status()}


@mcp.tool
def refdata_fetch(name: str, confirmed_by_user: bool = False, allow_unpinned: bool = False,
                  wait_seconds: float = 40) -> dict[str, Any]:
    """Download a reference dataset. Datasets over 100 MB require the user's confirmation: ask them,
    then call with confirmed_by_user=true. Datasets without a pinned checksum (pleioFDR files) also
    need allow_unpinned=true, again only after the user agrees."""
    args = ["refdata", "fetch", name] + (["--yes"] if confirmed_by_user else []) \
        + (["--allow-unpinned"] if allow_unpinned else [])
    return _submit(args, wait_seconds)


@mcp.tool
def backend_list() -> dict[str, Any]:
    """External tools (PLINK 2, GCTA, SMR, MAGMA, pleioFDR, R) and how to enable missing ones."""
    from genor import backends

    return {"platform": backends.platform_key(), "backends": backends.status()}


@mcp.tool
def backend_install(name: str, wait_seconds: float = 40) -> dict[str, Any]:
    """Download a pinned redistributable binary: plink2, gcta or smr."""
    return _submit(["backend", "install", name], wait_seconds)


@mcp.tool
def backend_add(name: str, path: str) -> dict[str, Any]:
    """Register a tool the user installed themselves (e.g. MAGMA, whose licence forbids bundling).
    path is the executable or its directory."""
    from genor import backends

    try:
        return {"status": "ok", **backends.add(name, path)}
    except GenorError as e:
        return _err(e)


@mcp.tool
def lookup_open_targets(kind: str, query: str, n: int = 10, entity: str | None = None) -> dict[str, Any]:
    """Open Targets Platform. kind: search (free text), target (Ensembl gene ID), disease (EFO/MONDO ID),
    variant (GRCh38 CHR_POS_REF_ALT: credible sets + L2G), coloc (studyLocusId)."""
    from genor.lookup import services as s

    fn = {"search": lambda: s.ot_search(query, entity, n), "target": lambda: s.ot_target(query, n),
          "disease": lambda: s.ot_disease_targets(query, n),
          "variant": lambda: s.ot_variant_credible_sets(query, n),
          "coloc": lambda: s.ot_study_locus_coloc(query, n)}
    if kind not in fn:
        return {"status": "error", "error": {"code": "bad_param", "message": f"kind must be one of {list(fn)}"}}
    try:
        return {"status": "ok", "result": fn[kind]()}
    except GenorError as e:
        return _err(e)


@mcp.tool
def lookup_variants(rsids: list[str], build: str = "GRCh38") -> dict[str, Any]:
    """Ensembl VEP: most severe consequence and genes for up to 200 rsIDs."""
    from genor.lookup import services as s

    try:
        return {"status": "ok", "variants": s.vep(rsids, build)}
    except GenorError as e:
        return _err(e)


@mcp.tool
def lookup_genes_in_region(region: str, build: str = "GRCh37") -> dict[str, Any]:
    """Genes overlapping chr:start-end (Ensembl), for positional mapping of loci."""
    from genor.lookup import services as s
    from genor.methods.coloc import parse_region

    try:
        c, a, b = parse_region(region)
        return {"status": "ok", "genes": s.genes_in_region(c, a, b, build)}
    except GenorError as e:
        return _err(e)


@mcp.tool
def lookup_string(genes: list[str], required_score: int = 700) -> dict[str, Any]:
    """STRING v12 interactions among genes and PPI enrichment p-value."""
    from genor.lookup import services as s

    try:
        return {"status": "ok", **s.string_network(genes, required_score=required_score)}
    except GenorError as e:
        return _err(e)


@mcp.tool
def lookup_enrichment(genes: list[str], background: list[str] | None = None,
                      sources: list[str] | None = None) -> dict[str, Any]:
    """g:Profiler functional enrichment (GO, Reactome, KEGG, HPA) with g:SCS correction."""
    from genor.lookup import services as s

    try:
        return {"status": "ok", **s.gprofiler(genes, background, sources)}
    except GenorError as e:
        return _err(e)


@mcp.tool
def report_build(run_dirs: list[str], out: str = "genor_report.md",
                 title: str = "GeNor analysis report") -> dict[str, Any]:
    """Combine run directories into a Markdown + HTML report with warnings, decisions and citations."""
    from genor.commands import report_build as rb

    ws = workspace()
    dirs = [Path(d) if Path(d).is_absolute() else ws / d for d in run_dirs]
    o = Path(out) if Path(out).is_absolute() else ws / out
    try:
        return rb(dirs, o, title)
    except GenorError as e:
        return _err(e)


def main() -> None:
    os.chdir(workspace())
    mcp.run(show_banner=False)


if __name__ == "__main__":
    main()
