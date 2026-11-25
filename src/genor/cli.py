# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""GeNor command-line interface.

Conventions (for agents):
- Output is one JSON object on stdout when ``--json`` is given (anywhere on the line)
  or when stdout is not a terminal. Bulk results go to files in the run directory.
- Exit codes: 0 ok, 2 actionable error (see ``error.code`` / ``error.hint``),
  3 a decision for the human is needed (``error.ask_user``), 1 unexpected failure.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Annotated, Any

import typer

from genor import RESEARCH_USE_NOTICE, __version__
from genor.errors import GenorError, NeedsUserDecision

JSON_MODE = False

app = typer.Typer(
    name="genor",
    help="GeNor: post-GWAS analysis from summary statistics, built for AI agents.\n\n"
         + RESEARCH_USE_NOTICE,
    no_args_is_help=True,
    add_completion=False,
    pretty_exceptions_enable=False,
)
sumstats_app = typer.Typer(help="Load, QC, lift over summary statistics; where to find GWAS.", no_args_is_help=True)
meta_app = typer.Typer(help="Meta-analysis (IVW fixed/random effects, sample-size weighted).", no_args_is_help=True)
ldsc_app = typer.Typer(help="LD score regression: SNP heritability and genetic correlation.", no_args_is_help=True)
clump_app = typer.Typer(help="LD clumping into independent loci.", no_args_is_help=True)
mr_app = typer.Typer(help="Two-sample Mendelian randomisation with STROBE-MR report.", no_args_is_help=True)
pleio_app = typer.Typer(help="condFDR / conjFDR with pleioFDR (Python port).", no_args_is_help=True)
finemap_app = typer.Typer(help="Fine-mapping (ABF, single causal variant).", no_args_is_help=True)
coloc_app = typer.Typer(help="Colocalisation (coloc ABF).", no_args_is_help=True)
magma_app = typer.Typer(help="MAGMA gene / gene-set / gene-property analysis (user-installed MAGMA).", no_args_is_help=True)
lookup_app = typer.Typer(help="Open Targets, Ensembl, STRING and g:Profiler lookups.", no_args_is_help=True)
report_app = typer.Typer(help="Build reports from run directories.", no_args_is_help=True)
refdata_app = typer.Typer(help="Reference data (LD scores, LD panels, chain files).", no_args_is_help=True)
backend_app = typer.Typer(help="External tools (PLINK 2, GCTA, MAGMA, pleioFDR, R).", no_args_is_help=True)
job_app = typer.Typer(help="Background jobs for long analyses.", no_args_is_help=True)
skills_app = typer.Typer(help="Agent skills: list, show, install for Claude/Codex/opencode.", no_args_is_help=True)

for sub, name in [(sumstats_app, "sumstats"), (meta_app, "meta"), (ldsc_app, "ldsc"),
                  (clump_app, "clump"), (mr_app, "mr"), (pleio_app, "pleio"),
                  (finemap_app, "finemap"), (coloc_app, "coloc"), (magma_app, "magma"),
                  (lookup_app, "lookup"), (report_app, "report"), (refdata_app, "refdata"),
                  (backend_app, "backend"), (job_app, "job"), (skills_app, "skills")]:
    app.add_typer(sub, name=name)


# ---------------------------------------------------------------- output


def emit(result: Any) -> None:
    if JSON_MODE:
        sys.stdout.write(json.dumps(result, indent=1, default=str) + "\n")
        return
    from rich.console import Console
    from rich.pretty import Pretty

    c = Console()
    if isinstance(result, dict):
        status = result.get("status")
        if status:
            c.print(f"[bold]{result.get('command', 'genor')}[/bold]: {status}")
        if result.get("run_dir"):
            c.print(f"run dir: {result['run_dir']}")
        body = {k: v for k, v in result.items() if k not in ("status", "command", "run_dir")}
        warnings = body.pop("warnings", None)
        asks = body.pop("ask_user", None)
        c.print(Pretty(body, max_length=20, max_string=300))
        if warnings:
            c.print("[yellow]warnings:[/yellow]")
            for w in warnings:
                c.print(f"  [{w['severity']}] {w['code']}: {w['message']}")
        if asks:
            c.print("[magenta]decisions for the human:[/magenta]")
            for a in asks:
                c.print(f"  ? {a['question']}")
    else:
        c.print(Pretty(result))


def _fail(err: GenorError) -> None:
    out = {"status": "needs_user_decision" if isinstance(err, NeedsUserDecision) else "error",
           "error": err.to_dict()}
    if JSON_MODE:
        sys.stdout.write(json.dumps(out, indent=1, default=str) + "\n")
    else:
        sys.stderr.write(f"error [{err.code}]: {err.message}\n")
        if err.hint:
            sys.stderr.write(f"hint: {err.hint}\n")
        if isinstance(err, NeedsUserDecision) and err.options:
            sys.stderr.write(f"options: {', '.join(err.options)}\n")
    raise SystemExit(3 if isinstance(err, NeedsUserDecision) else 2)


def _out(out: Path | None, default: str) -> Path:
    return out or Path("genor_runs") / default


def _kv(items: list[str] | None) -> dict[str, str] | None:
    if not items:
        return None
    out = {}
    for it in items:
        if "=" not in it:
            raise GenorError(f"--col expects canonical=FILECOLUMN, got '{it}'", code="bad_param")
        k, v = it.split("=", 1)
        out[k.strip()] = v.strip()
    return out


Out = Annotated[Path | None, typer.Option("--out", "-o", help="Run directory to write (created).")]

# ---------------------------------------------------------------- top level


@app.command()
def doctor() -> None:
    """Check platform, backends and reference data; say what each missing piece blocks."""
    from genor.commands import doctor as d

    emit(d())


@app.command()
def version() -> None:
    """Print the GeNor version."""
    emit({"genor": __version__})


@app.command()
def example(out: Annotated[Path, typer.Option("--out", "-o")] = Path("genor_example")) -> None:
    """Copy the bundled simulated GWAS pair (BMI_sim, MDD_sim) for tests and walkthroughs."""
    from genor.commands import example as ex

    emit(ex(out))


@app.command()
def runs() -> None:
    """List run directories under the current directory."""
    from genor.commands import list_runs

    emit(list_runs())


@app.command()
def mcp() -> None:
    """Start the MCP server on stdio (for Claude Desktop, Claude Code, Codex, opencode)."""
    from genor.mcp_server import main as m

    m()


# ---------------------------------------------------------------- sumstats


@sumstats_app.command("sources")
def sumstats_sources(name: Annotated[str | None, typer.Argument()] = None) -> None:
    """Where to find GWAS summary statistics (GeNor never downloads them itself)."""
    from genor.commands import sources

    emit(sources(name))


@sumstats_app.command("load")
def sumstats_load(
    file: Path,
    trait: Annotated[str, typer.Option(help="Short trait name, e.g. BMI.")],
    out: Out = None,
    build: Annotated[str | None, typer.Option(help="GRCh37 or GRCh38 (inferred if omitted).")] = None,
    ancestry: Annotated[str | None, typer.Option(help="EUR, EAS, AFR, AMR, SAS, MID, CSA, or multi.")] = None,
    n: Annotated[float | None, typer.Option(help="Total sample size if not in the file.")] = None,
    n_cases: Annotated[float | None, typer.Option(help="Number of cases (binary traits).")] = None,
    n_controls: Annotated[float | None, typer.Option(help="Number of controls (binary traits).")] = None,
    population_prevalence: Annotated[float | None, typer.Option(help="For liability-scale h2.")] = None,
    effect_type: Annotated[str | None, typer.Option(help="beta or or (odds ratio).")] = None,
    col: Annotated[list[str] | None, typer.Option(help="Column override canonical=FILECOL, repeatable.")] = None,
    maf: Annotated[float, typer.Option(help="Minimum MAF.")] = 0.01,
    info: Annotated[float, typer.Option(help="Minimum imputation INFO.")] = 0.8,
    drop_mhc: Annotated[bool, typer.Option(help="Remove chr6:25-34 Mb (default: keep and flag).")] = False,
    plots: Annotated[bool, typer.Option(help="Write Manhattan and QQ plots.")] = True,
) -> None:
    """Parse a summary statistics file (any common format, GWAS-VCF), standardise, QC, plot."""
    from genor.commands import sumstats_load as f

    emit(f(file, trait, _out(out, trait), build=build, ancestry=ancestry, n=n, n_cases=n_cases,
           n_controls=n_controls, population_prevalence=population_prevalence,
           effect_type=effect_type, columns=_kv(col), maf=maf, info=info, drop_mhc=drop_mhc,
           plots=plots))


@sumstats_app.command("liftover")
def sumstats_liftover(trait: Path, to: Annotated[str, typer.Option(help="GRCh37 or GRCh38")], out: Out = None) -> None:
    """Lift a loaded trait to another genome build (UCSC chain files)."""
    from genor.commands import sumstats_liftover as f

    emit(f(trait, to, _out(out, f"{trait.name}_{to}")))


# ---------------------------------------------------------------- meta


@meta_app.command("run")
def meta_run(traits: list[Path], out: Out = None,
             scheme: Annotated[str, typer.Option(help="ivw or samplesize")] = "ivw",
             name: Annotated[str, typer.Option(help="Trait name for the result.")] = "meta") -> None:
    """Meta-analyse two or more loaded traits."""
    from genor.commands import meta_run as f

    emit(f(traits, _out(out, f"meta_{name}"), scheme=scheme, trait_name=name))


# ---------------------------------------------------------------- ldsc


@ldsc_app.command("h2")
def ldsc_h2(trait: Path, out: Out = None,
            ref: Annotated[str, typer.Option(help="refdata name of LD scores")] = "ldsc-ldscores-eur",
            w: Annotated[str, typer.Option(help="refdata name of regression weights")] = "ldsc-weights-eur",
            pop_prev: Annotated[float | None, typer.Option(help="Population prevalence (binary traits).")] = None,
            intercept: Annotated[float | None, typer.Option(help="Constrain the intercept.")] = None) -> None:
    """SNP heritability with LD score regression."""
    from genor.commands import ldsc_h2 as f

    emit(f(trait, _out(out, f"h2_{trait.name}"), ref=ref, w=w, pop_prev=pop_prev, intercept=intercept))


@ldsc_app.command("rg")
def ldsc_rg(trait1: Path, trait2: Path, out: Out = None,
            ref: Annotated[str, typer.Option()] = "ldsc-ldscores-eur",
            w: Annotated[str, typer.Option()] = "ldsc-weights-eur",
            pop_prev1: Annotated[float | None, typer.Option(help="Population prevalence of trait1 (binary)")] = None,
            pop_prev2: Annotated[float | None, typer.Option(help="Population prevalence of trait2 (binary)")] = None) -> None:
    """Genetic correlation between two traits with LD score regression."""
    from genor.commands import ldsc_rg as f

    emit(f(trait1, trait2, _out(out, f"rg_{trait1.name}_{trait2.name}"), ref=ref, w=w,
           pop_prev1=pop_prev1, pop_prev2=pop_prev2))


# ---------------------------------------------------------------- clump


@clump_app.command("run")
def clump_run(trait: Path, out: Out = None,
              mode: Annotated[str, typer.Option(help="plink (LD panel) or distance (no LD; approximate)")] = "plink",
              preset: Annotated[str, typer.Option(help="loci (r2 0.1, 250 kb) or mr (r2 0.001, 10 Mb)")] = "loci",
              p1: Annotated[float | None, typer.Option(help="Index p-value threshold (preset default 5e-8)")] = None,
              r2: Annotated[float | None, typer.Option(help="LD r2 threshold (preset default)")] = None,
              kb: Annotated[int | None, typer.Option(help="Window in kb (preset default)")] = None,
              panel: Annotated[str | None, typer.Option(help="refdata name, e.g. g1000-eur")] = None) -> None:
    """Clump a trait into independent index variants."""
    from genor.commands import clump_run as f

    emit(f(trait, _out(out, f"clump_{trait.name}"), mode=mode, preset=preset, p1=p1, r2=r2, kb=kb, panel=panel))


# ---------------------------------------------------------------- MR


@mr_app.command("run")
def mr_run(
    exposure: Annotated[Path, typer.Option(help="Loaded exposure trait directory.")],
    outcome: Annotated[Path, typer.Option(help="Loaded outcome trait directory.")],
    out: Out = None,
    p_threshold: Annotated[float, typer.Option(help="Instrument p-value threshold.")] = 5e-8,
    clump_mode: Annotated[str, typer.Option(help="plink or distance")] = "plink",
    panel: Annotated[str | None, typer.Option(help="LD panel refdata name (default from ancestry, e.g. g1000-eur)")] = None,
    clump_r2: Annotated[float, typer.Option(help="Clumping r2 for instruments")] = 0.001,
    clump_kb: Annotated[int, typer.Option(help="Clumping window (kb)")] = 10_000,
    instruments: Annotated[Path | None, typer.Option(help="File whose first column lists instrument rsIDs.")] = None,
    exclude_mhc: Annotated[bool, typer.Option(help="Drop instruments in chr6:25-34 Mb (default: keep)")] = False,
    overlap: Annotated[str | None, typer.Option(help="Sample overlap stated by the user: none|partial|substantial|unknown")] = None,
    nboot: Annotated[int, typer.Option(help="Bootstrap replicates (weighted median/mode SE)")] = 1000,
    seed: Annotated[int, typer.Option(help="Random seed")] = 1,
) -> None:
    """Instrument selection, harmonisation, IVW/Egger/median/mode, diagnostics, plots, STROBE-MR."""
    from genor.commands import mr_run as f

    emit(f(exposure, outcome, _out(out, f"mr_{exposure.name}_{outcome.name}"), p_threshold=p_threshold,
           clump_mode=clump_mode, panel=panel, clump_r2=clump_r2, clump_kb=clump_kb,
           instruments=instruments, exclude_mhc=exclude_mhc, overlap=overlap, nboot=nboot, seed=seed))


# ---------------------------------------------------------------- pleio


@pleio_app.command("run")
def pleio_run(trait1: Path, trait2: Path, out: Out = None,
              stat: Annotated[str, typer.Option(help="condfdr or conjfdr")] = "condfdr",
              fdr: Annotated[float | None, typer.Option(help="Threshold (default 0.01 condFDR, 0.05 conjFDR)")] = None,
              randprune_n: Annotated[int, typer.Option()] = 500,
              seed: Annotated[int, typer.Option()] = 1,
              keep_mhc: Annotated[bool, typer.Option(help="Do not exclude the MHC")] = False,
              prepare_only: Annotated[bool, typer.Option(help="Write MAT files and config, don't run")] = False,
              ref_mat: Annotated[Path | None, typer.Option(help="Existing pleioFDR reference .mat (instead of refdata)")] = None,
              template: Annotated[Path | None, typer.Option(help="Existing 9545380.ref SNP template (instead of refdata)")] = None) -> None:
    """condFDR (trait1 conditional on trait2) or conjFDR (shared loci)."""
    from genor.commands import pleio_run as f

    emit(f(trait1, trait2, _out(out, f"{stat}_{trait1.name}_{trait2.name}"), stattype=stat, fdr=fdr,
           randprune_n=randprune_n, exclude_mhc=not keep_mhc, seed=seed, prepare_only=prepare_only,
           ref_mat=ref_mat, template=template))


# ---------------------------------------------------------------- locus level


@finemap_app.command("abf")
def finemap_abf(trait: Path, out: Out = None,
                region: Annotated[str | None, typer.Option(help="chr:start-end")] = None,
                lead: Annotated[str | None, typer.Option(help="Lead variant ID; region = lead +/- window")] = None,
                window_kb: Annotated[int, typer.Option()] = 500,
                p1: Annotated[float, typer.Option(help="Prior per variant")] = 1e-4,
                level: Annotated[float, typer.Option(help="Credible set level")] = 0.95) -> None:
    """Wakefield ABF fine-mapping with a credible set."""
    from genor.commands import finemap_abf as f

    emit(f(trait, _out(out, f"finemap_{trait.name}_{lead or region}".replace(":", "_")), region=region,
           lead=lead, window_kb=window_kb, p1=p1, level=level))


@coloc_app.command("abf")
def coloc_abf(trait1: Path, trait2: Path, out: Out = None,
              region: Annotated[str | None, typer.Option()] = None,
              lead: Annotated[str | None, typer.Option()] = None,
              window_kb: Annotated[int, typer.Option()] = 500,
              p1: Annotated[float, typer.Option()] = 1e-4,
              p2: Annotated[float, typer.Option()] = 1e-4,
              p12: Annotated[float, typer.Option()] = 1e-5) -> None:
    """Colocalisation of two traits in a region (coloc.abf)."""
    from genor.commands import coloc_run as f

    emit(f(trait1, trait2, _out(out, f"coloc_{trait1.name}_{trait2.name}"), region=region, lead=lead,
           window_kb=window_kb, p1=p1, p2=p2, p12=p12))


# ---------------------------------------------------------------- MAGMA


@magma_app.command("gene")
def magma_gene(trait: Path, gene_loc: Annotated[Path, typer.Option(help="MAGMA gene location file")],
               ref_bfile: Annotated[str, typer.Option(help="MAGMA reference prefix, e.g. g1000_eur")],
               out: Out = None, up_kb: int = 35, down_kb: int = 10) -> None:
    """MAGMA gene-based analysis."""
    from genor.commands import magma_gene as f

    emit(f(trait, _out(out, f"magma_{trait.name}"), gene_loc=gene_loc, ref_bfile=ref_bfile,
           window_up_kb=up_kb, window_down_kb=down_kb))


@magma_app.command("sets")
def magma_sets(genes_raw: Path, out: Out = None,
               set_annot: Annotated[Path | None, typer.Option(help="Gene-set file (e.g. MSigDB)")] = None,
               covar: Annotated[Path | None, typer.Option(help="Gene-property file (e.g. tissue/cell-type expression)")] = None) -> None:
    """MAGMA gene-set or gene-property (cell type / tissue) analysis."""
    from genor.commands import magma_sets as f

    emit(f(genes_raw, _out(out, "magma_sets"), set_annot=set_annot, covar=covar))


# ---------------------------------------------------------------- lookups


@lookup_app.command("search")
def lookup_search(term: str, entity: Annotated[str | None, typer.Option(help="target | disease | variant | study")] = None,
                  n: Annotated[int, typer.Option(help="Number of hits")] = 5) -> None:
    """Search Open Targets for a gene, disease, variant or study ID."""
    from genor.lookup import services as s

    emit({"status": "ok", "hits": s.ot_search(term, entity, n)})


@lookup_app.command("target")
def lookup_target(ensembl_id: str, n: int = 10) -> None:
    """Open Targets target summary with top associated diseases and tractability."""
    from genor.lookup import services as s

    emit({"status": "ok", "target": s.ot_target(ensembl_id, n)})


@lookup_app.command("disease")
def lookup_disease(efo_id: str, n: int = 25) -> None:
    """Open Targets top targets for a disease (EFO/MONDO ID)."""
    from genor.lookup import services as s

    emit({"status": "ok", "disease": s.ot_disease_targets(efo_id, n)})


@lookup_app.command("variant")
def lookup_variant(variant_id: Annotated[str, typer.Argument(help="GRCh38 CHR_POS_REF_ALT")], n: int = 10) -> None:
    """Open Targets credible sets and L2G genes that include a variant."""
    from genor.lookup import services as s

    emit({"status": "ok", "variant": s.ot_variant_credible_sets(variant_id, n)})


@lookup_app.command("coloc")
def lookup_coloc(study_locus_id: str, n: int = 20) -> None:
    """Open Targets colocalisations for a credible set (studyLocusId)."""
    from genor.lookup import services as s

    emit({"status": "ok", "credible_set": s.ot_study_locus_coloc(study_locus_id, n)})


@lookup_app.command("vep")
def lookup_vep(rsids: list[str], build: str = "GRCh38") -> None:
    """Ensembl VEP consequences and genes for up to 200 rsIDs."""
    from genor.lookup import services as s

    emit({"status": "ok", "variants": s.vep(rsids, build)})


@lookup_app.command("genes")
def lookup_genes(region: Annotated[str, typer.Option(help="chr:start-end")],
                 build: Annotated[str, typer.Option()] = "GRCh37") -> None:
    """Genes overlapping a region (Ensembl)."""
    from genor.lookup import services as s
    from genor.methods.coloc import parse_region

    c, a, b = parse_region(region)
    emit({"status": "ok", "build": build, "genes": s.genes_in_region(c, a, b, build)})


@lookup_app.command("string")
def lookup_string(genes: list[str], required_score: int = 700) -> None:
    """STRING v12 protein interaction network and PPI enrichment for a gene list."""
    from genor.lookup import services as s

    emit({"status": "ok", **s.string_network(genes, required_score=required_score)})


@lookup_app.command("enrich")
def lookup_enrich(genes: list[str],
                  background: Annotated[Path | None, typer.Option(help="File with background genes, one per line")] = None,
                  sources: Annotated[str | None, typer.Option(help="Comma-separated, e.g. GO:BP,REAC")] = None) -> None:
    """g:Profiler enrichment (GO, Reactome, KEGG, HPA)."""
    from genor.lookup import services as s

    bg = background.read_text().split() if background else None
    emit({"status": "ok", **s.gprofiler(genes, bg, sources.split(",") if sources else None)})


# ---------------------------------------------------------------- report


@report_app.command("build")
def report_build(run_dirs: list[Path],
                 out: Annotated[Path, typer.Option("--out", "-o")] = Path("genor_report.md"),
                 title: str = "GeNor analysis report") -> None:
    """Combine run directories into report.md + report.html with all citations and warnings."""
    from genor.commands import report_build as f

    emit(f(run_dirs, out, title))


# ---------------------------------------------------------------- refdata


@refdata_app.command("list")
def refdata_list() -> None:
    """Reference datasets, sizes and whether they are downloaded."""
    from genor import refdata

    emit({"status": "ok", "datasets": refdata.status()})


@refdata_app.command("fetch")
def refdata_fetch(names: list[str],
                  yes: Annotated[bool, typer.Option("--yes", "-y", help="Confirm large downloads")] = False,
                  allow_unpinned: Annotated[bool, typer.Option(help="Allow datasets without a pinned checksum")] = False) -> None:
    """Download reference datasets (checksum-verified)."""
    from genor import refdata

    emit({"status": "ok", "fetched": [refdata.fetch(n, yes=yes, allow_unpinned=allow_unpinned) for n in names]})


@refdata_app.command("remove")
def refdata_remove(name: str) -> None:
    """Delete a downloaded dataset."""
    from genor import refdata

    refdata.remove(name)
    emit({"status": "ok", "removed": name})


# ---------------------------------------------------------------- backend


@backend_app.command("list")
def backend_list() -> None:
    """External tools and how to enable the missing ones."""
    from genor import backends

    emit({"status": "ok", "platform": backends.platform_key(), "backends": backends.status()})


@backend_app.command("install")
def backend_install(name: str) -> None:
    """Download a pinned, redistributable binary (plink2, gcta, smr)."""
    from genor import backends

    emit({"status": "ok", **backends.install(name)})


@backend_app.command("add")
def backend_add(name: str, path: Path) -> None:
    """Register a tool you installed yourself (e.g. MAGMA)."""
    from genor import backends

    emit({"status": "ok", **backends.add(name, path)})


@backend_app.command("remove")
def backend_remove(name: str) -> None:
    """Forget a registered tool."""
    from genor import backends

    backends.remove(name)
    emit({"status": "ok", "removed": name})


# ---------------------------------------------------------------- jobs


@job_app.command("submit", context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def job_submit(ctx: typer.Context) -> None:
    """Run any genor command in the background: genor job submit -- ldsc rg A B --out R"""
    from genor import jobs

    args = [a for a in ctx.args if a != "--"]
    if not args:
        raise GenorError("Nothing to run.", code="bad_param", hint="genor job submit -- <command> ...")
    emit({"status": "ok", **jobs.submit(args)})


@job_app.command("status")
def job_status(job_id: str) -> None:
    """Job state (running / done / failed / cancelled)."""
    from genor import jobs

    emit({"status": "ok", "job": jobs.status(job_id)})


@job_app.command("result")
def job_result(job_id: str) -> None:
    """The finished job's result JSON."""
    from genor import jobs

    emit(jobs.result(job_id))


@job_app.command("logs")
def job_logs(job_id: str, n: int = 50) -> None:
    """Last lines of the job's stdout/stderr."""
    from genor import jobs

    emit({"status": "ok", **jobs.logs(job_id, n)})


@job_app.command("cancel")
def job_cancel(job_id: str) -> None:
    """Stop a running job."""
    from genor import jobs

    emit({"status": "ok", **jobs.cancel(job_id)})


@job_app.command("list")
def job_list(limit: int = 20) -> None:
    """Recent jobs."""
    from genor import jobs

    emit({"status": "ok", "jobs": jobs.list_jobs(limit)})


# ---------------------------------------------------------------- skills


@skills_app.command("list")
def skills_list() -> None:
    """Skills shipped with GeNor."""
    from genor import skills

    emit({"status": "ok", "skills": skills.list_skills()})


@skills_app.command("show")
def skills_show(topic: Annotated[str | None, typer.Argument()] = None) -> None:
    """Print a skill (or a reference file: <skill>/references/<file>.md)."""
    from genor import skills

    r = skills.guidance(topic)
    if not JSON_MODE and "text" in r:
        sys.stdout.write(r["text"])
        return
    emit(r)


@skills_app.command("install")
def skills_install(agent: Annotated[str, typer.Option(help="claude | codex | opencode")] = "claude",
                   scope: Annotated[str, typer.Option(help="project | user")] = "project",
                   dest: Annotated[Path | None, typer.Option(help="Explicit target directory")] = None,
                   force: bool = False) -> None:
    """Copy GeNor's skills to .claude/skills or .agents/skills."""
    from genor import skills

    emit(skills.install(agent, scope, dest, force))


# ---------------------------------------------------------------- entry point


def main() -> None:
    global JSON_MODE
    argv = sys.argv[1:]
    if "--json" in argv:
        argv = [a for a in argv if a != "--json"]
        JSON_MODE = True
    elif os.environ.get("GENOR_JSON") == "1" or not sys.stdout.isatty():
        JSON_MODE = True
    if argv and argv[0] == "mcp":
        JSON_MODE = False
    try:
        app(args=argv, prog_name="genor", standalone_mode=False)
    except GenorError as e:
        _fail(e)
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001
        import click

        if isinstance(e, click.exceptions.ClickException):
            e.show()
            raise SystemExit(e.exit_code) from e
        if isinstance(e, click.exceptions.Abort):
            raise SystemExit(1) from e
        if JSON_MODE:
            sys.stdout.write(json.dumps({"status": "error", "error": {
                "code": "internal_error", "message": f"{type(e).__name__}: {e}",
                "hint": "This is a GeNor bug or an unexpected input. Rerun without --json for a traceback."}}) + "\n")
            raise SystemExit(1) from e
        raise


if __name__ == "__main__":
    main()
