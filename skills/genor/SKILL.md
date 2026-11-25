---
name: genor
description: Start here for any GeNor task. Covers GeNor's policy (GWAS summary statistics only, research use only), the session start (genor doctor), the standard post-GWAS pipeline, how to read GeNor's JSON output and exit codes, background jobs, when to stop and ask the human, and which GeNor skill handles which task (QC, meta-analysis, LDSC, MR, condFDR, fine-mapping, coloc, MAGMA, Open Targets, reports).
---

# GeNor

GeNor runs post-GWAS analyses on GWAS summary statistics: loading and QC, meta-analysis,
LDSC heritability and genetic correlation, clumping, two-sample Mendelian randomisation,
condFDR/conjFDR, ABF fine-mapping, colocalisation, MAGMA, Open Targets / Ensembl / STRING /
g:Profiler lookups, and reports. Every command writes a run directory with `manifest.json`
(inputs with SHA-256, parameters, versions, reference data, warnings, citations) and `result.json`.

## Policy (enforced in code; do not work around it)

- **Summary statistics only.** GeNor refuses individual-level data (PLINK `.bed/.pgen/.bgen`,
  `.ped`, VCFs with per-sample genotypes, BAM/CRAM) with `individual_level_refused`. Do not try
  to convert such files. Tell the user to run their GWAS elsewhere and bring the summary statistics.
- **Users bring their own GWAS.** GeNor never downloads GWAS files and never handles tokens or
  credentials (OpenGWAS JWTs, FinnGen forms, PGC terms). Point users to
  `genor sumstats sources` (MCP: `sumstats_sources`); they download the file and give you a local path.
- **Research use only.** Not a medical device; no clinical decisions, diagnosis or commercial services.
  Say so if the user frames results clinically.

## Running genor

- Claude Code plugin: `genor` is on PATH (the plugin's `bin/genor` runs it with `uvx`); the MCP
  tools are also available.
- Elsewhere: `genor` if installed (`uv tool install genor`), else `uvx genor ...`. Before the PyPI
  release use `uvx --from git+https://github.com/amanzadi/GeNor genor ...`.
- MCP-only clients (Claude Desktop, Codex desktop): use the MCP tools named below.

## Always start with

```bash
genor doctor --json        # MCP: doctor
```

Read `backends`, `refdata`, `works_without_downloads` and `blocked_until_setup`. Do not promise an
analysis that `doctor` says is blocked until the user agrees to the setup step it names.
For a walkthrough without real data: `genor example --out genor_example` (MCP: `example_data`).

## Standard pipeline order

1. Load + QC each trait: `genor sumstats load` (skill `sumstats-qc`). Lift over if builds differ.
2. Optional: meta-analyse cohorts of one trait (skill `meta-analysis`).
3. LDSC h2 per trait, then rg between traits (skill `heritability-rg`). The cross-trait intercept
   is evidence on sample overlap for step 5.
4. Clump into independent loci: `genor clump run <trait>` (default preset `loci`: r2 0.1, 250 kb).
5. Then, as the question requires: MR (skill `mendelian-randomisation`), condFDR/conjFDR
   (skill `pleiotropy-condfdr`), locus work (fine-mapping, colocalisation, gene mapping).
6. Report: `genor report build <run_dir> ... --out report.md` (skill `interpretation-and-reporting`).

Analyses take the trait **directory** written by `sumstats load` (`--out`), never the raw file.

## Reading the JSON envelope

Every command prints one JSON object (`--json`, or automatically when stdout is not a terminal):

| Field | Meaning |
|---|---|
| `status` | `ok`, `error`, or `needs_user_decision` |
| `run_dir`, `files` | where outputs are; read files from here rather than re-running |
| `summary` | the short numeric result |
| `warnings[]` | `{code, message, severity, details}`; severity is `info`, `warning` or `serious` |
| `ask_user[]` | `{question, options, why}`: the run finished, but you must put this to the human before interpreting |
| `next_steps[]` | suggested follow-up commands (check the decision points before running them) |

Exit codes: `0` ok; `2` actionable error, read `error.code` and `error.hint` and act on the hint;
`3` needs a human decision, read `error.ask_user.question` / `options`, ask, then rerun with the
answer; `1` unexpected failure (`internal_error`), rerun without `--json` for a traceback and report it.

Report every `serious` warning to the user in plain words. Never hide a warning to make a result look cleaner.

## Long runs: jobs

```bash
genor job submit -- ldsc rg ex/bmi ex/mdd --out runs/rg   # returns job_id
genor job status <job_id>     # running | done | failed | cancelled
genor job result <job_id>     # the command's JSON once finished
genor job logs <job_id>; genor job cancel <job_id>; genor job list
```

MCP tools (`sumstats_load`, `run_analysis`, `refdata_fetch`, `backend_install`) already run as jobs:
they wait up to `wait_seconds` (default 40, capped at 50) and otherwise return a `job_id`. Poll
`job_status`, then `job_result`. Give MCP tools absolute paths (relative paths resolve against
`$GENOR_WORKDIR` or the server's working directory).

## MCP equivalents

`run_analysis(method, params)`: methods `sumstats.load`, `sumstats.liftover`, `meta`, `ldsc.h2`,
`ldsc.rg`, `clump`, `mr`, `pleio`, `finemap.abf`, `coloc.abf`, `magma.gene`, `magma.sets`;
params are CLI option names with underscores (`list_methods`). Also `get_guidance(topic)`,
`read_run`, `refdata_fetch(name, confirmed_by_user, allow_unpinned)`, `backend_add(name, path)`,
`lookup_*`, `report_build`. Set `confirmed_by_user` / `allow_unpinned` only after the user agrees.

## Stop and ask the human

Never guess these. Ask, record the answer, and pass it as an option:

- **Ancestry** of each GWAS sample (`--ancestry`); it picks the LD reference.
- **Genome build** when `build_inferred` is false or `build_conflict` fires.
- **Sample size** and, for binary traits, **case and control counts** if the file lacks them.
- **Population prevalence** for liability-scale h2 of a binary trait (not the sample prevalence).
- **Exposure and outcome definitions** for MR, and the direction of the question.
- **Sample overlap** between GWAS (`--overlap none|partial|substantial|unknown`).
- **Downloads over 100 MB** (`confirm_download`): state the size, then fetch with `--yes`
  (MCP: `confirmed_by_user=true`). Unpinned datasets (`unpinned_refdata`) also need the user's OK.
- **Which loci to follow up** after clumping, condFDR or MR.
- **Installing user-licensed tools** such as MAGMA: the user downloads it, accepts the licence,
  and registers it with `genor backend add magma <path>`.
- Relaxing thresholds (MR instrument p-value, FDR levels) or using approximate fallbacks
  (distance clumping).

## Which skill for which task

| Task | Skill |
|---|---|
| Load a file, column mapping, OR vs beta, build, liftover, QC report | `sumstats-qc` |
| Combine cohorts of one trait (IVW, sample-size, heterogeneity) | `meta-analysis` |
| SNP heritability, genetic correlation, intercepts | `heritability-rg` |
| Causal effect of exposure on outcome, STROBE-MR | `mendelian-randomisation` |
| Shared or conditional loci across two traits (condFDR / conjFDR) | `pleiotropy-condfdr` |
| Credible sets in a locus (ABF) | `fine-mapping` |
| Shared causal variant between two traits in a region | `colocalisation` |
| Gene, gene-set and cell-type analysis (MAGMA), STRING, enrichment | `gene-pathway-celltype` |
| Variant-to-gene and target evidence from Open Targets | `open-targets-triage` |
| Writing up results, building the report, wording limits | `interpretation-and-reporting` |
| Background on GWAS statistics and terms | `gwas-concepts` |

Not in GeNor yet: SuSiE fine-mapping, COJO runs, SMR/HEIDI runs, MR-PRESSO, multivariable MR,
LD proxies for missing instruments, partitioned/stratified LDSC, cross-ancestry rg (Popcorn),
MiXeR. Say so instead of improvising a substitute.
