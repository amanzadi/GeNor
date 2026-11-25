---
name: interpretation-and-reporting
description: Build and write up GeNor results - genor report build (Markdown + HTML from run manifests with warnings, human decisions and citations), how to state effect sizes with confidence intervals and scales, multiple testing, causal language, reproducibility via manifest.json and SHA-256, STROBE-MR, and the research-use notice. Use when the user asks for a report, summary, write-up, methods or results text, figure legends, or asks what a result means and how strongly it can be claimed.
---
# Interpretation and reporting

Every GeNor analysis writes a run directory with `manifest.json` (inputs with SHA-256, parameters,
versions, reference data, backends, warnings, decisions for the human, outputs, citations) and
`result.json`. `genor report build` turns one or more run directories into `report.md` and
`report.html`. The agent's job is to add the interpretation the tool cannot: correct scales,
caveats, and wording that does not overstate.

## When to use

- At the end of any workflow, or when the user asks "what does this mean?", "write the methods",
  "summarise the results", or "make a report".
- Before sharing results: check that every warning and every `ask_user` entry is handled.

## Inputs and prerequisites

- Run directories from GeNor commands (each has `manifest.json`). List them with `genor runs`
  (MCP: `list_runs`); inspect one with MCP `read_run(run_dir)` or by reading `result.json`.
- Answers from the human to every `ask_user` question raised by the runs.
- For lookups (Open Targets, VEP, g:Profiler, STRING): the outputs you saved and the query date.
  Lookups do not create run directories.

## Steps

1. Collect run directories: `genor runs --json`.
2. Review each run's `warnings` (severity `info`, `warning`, `serious`) and `ask_user`. Resolve
   open decisions with the user before interpreting. Never answer them yourself.
3. Build the report. Put the output file in a directory that contains the run directories so
   figure links are relative and the HTML is portable:
   ```bash
   genor report build runs/bmi runs/clump_bmi runs/fm_bmi_rs796513 runs/coloc_bmi_mdd_rs2237324 \
     --out runs/report.md --title "BMI and MDD: post-GWAS analysis" --json
   ```
   MCP: `report_build(run_dirs=[...], out="runs/report.md", title="...")`.
   Output: `markdown`, `html`, `n_runs`, `n_citations`.
4. Add an interpretation section in your own words (the report contains numbers, not
   conclusions). Use the templates in
   [references/writing-results.md](references/writing-results.md).
5. Add lookup results and their citations by hand (see Reading the output).
6. For MR runs, point the user to the STROBE-MR checklist in the report and list the items marked
   **HUMAN** that they must complete (see `mendelian-randomisation`).
7. Give the user the paths to `report.md` and `report.html`, and the list of open caveats.

## Defaults and why

- `--out genor_report.md`, `--title "GeNor analysis report"`. HTML is written next to the
  Markdown with the same stem.
- Every report starts with GeNor's research-use notice: "GeNor is research software. It is not a
  medical device and must not be used for clinical decisions, diagnosis, or commercial services."
  Keep it.
- Citations are deduplicated across runs, from the keys each run recorded (GeNor itself, plus the
  methods and data used, e.g. LDSC, coloc, Wakefield ABF, MAGMA, PLINK 2, 1000 Genomes,
  STROBE-MR).
- Reporting standards: STROBE-MR for MR (Skrivankova et al. 2021); for GWAS-derived analyses in
  general, the STREGA extension of STROBE is a useful checklist.

## Decision points

- **Ask the user which runs belong in the report** (exploratory and failed runs are included if
  passed; a failed run shows its error).
- **Ask the user to answer every open `ask_user`** (build, ancestry, N, population prevalence,
  sample overlap). Report unanswered ones as open limitations.
- **Ask the user for the facts only they know**: phenotype definitions, cohorts, consent and data
  access terms, funding, pre-registration, conflicts of interest.
- **Ask the user how to handle multiple testing** across the analyses they ran (see
  [references/writing-results.md](references/writing-results.md)).

## Reading the output

Each run becomes a section (layout in [references/build-output-layout.md](references/build-output-layout.md)):
heading with the command, run directory, status, GeNor version and timestamp; inputs with the
first 12 characters of the SHA-256; reference data and tools with versions; parameters; a result
table (flattened summary, first 60 rows); warnings as `[severity] code: message`; "Decisions for
the human"; figures (PNG outputs); the STROBE-MR table for MR runs; the error if the run failed.
A "References" list ends the report.

Not in the report automatically:
- Lookup results and their citations. Add them, e.g. "Open Targets Platform (queried 2026-11-25;
  Buniello et al. 2025)", "Ensembl VEP (McLaren et al. 2016)", "g:Profiler version from
  `gprofiler_version` (Kolberg et al. 2023)", "STRING v12.0, score >= 700 (Szklarczyk et al. 2023)".
- Interpretation, limitations, and anything requiring judgement.

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `not_a_run_dir` | Directory has no `manifest.json` | Pass directories written by GeNor commands (a trait directory from `sumstats load` is a run directory). |
| Broken images in HTML | Report written outside the run directories' parent, so image paths are absolute | Use `--out` inside the common parent directory. |
| Run section shows status `error` | That run failed; its error code is shown | Fix and rerun, or explain in the text why it is included. |

## What not to conclude

- Do not write "X causes Y" from association, colocalisation, or a single MR estimate. Use the
  wording table in [references/writing-results.md](references/writing-results.md).
- Do not present posterior probabilities (fine-mapping PP, coloc PP.H4) as p-values, or p-values
  as the probability that a hypothesis is true.
- Do not drop warnings because the main result looks good. A `serious` warning must appear in the
  limitations.
- Do not report simulated example data (BMI_sim, MDD_sim) as biology.
- Nothing in a GeNor report supports clinical decisions about individuals.

## Related skills

- `genor` (router), `mendelian-randomisation` (STROBE-MR), `heritability-rg`, `fine-mapping`,
  `colocalisation`, `gene-pathway-celltype`, `open-targets-triage`, `pleiotropy-condfdr`,
  `meta-analysis`, `gwas-concepts` (explaining terms to the user).
- References: [build-output-layout.md](references/build-output-layout.md),
  [writing-results.md](references/writing-results.md),
  [reproducibility.md](references/reproducibility.md).
