# GeNor

**Post-GWAS analysis from summary statistics, built for AI agents.**

GeNor lets Claude (and Codex, opencode, or you at a terminal) take GWAS summary statistics from
QC to genetic correlation, Mendelian randomisation, pleiotropy, fine-mapping, colocalisation and
gene-level follow-up, with a reproducible record of every step. It installs with nothing but
[uv](https://docs.astral.sh/uv/).

> [!IMPORTANT]
> GeNor is research software for academic use. It is not a medical device and must not be used
> for clinical decisions, diagnosis, or commercial services.

- **Summary statistics only.** GeNor refuses individual-level genotype data. You bring your own
  GWAS files; GeNor tells you where to find them (`genor sumstats sources`) but never downloads
  GWAS or handles tokens.
- **Agent-first.** Every command prints one JSON object with a short summary, structured
  warnings, and the questions the agent must put to you (ancestry, sample overlap, large
  downloads...). Bulk results go to files, not into the agent's context.
- **Guard rails in code.** Build and ancestry mismatches against an LD panel are errors;
  palindromic SNPs, allele flips, MHC, inflation and sample overlap are checked and reported.
- **Reproducible.** Each run writes `manifest.json`: input checksums, parameters, tool and data
  versions, warnings, and the citations to use.
- **Validated.** Native reimplementations are tested against the reference tools
  (TwoSampleMR, coloc, R `meta`, original LDSC). See [validation](docs/validation.md).

## Quick start

### Claude Code (recommended)

```text
/plugin marketplace add amanzadi/Genor
/plugin install genor@genor
```

This installs the GeNor skills, the MCP server, slash commands (`/genor:pipeline`, `/genor:mr`,
`/genor:rg`, `/genor:doctor`, `/genor:example`) and puts `genor` on Claude's PATH. Then ask, in a
folder with your files:

> Estimate the genetic correlation and causal effect of BMI on depression from `bmi.tsv.gz`
> and `mdd.tsv.gz`.

Claude checks the setup, asks you for what it cannot know (ancestry, case/control counts,
population prevalence, sample overlap, permission for downloads), runs the analyses and writes
a report.

### Claude Desktop, Codex, opencode

All use the same MCP server (`genor mcp`). Copy-paste configs: [docs/install.md](docs/install.md).

### Command line

```bash
uvx --from git+https://github.com/amanzadi/Genor genor doctor      # until the PyPI release
# after release: uvx genor doctor   or   uv tool install genor
```

## A first run on the bundled example

GeNor ships a simulated GWAS pair (not real data): BMI_sim (quantitative, N = 300,000) and
MDD_sim (60,000 cases / 140,000 controls) on ~98,000 HapMap3 SNPs, with a planted causal
effect of 0.15 and a genetic correlation of 0.15.

```bash
genor example --out ex
genor sumstats load ex/bmi_sim.tsv.gz --trait BMI_sim --ancestry EUR --out runs/bmi
genor sumstats load ex/mdd_sim.tsv.gz --trait MDD_sim --ancestry EUR \
      --population-prevalence 0.15 --out runs/mdd
genor refdata fetch ldsc-ldscores-eur ldsc-weights-eur          # 25 MB
genor ldsc rg runs/bmi runs/mdd --out runs/rg
genor backend install plink2 && genor refdata fetch g1000-eur --yes   # 288 MB LD panel
genor mr run --exposure runs/bmi --outcome runs/mdd --overlap none --out runs/mr
genor report build runs/bmi runs/mdd runs/rg runs/mr --out report.md
```

What comes back (truth in parentheses):

| | Estimate |
|---|---|
| LDSC rg | 0.142 (SE 0.032) (0.15) |
| MDD h2, liability scale | 0.092 |
| MR, IVW (141 instruments) | 0.136 (SE 0.018), OR 1.15 per unit (0.15) |
| Weighted median / weighted mode | 0.147 / 0.152 |
| Flags raised | heterogeneity, directional pleiotropy (3 pleiotropic instruments were planted) |

The MR run also writes scatter, forest, funnel and leave-one-out plots and a
[STROBE-MR](https://doi.org/10.1001/jama.2021.18236) checklist with the items GeNor can fill
from the run, the ones it can partly fill, and the ones only the authors can write.

## What is included

| Area | Method | How |
|---|---|---|
| Sumstats | load any common format (incl. GWAS-VCF, PLINK 2, REGENIE, SAIGE, METAL, FinnGen), standardise to GWAS-SSF columns, infer build, QC, liftover, Manhattan/QQ | native |
| Meta-analysis | fixed-effect IVW, DerSimonian-Laird random effects, Q, I², tau², sample-size weighted Z | native (METAL cannot be redistributed) |
| Heritability, rg | LDSC h2 (observed and liability scale) and genetic correlation | vendored LDSC regression code + Zenodo reference files |
| Loci | LD clumping against 1000 Genomes (EUR, EAS); distance-only fallback | PLINK 2 (auto-downloaded) |
| Mendelian randomisation | IVW, MR-Egger (+ intercept, I²_GX), weighted median, weighted/simple mode, Wald ratios, leave-one-out, Cochran's Q, Steiger, F statistics, plots, STROBE-MR | native |
| Pleiotropy | condFDR / conjFDR | [pleioFDR Python port](https://github.com/espenhgn/pleiofdr/tree/python-port) (optional, Python ≥ 3.13) |
| Fine-mapping | Wakefield ABF, credible sets | native |
| Colocalisation | coloc ABF (H0–H4) | native |
| Gene analysis | MAGMA gene, gene-set, gene-property | MAGMA, installed by the user (licence forbids bundling) |
| Lookups | Open Targets Platform (targets, diseases, credible sets, L2G, coloc), Ensembl VEP and genes in region, STRING v12, g:Profiler | public APIs, no keys |
| Reports | Markdown + HTML combining runs, warnings, decisions and references | native |

Not yet in GeNor: SuSiE-RSS, COJO runs, SMR/HEIDI runs, QTL API lookups, MR-PRESSO, MVMR,
partitioned LDSC, LAVA, MiXeR, PRS weights. See the [roadmap](docs/roadmap.md).

## How agents use it

```
skills/   knowledge: when to use each method, defaults and why, failure signatures,
          and the points where the agent must stop and ask you
genor     CLI: every command takes --json, writes a run directory, exits 0 / 2 (fixable error) /
          3 (needs your decision)
genor mcp MCP server (~25 tools): the same commands for clients without a shell;
          long analyses return a job id to poll
```

Skills: `genor` (start here), `sumstats-qc`, `meta-analysis`, `heritability-rg`,
`mendelian-randomisation`, `pleiotropy-condfdr`, `fine-mapping`, `colocalisation`,
`gene-pathway-celltype`, `open-targets-triage`, `interpretation-and-reporting`, `gwas-concepts`.
Install them for other agents with `genor skills install --agent codex|opencode|claude`.

## Documentation

- [Installation for each agent](docs/install.md)
- [CLI reference](docs/cli.md) and [MCP tools](docs/mcp.md)
- [Methods, defaults and validation](docs/methods.md), [validation details](docs/validation.md)
- [Reference data and external tools](docs/reference-data.md)
- [Where to find GWAS summary statistics](docs/data-sources.md)
- [Architecture](docs/architecture.md), [roadmap](docs/roadmap.md), [contributing](CONTRIBUTING.md)

## Platforms

Linux and macOS (x86-64 and Apple silicon). On Windows, use WSL. Python ≥ 3.11 (uv installs it
for you); condFDR/conjFDR needs Python ≥ 3.13.

## Citing

Cite GeNor (see [CITATION.cff](CITATION.cff)) **and** the methods and data you used. Every run's
`manifest.json` lists the citations for that run, and `genor report build` prints them.

## Licence

GPL-3.0-or-later (GeNor vendors GPL-3 LDSC code). Reference data and external tools keep their
own licences, listed in `genor refdata list` and `genor backend list`. Some resources are
non-commercial; GeNor is intended for academic research.
