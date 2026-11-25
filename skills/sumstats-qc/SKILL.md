---
name: sumstats-qc
description: Load and quality-control a GWAS summary statistics file with GeNor (genor sumstats load / liftover). Use when a user gives you a sumstats file in any format (TSV, CSV, gz, zip, GWAS-VCF, PLINK, REGENIE, SAIGE, METAL, FinnGen, Pan-UKB, GWAS Catalog), when columns are misdetected, for OR vs beta, genome build inference and liftover, ancestry, sample size, QC filters, and reading lambda_GC, lambda_1000, mean chi2, palindromic and MHC counts.
---

# Summary statistics: load and QC

## When to use

- First step for every trait, before LDSC, clumping, MR, condFDR, fine-mapping or coloc.
- When a later command fails with `not_a_trait_dir`, `build_mismatch`, `n_missing`,
  `cannot_match_ldscores` or allele-mismatch warnings: the fix is usually here.

## Inputs and prerequisites

- A local summary statistics file the user downloaded (GeNor never downloads GWAS; see
  `genor sumstats sources` / MCP `sumstats_sources`). Individual-level files are refused.
- Facts only the user can give: ancestry, sample size (or case and control counts), population
  prevalence for binary traits, and the build if inference fails.
- No reference data is needed to load. Liftover needs a chain file (0.2 to 1.2 MB).

## Steps

1. Look at the header first (`zcat file.gz | head -3`) and check the provider's README for which
   allele is the effect allele, whether N is in the file, and the build.
2. Load:

```bash
genor sumstats load bmi.tsv.gz --trait BMI --ancestry EUR --out runs/bmi --json
genor sumstats load mdd.tsv.gz --trait MDD --ancestry EUR \
    --n-cases 170756 --n-controls 329443 --population-prevalence 0.15 --out runs/mdd --json
```

   MCP: `sumstats_load(file, trait, out, build, ancestry, n, n_cases, n_controls,
   population_prevalence, effect_type, columns={"effect_allele": "A1", ...})`. To change `maf`,
   `info` or `drop_mhc` over MCP, use `run_analysis("sumstats.load", {...})`.

3. Check `summary.columns_detected` (file column -> canonical name). If anything is wrong or
   missing, rerun with overrides: `--col effect_allele=ALT --col other_allele=REF` (repeatable;
   canonical names and format recipes in [references/formats-and-columns.md](references/formats-and-columns.md)).
4. Read `summary.qc`, `warnings` and `ask_user` (see Reading the output). Resolve every `ask_user`.
5. If the build differs from the reference you will use (GeNor's LD references are GRCh37):

```bash
genor refdata fetch chain-hg38-to-hg19          # 1.2 MB; chain-hg19-to-hg38 for the other way
genor sumstats liftover runs/finngen_t2d --to GRCh37 --out runs/finngen_t2d_b37 --json
```

   MCP: `refdata_fetch("chain-hg38-to-hg19")`, `run_analysis("sumstats.liftover", {"trait": ..., "to": "GRCh37", "out": ...})`.

Outputs in the run directory: `sumstats.parquet` (standardised, QC'd table), `sumstats.json`
(trait metadata), `manhattan.png`, `qq.png` (skip with `--no-plots`), `manifest.json`.
Pass the **directory** to every later command.

## Defaults and why

| Filter | Default | Why |
|---|---|---|
| MAF | `--maf 0.01` (only if an allele frequency column exists) | Rare-variant estimates are noisy and a ~500-person 1000G panel cannot tag them for LD-based methods. |
| INFO | `--info 0.8` (only if present) | Poorly imputed variants add noise and spurious signals; 0.8 is a common GWAS cut. LDSC applies its own stricter INFO 0.9 internally. |
| P, SE, beta | drop P outside (0, 1], SE <= 0, non-finite beta | Malformed rows. |
| Alleles | drop missing or identical alleles | Cannot be harmonised. |
| Chromosomes | keep 1-22, X=23, Y=24, XY=25, MT=26 | Others are contigs or malformed. |
| Duplicates | keep the smallest P per ID and per chr:pos:alleles | One row per variant. |
| MHC (chr6:25-34 Mb) | kept and flagged; `--drop-mhc` removes it | Extreme LD; LDSC weights already exclude it, MR and condFDR have their own switches. |
| Palindromic (A/T, C/G) | kept and flagged | Resolved later by frequency in harmonisation (dropped when MAF > 0.42 or no frequency). |

Every removed row is counted in `qc.removed`; nothing is silently dropped.

Effect sizes: an odds-ratio column is converted to log-odds (`or_to_logor`). If a `beta` column is
always positive with median near 1, GeNor stops with exit 3 and asks whether it is an OR. Z is
recomputed as beta/SE when both exist. With beta and P but no SE, SE is derived (`se_from_p`).
A Z-only file loads, but it has no beta/SE, so MR cannot use it.

## Decision points

- **Ask the user** for ancestry when not given (`ancestry_unknown`). Accepted: EUR, EAS, AFR, AMR,
  SAS, MID, CSA, or `multi`. Do not infer it from the consortium name alone.
- **Ask the user** for the build when `build_inferred` is false (`build_unknown`) or the stated
  build disagrees with positions (`build_conflict`). Pass `--build GRCh37|GRCh38`.
- **Ask the user** for N, or case and control counts, when `n_missing` fires. `--n` overrides any
  N column in the file. For case-control traits give `--n-cases` and `--n-controls`; their sum
  becomes N and their ratio the sample prevalence.
- **Ask the user** for the population prevalence of a binary trait (`--population-prevalence`).
  It is the prevalence in the population, not the case fraction in the GWAS.
- **Ask the user** when the file's N is an effective N (Neff, common in PGC meta-analyses). LDSC
  liability conversion then needs sample prevalence 0.5; one way is `--n-cases` and
  `--n-controls` each set to Neff/2. Confirm before doing it.
- **Ask the user** whether the effect column is OR or beta when GeNor asks (exit 3), unless the
  README states it; then pass `--effect-type or|beta`.

## Reading the output

Key fields in `summary`: `n_variants`, `build`, `build_inferred`, `ancestry`, `n`, `n_cases`,
`n_controls`, `trait_type` (`binary` if case counts or ORs were given), `lambda_gc`, `mean_chi2`,
`n_genome_wide_sig`, and `qc` (`n_input`, `removed{...}`, `flags{non_snv, palindromic,
palindromic_ambiguous_freq, mhc}`, `n_output`, `lambda_1000`, `min_p`).

Rough reading (details and formulas in [references/qc-metrics.md](references/qc-metrics.md)):

- `lambda_gc` grows with N and polygenicity. A large GWAS with lambda_GC 1.2-1.5 is usually fine;
  use the LDSC intercept, not lambda_GC, to judge confounding (skill `heritability-rg`).
- `lambda_gc` < 0.9: the file was probably GC-corrected or is a hit list.
- `mean_chi2` < 1.02: too little signal for LDSC.
- `n_output` < 50,000: not genome-wide. LDSC and condFDR need genome-wide files.
- Many `palindromic_ambiguous_freq` or no frequency column: expect palindromic losses in MR.

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `individual_level_refused` | genotype file (.bed/.pgen/.bgen/.vcf with samples) | Refuse; ask for summary statistics. |
| `gwas_vcf_multi` | GWAS-VCF with several traits | Split with `bcftools view -s <trait>`, load each. |
| `columns_not_detected` | no effect/other allele column found | Read `error.hint` (lists file columns); add `--col`. |
| `bad_column_override` | `--col` names a column not in the header | Fix spelling (case-sensitive, exact header). |
| `no_signed_statistic` | no beta+SE, beta+P or Z | File lacks an effect column; check the right file/columns. |
| `needs_user_decision` (exit 3) | beta looks like OR; or build unknown in liftover | Ask, then rerun with `--effect-type` / `--build`. |
| `zip_multiple_members`, `empty_file`, `file_not_found` | input problems | Unzip / correct path. |
| `bad_build` | `--build`/`--to` is not GRCh37/GRCh38 (e.g. hg18) | Use GRCh37, GRCh38, hg19, hg38, b37 or b38; never guess. |
| `vcf_not_gwas_vcf` | VCF without FORMAT/sample columns (sites-only) | Need a GWAS-VCF with ES/SE/LP fields, or a tab-separated table. |
| `unsupported_compression` | `.zst` file | Decompress (`zstd -d`) or recompress with gzip. |
| `nothing_to_do` | liftover target equals current build | Skip liftover. |
| `refdata_missing` | chain file not downloaded | `genor refdata fetch chain-...` (small; still tell the user). |
| `not_a_trait_dir` (later commands) | raw file passed instead of loaded directory | Pass the `--out` directory. |

Warnings: `alleles_from_id` (serious: alleles parsed from IDs like 1:123:A:G as REF:ALT with ALT
as effect; confirm with the README), `position_from_id`, `p_z_inconsistent` (serious: P and
beta/SE disagree; SE on another scale or P GC-corrected), `build_unknown`, `build_conflict`,
`n_missing`, `ancestry_unknown`, `qc_removed_many` (serious: >30% removed; usually a mis-mapped
frequency or INFO column), `inflation_high`, `deflation`, `few_variants`, `no_eaf`, `mhc_present`,
`nonpositive_or`, `se_from_p`, `z_recomputed`, `p_from_z`, `column_note`,
`liftover_many_unmapped` (>5% did not map), `p_underflow` (P below 1e-308 read as 0; set to 1e-300).

A silent failure to watch for: the wrong allele mapped as effect allele. Nothing errors, but every
sign is flipped. Symptoms later: harmonisation `allele_mismatch`, MR or rg with an implausible sign.
Check `columns_detected` against the README every time.

## What not to conclude

- A high lambda_GC alone does not show confounding; a low one does not prove a clean GWAS.
- `build_inferred` is evidence from ~300 marker SNPs, not a guarantee; it can fail on sparse files.
- QC passing does not validate phenotype definition, ancestry labels or the effect allele.

## Related skills

`genor` (router), `meta-analysis`, `heritability-rg`, `mendelian-randomisation`,
`pleiotropy-condfdr`, `fine-mapping`, `colocalisation`, `gwas-concepts`.
