---
name: heritability-rg
description: SNP heritability (h2) and genetic correlation (rg) with LD score regression in GeNor (genor ldsc h2 / ldsc rg). Use when the user asks for heritability, liability-scale h2 of a disease, genetic correlation between traits, the LDSC intercept or ratio (confounding vs polygenicity), or evidence of sample overlap before Mendelian randomisation.
---

# Heritability and genetic correlation (LDSC)

## When to use

- How much trait variance common SNPs explain (h2), on the observed or liability scale.
- Whether inflation is confounding or polygenicity (intercept, ratio).
- Genetic correlation between two traits, and the cross-trait intercept as evidence of sample
  overlap for MR.

## Inputs and prerequisites

- Genome-wide traits loaded with `genor sumstats load` (not hit lists), with N. For binary traits,
  case and control counts and the population prevalence.
- Same ancestry as the LD scores. GeNor's registry has **EUR only** (1000G Phase 3, GRCh37,
  HapMap3 SNPs). Matching is by rsID when more than half the IDs are rsIDs; otherwise positions
  must be GRCh37.
- Reference data (about 25 MB in total; tell the user before fetching):

```bash
genor refdata fetch ldsc-ldscores-eur ldsc-weights-eur      # MCP: refdata_fetch(name) per dataset
```

## Steps

```bash
genor ldsc h2 runs/bmi --out runs/h2_bmi --json
genor ldsc h2 runs/mdd --pop-prev 0.15 --out runs/h2_mdd --json
genor ldsc rg runs/bmi runs/mdd --pop-prev2 0.15 --out runs/rg_bmi_mdd --json
```

MCP: `run_analysis("ldsc.h2", {"trait": ..., "pop_prev": 0.15, "out": ...})`,
`run_analysis("ldsc.rg", {"trait1": ..., "trait2": ..., "pop_prev1": ..., "pop_prev2": ..., "out": ...})`.
`ldsc rg` estimates both h2 values too, so a separate `h2` run is optional when you need both.

Options: `--intercept 1` constrains the h2 intercept (only when the user has reason to; see
below). `--ref` / `--w` take other refdata names, but only the EUR sets exist today.

## Defaults and why

Following the reference LDSC implementation (Bulik-Sullivan et al. 2015, Nat Genet 47:291;
rg: Bulik-Sullivan et al. 2015, Nat Genet 47:1236):

- Munging as `munge_sumstats.py`: SNVs only, strand-ambiguous SNPs removed, INFO >= 0.9 and
  MAF > 0.01 when present, SNPs with N < 90th percentile / 1.5 removed (`ldsc_low_n_removed`).
- HapMap3 SNPs via the LD scores; MHC excluded through the regression weights.
- Free intercept, two-step estimator (chi2 cutoff 30), 200 block jackknife for SEs.
- h2 refers to SNPs with MAF > 5% (`M_5_50`), the LDSC convention.
- Liability scale (Lee et al. 2011, AJHG 88:294): h2_liab = h2_obs x K^2 (1-K)^2 / (P (1-P) z^2),
  with K the population prevalence, P the sample prevalence (cases / (cases + controls)) and z
  the standard normal density at the liability threshold. GeNor needs both P (from
  `--n-cases/--n-controls` at load time) and K (`--pop-prev`, or `--population-prevalence` at load).

## Decision points

- **Ask the user** for the population prevalence K of each binary trait. Do not use the GWAS case
  fraction. If unsure, report the observed-scale h2 and liability h2 under the user's stated range of K.
- **Ask the user** before constraining the intercept (`--intercept 1`). It narrows SEs only if
  there is truly no confounding and no sample overlap; it is wrong for GC-corrected files.
- **Ask the user** what to do for non-EUR GWAS: GeNor has no non-EUR LD scores yet, and a
  different-ancestry panel gives biased estimates. Stop rather than substitute.
- For rg between different ancestries GeNor warns (`rg_cross_ancestry`); a cross-ancestry method
  such as Popcorn is needed and is not in GeNor yet.

## Reading the output

`ldsc h2` summary: `h2_observed`, `h2_observed_se`, `h2_z`, `intercept`, `intercept_se`,
`lambda_gc`, `mean_chi2`, `n_snps`, `ratio`, `ratio_se` (when mean chi2 > 1), and for binary traits
`h2_liability`, `h2_liability_se`, `sample_prevalence`, `population_prevalence`.

`ldsc rg` summary: `rg`, `rg_se`, `z`, `p`, `gencov`, `gencov_se`, `gencov_intercept`,
`gencov_intercept_se`, `gencov_intercept_p`, `n_snps`, `h2_trait1{...}`, `h2_trait2{...}`.

Interpretation (detail in [references/ldsc-interpretation.md](references/ldsc-interpretation.md)):

- **Intercept** near 1: inflation is polygenic signal. Above 1: confounding (stratification,
  cryptic relatedness), or attenuation in very large GWAS. Below 1: probably GC-corrected.
- **Ratio** = (intercept - 1) / (mean chi2 - 1): share of inflation not explained by polygenicity.
  Near 0 is good; 0.1-0.2 is common in large GWAS; high ratios warrant caution.
- **h2 z < 4**: h2 is too imprecise for a meaningful rg with this trait (atlas recommendation).
- **mean chi2 < 1.02**: LDSC is unreliable; do not report h2 or rg as findings.
- **rg** is the same on observed and liability scales; prevalences only affect the h2 entries.
- **rg outside [-1, 1]** (GeNor flags |rg| > 1.2): almost always a low or noisy h2. Do not interpret.
- **Cross-trait intercept** (`gencov_intercept`) is approximately N_s x rho_p / sqrt(N1 N2), where
  N_s is the number of overlapping participants and rho_p their phenotypic correlation. Clearly
  non-zero (GeNor flags |z| > 3 as `sample_overlap_signal`; `gencov_intercept_p` gives the test):
  evidence of sample overlap (or shared confounding). Pass this to the MR decision on `--overlap`.
  Near zero does not prove no overlap (rho_p may be small).

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `refdata_missing` | LD scores / weights not downloaded | Tell the user (about 13 MB each), then `genor refdata fetch ...`. |
| `needs_user_decision` (exit 3) | ancestry unknown or `multi` | Ask; reload with `--ancestry`, or accept EUR for a mostly-EUR multi-ancestry GWAS. |
| `ancestry_mismatch` | GWAS ancestry is not EUR | No matching LD scores in GeNor; stop and say so. |
| `n_missing` | trait has no N | Reload with `--n` or `--n-cases/--n-controls`. |
| `cannot_match_ldscores` | no rsIDs and positions not GRCh37 | Lift over to GRCh37, or load a file with rsIDs. |
| `ldsc_too_few_shared` | < 1,000 SNPs shared between traits in rg | Check that both are genome-wide and use rsIDs. |
| `refdata_corrupt`, `refdata_no_m` | incomplete LD score files | `genor refdata remove <name>` then fetch again. |

Warnings: `ldsc_low_chi2` (serious), `ldsc_h2_weak` (h2 z < 4), `ldsc_intercept_high` (> 1.1,
serious), `ldsc_intercept_low` (< 0.95), `ldsc_few_snps` (< 200,000 SNPs overlap the LD scores;
LDSC expects about 1.2M HapMap3 SNPs), `ldsc_low_n_removed`, `rsid_from_position`,
`rg_cross_ancestry` (serious), `rg_negative_h2` (h2 <= 0 for one trait; rg undefined),
`rg_out_of_bounds` (serious), `sample_overlap_signal`.

The `ask_user` "What is the population prevalence of <trait>?" also appears when K was given
but case and control counts were not loaded; liability conversion needs both.

## What not to conclude

- h2 is SNP heritability of common HapMap3 variants, not twin or total heritability.
- A non-zero rg is not causation and says nothing about direction (that is MR's question).
- rg is genome-wide; local shared signals can exist with rg near 0 (see `pleiotropy-condfdr`).
- A significant rg with tiny h2 or out-of-bounds values is not a finding.
- Partitioned or stratified LDSC (enrichment by annotation, cell type) is not in GeNor yet.

## Related skills

`sumstats-qc`, `mendelian-randomisation` (uses the cross-trait intercept), `pleiotropy-condfdr`,
`gwas-concepts`, `interpretation-and-reporting`.
