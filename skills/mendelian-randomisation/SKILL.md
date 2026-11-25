---
name: mendelian-randomisation
description: Two-sample Mendelian randomisation from GWAS summary statistics with GeNor (genor mr run): instrument selection and clumping, F statistics, allele harmonisation, IVW (multiplicative random effects), MR-Egger with intercept and I2_GX, weighted median, weighted and simple mode, Steiger directionality, heterogeneity, leave-one-out, binary exposures and outcomes, sample overlap, and the STROBE-MR checklist. Use when the user asks whether an exposure causally affects an outcome, or for MR, causal inference, instruments or pleiotropy-robust estimates.
---

# Mendelian randomisation (two-sample)

## When to use

The user asks whether exposure X has a causal effect on outcome Y, and has GWAS summary statistics
for both. MR estimates the effect of lifelong, genetically predicted differences in X on Y, valid
only under three instrumental-variable (IV) assumptions:

1. **Relevance**: the variants are associated with X (checked: p-value threshold, F statistics).
2. **Independence**: no confounding of the variant-outcome association (not testable; supported
   by same-ancestry, PC-adjusted GWAS).
3. **Exclusion restriction**: variants affect Y only through X, no horizontal pleiotropy (probed
   with heterogeneity, Egger intercept, robust estimators; never proven).

## Inputs and prerequisites

- Exposure and outcome loaded with `genor sumstats load` (skill `sumstats-qc`), ideally the same
  ancestry and with beta, SE, allele frequencies and N.
- For LD clumping (default): PLINK 2 and the 1000G panel matching the exposure ancestry.

```bash
genor backend install plink2                  # MCP: backend_install("plink2")
genor refdata fetch g1000-eur --yes           # 275 MB: ask the user first. g1000-eas for EAS.
```

  The panels are GRCh37 and matched by rsID: the exposure must be on GRCh37 and have rsIDs.
- Run `genor ldsc rg` first if you can: its cross-trait intercept is evidence on sample overlap
  (skill `heritability-rg`).

## Steps

1. Agree with the user which trait is the exposure and which the outcome, the units, and the
   sample overlap (see Decision points).
2. Run:

```bash
genor mr run --exposure runs/bmi --outcome runs/mdd --overlap none --out runs/mr_bmi_mdd --json
```

   MCP: `run_analysis("mr", {"exposure": ..., "outcome": ..., "overlap": "none", "out": ...})`.
   MR runs with bootstraps can exceed the MCP wait; poll `job_status` / `job_result`.
3. Read warnings and `ask_user` before the estimates. Inspect `scatter.png`, `leave_one_out.png`,
   `funnel.png`, `forest.png`.
4. Optionally check reverse direction by swapping `--exposure` and `--outcome` (a separate run).
5. Hand `strobe_mr.md` to the user; it marks the items only they can write.
6. `genor report build runs/mr_bmi_mdd --out runs/mr_bmi_mdd/report.md`.

Options: `--p-threshold 5e-8`, `--clump-mode plink|distance`, `--panel g1000-eur`,
`--clump-r2 0.001`, `--clump-kb 10000`, `--instruments file` (first column = rsIDs; skips the
threshold and clumping), `--exclude-mhc`, `--overlap none|partial|substantial|unknown`,
`--nboot 1000`, `--seed 1`.

## Defaults and why

| Setting | Default | Why |
|---|---|---|
| Instrument threshold | p < 5e-8 | Genome-wide significance limits weak and false-positive instruments. |
| Clumping | r2 < 0.001 within 10 Mb, 1000G panel | Near-independent instruments (TwoSampleMR defaults; Hemani et al. 2018). |
| Palindromes | resolved by allele frequency; dropped if MAF > 0.42 or frequencies missing | Strand cannot be inferred near 0.5 (TwoSampleMR `action = 2`). |
| Main estimator | IVW, multiplicative random effects | Efficient with valid instruments; random-effects SEs absorb heterogeneity (Burgess et al. 2013). |
| Sensitivity | MR-Egger, weighted median, weighted mode, simple mode, IVW fixed effects | Each is consistent under different pleiotropy assumptions (Bowden 2015, 2016; Hartwig 2017). |
| Bootstrap | 1000 replicates, seed 1 | SEs for median and mode estimators; seed for reproducibility. |
| Missing instruments in outcome | dropped (no LD proxies) | Proxies are not in GeNor yet. |

Estimator details, formulas and assumptions: [references/estimators-and-diagnostics.md](references/estimators-and-diagnostics.md).

## Decision points

- **Ask the user** for the exposure and outcome definitions: which GWAS, what units (SD, log-odds,
  raw units), the direction of the question, and whether a bidirectional analysis is wanted.
- **Ask the user** about sample overlap (`--overlap`). If unstated, GeNor still runs and returns
  an `ask_user`. Show them the LDSC cross-trait intercept as evidence. Overlap with weak
  instruments biases estimates towards the confounded observational association.
- **Ask the user** before relaxing `--p-threshold` (for example to 5e-6) when there are few or no
  instruments (`mr_no_instruments`, `mr_few_instruments`). Relaxing adds weak instruments: bias
  towards the null without overlap, towards the observational estimate with overlap.
- **Ask the user** before using `--clump-mode distance` when the LD panel is unavailable or the
  ancestry has no panel. Distance pruning ignores LD, so instruments may be correlated and SEs too
  small (`distance_clumping`, serious). Label results as approximate.
- **Ask the user** before downloading the 1000G panel (275-308 MB).
- **Ask the user** whether to exclude the MHC (`--exclude-mhc`) when instruments fall in it
  (`mhc_present` at load).
- **Ask the user** before using user-supplied instruments (`--instruments`): GeNor does not
  re-check their threshold or independence.

## Reading the output

`summary` keys: `n_candidates` (exposure SNPs passing the threshold), `clump` (mode, r2, kb, panel,
`n_index`), `n_instruments_clumped`, `harmonisation` (`aligned`, `swapped`, `strand_flipped`,
`palindromic_dropped`, `allele_mismatch_dropped`, `n_kept`), `n_missing_in_outcome`, `nsnp`,
`estimates[]` (`method`, `nsnp`, `b`, `se`, `pval`), `main_or` and `main_or_95ci` (binary outcome),
`heterogeneity_ivw` (Cochran's Q), `heterogeneity_egger` (Rucker's Q'), `pleiotropy_egger`
(intercept, se, pval), `I2_GX`, `F` (`F_mean`, `F_min`, `n_weak_F_lt_10`), `steiger`
(`r2_exposure`, `r2_outcome`, `correct_causal_direction`, `steiger_pval`),
`direction_consistent_across_methods`, `sample_overlap`, `strobe_mr` (counts of filled / partial / human items).

Files: `harmonised.tsv`, `mr_results.tsv`, `single_snp.tsv` (Wald ratios), `leave_one_out.tsv`,
`scatter.png`, `forest.png`, `funnel.png`, `leave_one_out.png`, `strobe_mr.md`, `strobe_mr.json`.

Scale:

- Estimates are per unit of the exposure's GWAS scale (often per SD) on the outcome's scale.
- **Binary outcome**: `b` is a log odds ratio; `main_or = exp(b)` is the OR per unit of exposure.
- **Binary exposure** (`binary_exposure`): `b` is per unit increase in log-odds of liability to
  the exposure, not the effect of having the condition. Multiply by ln 2 = 0.693 for the effect per
  doubling of the odds (Burgess and Labrecque 2018).

How to read a result set:

- Main answer: the IVW row. Then ask whether Egger, weighted median and mode agree in direction
  and roughly in size (`direction_consistent_across_methods`). Agreement across estimators with
  different assumptions is the strongest evidence MR can give; disagreement is a finding to report.
- Heterogeneity (Q p < 0.05): instruments disagree more than chance. Expected with many SNPs; it
  is why IVW uses random effects. Look at outliers in the scatter and leave-one-out plots.
- Egger intercept p < 0.05: directional pleiotropy; IVW may be biased; give weight to Egger,
  median and mode. Egger is only interpretable when I2_GX is high (> 0.9).
- Steiger: `correct_causal_direction` false means instruments explain more variance in the outcome
  than the exposure: reverse causation or invalid instruments.

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `mr_no_instruments` | no exposure SNP below the threshold | **Ask the user** about relaxing; or the exposure GWAS is underpowered. |
| `refdata_missing` | 1000G panel not downloaded | Ask; `genor refdata fetch g1000-eur --yes`. |
| `confirm_download` | fetch without `--yes` for a > 100 MB dataset | Ask; rerun with `--yes` (MCP `confirmed_by_user=true`). |
| `backend_missing` | PLINK 2 not installed | `genor backend install plink2`. |
| `build_mismatch` | exposure not GRCh37 for the panel | `genor sumstats liftover <exposure> --to GRCh37`. |
| `needs_user_decision` (exit 3) | exposure ancestry or build unknown, or multi-ancestry | Ask; reload with `--ancestry` / `--build`. |
| `ancestry_mismatch` | no panel for this ancestry | Ask about `--clump-mode distance` or stop. |
| `clump_needs_rsid` | exposure has no rsIDs | Use `--clump-mode distance` (ask) or a file with rsIDs. |
| `cannot_match_variants` | no shared rsIDs or positions with the outcome | Check builds and IDs of both traits. |
| `backend_failed` | PLINK 2 exited non-zero | Read `error.hint` and the log path in `error.details`. |

Warnings: `mr_ancestry_mismatch` (serious), `distance_clumping` (serious),
`instruments_missing_in_outcome`, `mr_few_instruments` (< 3, serious: no Egger, median or mode),
`weak_instruments` (F < 10), `mr_heterogeneity`, `mr_directional_pleiotropy` (serious),
`egger_low_i2gx` (I2_GX < 0.9), `steiger_reverse` (serious), `binary_exposure`,
`palindromic_dropped_no_freq`, `allele_mismatch` (serious: > 5% of shared variants with
incompatible alleles; build or ID mismatch), `chromosome_not_in_panel` (e.g. chrX candidates skipped).

## What not to conclude

- MR does not prove causation. Say "evidence consistent with a causal effect, under the IV
  assumptions", and name the diagnostic that would most undermine it.
- A non-significant Egger intercept does not show absence of pleiotropy; it has low power.
- The estimate is a lifelong effect of genetically predicted exposure, not the effect of an
  intervention, a drug, or a change in adulthood. Do not translate it into clinical advice.
- A null MR result with weak instruments or a small outcome GWAS is not evidence of no effect.
- Results do not transfer across ancestries without testing.
- Not in GeNor yet: MR-PRESSO, multivariable MR, LD proxies, Steiger filtering of SNPs (only the
  directionality test is run), MR-RAPS, CAUSE, radial MR, SMR/HEIDI with QTLs. Do not claim
  these were done.

## STROBE-MR

`strobe_mr.md` covers all 20 items: GeNor fills items it knows from the run (methods, numbers,
estimates, diagnostics), partially fills others, and marks the rest **HUMAN**. Never write HUMAN
items (rationale, objectives, interpretation, funding, conflicts) yourself without the user's input.
Item list and GeNor's status for each: [references/strobe-mr.md](references/strobe-mr.md).

## Related skills

`sumstats-qc`, `heritability-rg` (overlap evidence), `colocalisation` (locus-level check of a
shared causal variant for a key instrument), `open-targets-triage`, `interpretation-and-reporting`,
`gwas-concepts`.
