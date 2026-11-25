# MR estimators and diagnostics in GeNor

Notation per instrument j: exposure effect bx_j (SE sx_j), outcome effect by_j (SE sy_j),
Wald ratio by_j / bx_j. All outcome effects are aligned to the exposure effect allele.
GeNor's estimators are validated against TwoSampleMR 0.7.x (`tests/golden/`).

## Estimators

| Method | What GeNor computes | Valid when | Notes |
|---|---|---|---|
| Wald ratio | by/bx, SE = sy/abs(bx) | one instrument | Used when only one SNP survives. |
| IVW (main) | weighted regression of by on bx through the origin, weights 1/sy^2; SE multiplied by the residual SE when it exceeds 1 (multiplicative random effects) | all instruments valid, or pleiotropy balanced with mean zero | Burgess et al. 2013. Normal-based p. |
| IVW (fixed effects) | same point estimate, SE not inflated | no heterogeneity | Report only alongside the main IVW. |
| MR-Egger | regression of by on bx with intercept (SNPs oriented so bx > 0), same weights | InSIDE: pleiotropic effects independent of instrument strength | Bowden et al. 2015. t-distribution, k-2 df. Imprecise; needs >= 3 SNPs and high I2_GX. |
| Weighted median | median of Wald ratios with inverse-variance weights (second-order); bootstrap SE | at least 50% of the weight from valid instruments | Bowden et al. 2016. |
| Weighted mode | mode of the smoothed Wald-ratio distribution, weights from second-order SEs; bootstrap SE (MAD), bandwidth phi = 1 | the largest group of instruments with the same ratio is valid (ZEMPA) | Hartwig et al. 2017. t-distribution, k-1 df. |
| Simple mode | as above, unweighted | ZEMPA | Less efficient. |

Sensitivity estimators need at least 3 instruments (`mr_few_instruments` otherwise).

## Diagnostics

- **F statistic** per SNP: (bx/sx)^2. F < 10 is the conventional weak-instrument flag
  (`weak_instruments`). Mean F also bounds weak-instrument bias: in one-sample or fully
  overlapping two-sample MR the relative bias towards the observational estimate is about 1/F.
- **Cochran's Q** (IVW) with k-1 df and **Rucker's Q'** (Egger) with k-2 df. Q - Q' large means
  Egger fits better, consistent with directional pleiotropy.
- **Egger intercept**: the average pleiotropic effect when bx = 0. p < 0.05 raises
  `mr_directional_pleiotropy`.
- **I2_GX** (Bowden et al. 2016, Int J Epidemiol 45:1961): heterogeneity of bx relative to its
  uncertainty. Below 0.9 Egger is diluted towards the null (NOME violation); flagged
  `egger_low_i2gx`. Simulation extrapolation correction is not in GeNor yet.
- **Steiger directionality** (Hemani et al. 2017): compares variance explained in exposure and
  outcome, with r^2 per SNP from F/(N-2+F) using b, SE and N for both traits. For binary traits
  this is an observed-scale approximation, not the liability-scale r^2 that TwoSampleMR derives
  from log ORs and prevalence; read it cautiously for binary outcomes. Requires N in both traits.
- **Leave-one-out**: IVW refitted without each SNP (`leave_one_out.tsv`, `.png`). A single SNP
  that moves the estimate across the null is a finding; look it up (skill `open-targets-triage`).
- **Single-SNP Wald ratios**: `single_snp.tsv`, forest plot; outliers suggest pleiotropy.
- **Funnel plot**: asymmetry suggests directional pleiotropy.

## Sample overlap

Two-sample MR assumes independent samples. With overlap, weak-instrument bias points towards the
confounded observational estimate, roughly in proportion to the overlap fraction and 1/F
(Burgess, Davies and Thompson 2016, Genet Epidemiol 40:597). With strong instruments (high mean F)
the bias is small even with substantial overlap. Without overlap, weak-instrument bias is towards
the null. GeNor records the user's statement in `sample_overlap` and STROBE-MR item 4 needs the
human to describe it. Evidence: the LDSC cross-trait intercept (`gencov_intercept`).

## What changes the instrument set

- Instruments absent from the outcome are dropped (no LD proxies); count in `n_missing_in_outcome`.
- Palindromic SNPs with MAF > 0.42, or without frequencies, are dropped.
- Allele mismatches are dropped (a large count points to a build or ID mismatch).
- `--exclude-mhc` drops instruments in chr6:25-34 Mb.
- With PLINK clumping, candidates absent from the 1000G panel are lost (summary key
  `n_candidates_not_in_panel_or_clumped` counts them together with SNPs clumped away).

## Binary traits

- Outcome binary: b = log OR per unit exposure. OR = exp(b); 95% CI exp(b +/- 1.96 SE).
- Exposure binary: b per unit log-odds of liability. For effect per doubling of exposure odds,
  multiply b by ln 2. The estimate reflects liability to the exposure, not the condition itself,
  and the exclusion restriction is harder to defend (Burgess and Labrecque 2018, Eur J Epidemiol 33:947).
