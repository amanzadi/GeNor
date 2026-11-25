# Statistical principles

## Association tests

- **Quantitative traits**: linear regression of trait on allele dosage (0, 1, 2) with covariates
  (age, sex, principal components). Beta = change in trait per effect allele.
- **Binary traits**: logistic regression; beta = log odds ratio per effect allele, OR = exp(beta).
- **Mixed models** (BOLT-LMM, REGENIE, SAIGE, fastGWA): account for relatedness and structure and
  increase power. Betas from some mixed models are on a transformed scale (e.g. BOLT-LMM on binary
  traits gives linear-scale betas that need conversion to log-OR).
- **Summary statistics**: per variant beta, SE, Z = beta / SE, p-value, N, allele frequency.
  GeNor derives missing pieces where possible (e.g. SE from beta and P: warning `se_from_p`).

## Odds ratios and scales

- OR = exp(beta); 95% CI = exp(beta +/- 1.96 x SE).
- OR is not relative risk; they are similar only when the outcome is rare.
- Per-allele effects are small for common variants (OR typically 1.02-1.2).

## Multiple testing and 5e-8

About one million independent common variants exist in European genomes (fewer independent tests
than variants, because of LD). Bonferroni 0.05 / 1e6 = 5e-8. This is the genome-wide significance
convention. African-ancestry and rare-variant analyses have more independent tests; stricter
thresholds are sometimes used. Sub-threshold variants are not "null"; many are true but
underpowered.

## Power

Power to detect a variant depends on the variance it explains: for a quantitative trait,
2 f (1 - f) beta^2 (in trait-SD units), times N. The non-centrality of the chi-square test is
about N x variance explained. A variant explaining 0.01% of variance needs N of roughly
400,000 for 80% power at 5e-8. Case-control power depends on the effective sample size
4 / (1/N_cases + 1/N_controls).

## Winner's curse

Variants selected because they passed a threshold tend to have overestimated effects in the
discovery sample, especially near the threshold. Consequences: replication effects are smaller;
MR instruments selected and estimated in the same sample bias towards the observational estimate
with sample overlap; PRS weights are inflated. Use independent samples for selection and
estimation where possible.

## Genomic inflation: lambda GC and mean chi-square

- **lambda GC** = median observed chi-square / 0.4549 (the median of a 1-df chi-square). Under no
  association anywhere it is 1.
- In a polygenic trait with large N, lambda rises well above 1 from true signal alone. It
  increases with sample size, so GeNor also reports `lambda_1000`, computed as
  1 + (lambda - 1) x 1000 / median N (lambda rescaled to a sample of 1,000); for case-control
  studies GeNor uses 1 + (lambda - 1) x 500 x (1/N_cases + 1/N_controls).
- **mean chi-square**: average of Z^2. LDSC needs mean chi-square above about 1.02 to be reliable
  (warning `ldsc_low_chi2`).
- lambda < 1 (warning `deflation`) suggests the file was already GC-corrected or filtered.

## LDSC intercept and ratio

LD score regression regresses each variant's chi-square on its LD score (how much variation it
tags). True polygenic signal increases chi-square in proportion to LD score (slope). Confounding
(stratification, cryptic relatedness) inflates chi-square regardless of LD score (intercept).

- **Intercept** about 1: inflation is mostly polygenicity. Above about 1.1 (warning
  `ldsc_intercept_high`): suggests confounding or misspecified N. Below 1 (`ldsc_intercept_low`):
  often prior GC correction.
- **Ratio** = (intercept - 1) / (mean chi-square - 1): share of inflation not explained by
  polygenicity. Large in absolute terms only when mean chi-square is clearly above 1.
- Note: with very large N, the intercept can be somewhat above 1 because of attenuation and
  model mis-specification, not only confounding.

## Bayes factors and posterior probabilities

Fine-mapping and coloc use Bayes factors: how much more likely the data are under "variant
associated" than "no association", given a prior effect size distribution. Combined with prior
probabilities (p1, p12) they give posteriors (PP). PPs depend on the prior, are not p-values, and
are calibrated only if the model (e.g. one causal variant) holds.

## Explaining to a user

"lambda GC of 1.18 with N = 500,000 is expected for a polygenic trait. The LDSC intercept of 1.02
says almost none of that inflation looks like confounding."
