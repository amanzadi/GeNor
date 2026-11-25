# Wakefield ABF fine-mapping: method details

This is what `genor finemap abf` computes (source: `src/genor/methods/coloc.py`, a port of
`coloc::finemap.abf` from coloc 5.x).

## Model

For each variant j in the region, with estimated effect `beta_j` and variance `V_j = SE_j^2`:

- Under H0 (no association) `beta_j ~ N(0, V_j)`.
- Under H1 (variant j causal) `beta_j ~ N(0, V_j + W)`, where `W` is the prior variance of the true
  effect.

Wakefield's approximate Bayes factor (log scale) is

```
r_j    = W / (W + V_j)
lABF_j = 0.5 * ( log(1 - r_j) + r_j * z_j^2 ),   z_j = beta_j / SE_j
```

With exactly one causal variant (or none), the posterior for variant j is

```
PP_j    = p1 * ABF_j / ( (1 - n*p1) + sum_k p1 * ABF_k )
PP_null = (1 - n*p1) / ( (1 - n*p1) + sum_k p1 * ABF_k )
```

where `n` is the number of variants in the region and `p1` is the per-variant prior (default
1e-4). The sum of the variant PPs plus `PP_null` is 1.

Only marginal statistics are used. Under a single causal variant the marginal z-scores carry all
the information needed, so no LD matrix is required. With two or more causal variants that
reasoning fails.

## The prior effect size W = sd_prior^2

- Quantitative trait: `sd_prior = 0.15 * sdY`. The 0.15 is on the standardised scale (trait SD = 1).
- Case-control trait: `sd_prior = 0.2` on the log-odds-ratio scale, so a priori 95% of causal
  effects have OR between about exp(-0.39) = 0.68 and exp(0.39) = 1.48.

These are the coloc defaults (Wakefield 2009; Giambartolomei et al. 2014). Larger W favours
variants with large z; smaller W shrinks differences between variants.

## sdY: the standard deviation of a quantitative trait

`beta` must be on the trait's measurement scale for the prior to make sense. If betas are per SD
(standardised or inverse-normal transformed), sdY = 1. If they are in raw units (kg/m2, mmol/L),
sdY is the trait SD in those units.

GeNor estimates sdY like `coloc:::sdY.est`: for each variant, `Var(beta_j) ~ sdY^2 / (2 N f_j (1 - f_j))`,
so GeNor regresses `2 N f (1 - f)` on `1 / V` through the origin and takes the square root of the
slope. It needs `effect_allele_frequency` and `n` for every variant. The estimate is reported as
`sdY_estimated`. If it cannot be estimated, GeNor uses sdY = 1 and warns `sdy_assumed`.

Checks:
- `sdY_estimated` close to 1 for a trait the user says is standardised: consistent.
- Far from 1 for a standardised trait: N or allele frequencies may be wrong (e.g. N is the total
  meta-analysis N while many variants have a smaller per-variant N).

## P-value route (no beta/SE)

If beta/SE are missing, GeNor reconstructs `V_j` from MAF and N:

- Quantitative: `V_j = 1 / (2 N f (1 - f))` (assumes sdY = 1), `sd_prior = 0.15`.
- Case-control: `V_j = 1 / (2 N f (1 - f) s (1 - s))`, with `s` the case fraction, `sd_prior = 0.2`.

and uses `|z_j| = Phi^-1(1 - p_j / 2)`. This route is less accurate: N is often approximate and
the sign is lost (sign does not matter for ABF). This route assumes sdY = 1 and raises
`sdy_assumed`.

## Credible sets

GeNor sorts variants by PP and keeps the smallest top set whose cumulative PP reaches `level`.
The set is built from PPs renormalised without the null (PP / (1 - PP_null)), i.e. conditional
on one causal variant being in the region; `finemap.tsv` keeps the unconditional PPs.

- If `PP_null` is near 0 (typical for genome-wide significant loci), the two coincide.
- If `PP_null` is high, the set is still computed but the region probably has no signal; say so.
- GeNor refuses `n_variants x p1 >= 1` (`bad_param`): the null prior would be negative.

Interpretation: under the model, the causal variant is in the set with probability >= level. This
calibration holds only if the model is right (one causal variant, typed, correct prior).

## Effect of region size

- `n * p1` is the prior probability that the region contains a causal variant. At the defaults a
  1 Mb window with 2,000 variants has prior 0.2; 10,000 variants would give prior 1 and a null prior
  of 0 (GeNor does not guard against this; the log of a non-positive number gives undefined PPs).
- Wider windows add more variants in weak LD with the lead. They rarely change the credible set
  for a single strong signal, but they increase the chance of including a second independent signal.
- Too narrow a window can exclude the causal variant if LD extends far (common near the MHC or
  in regions of low recombination). The credible set is always conditional on the causal variant
  being inside the region.

## Coverage of variants

The credible set can only contain variants in the summary statistics after QC (`sumstats load`
filters MAF < 0.01 and INFO < 0.8 by default). Rare causal variants and variants poorly imputed
in the GWAS are invisible. Example data in GeNor is HapMap3 only (about 98k variants), so its
credible sets are illustrative.

## References

- Wakefield J. Bayes factors for genome-wide association studies: comparison with P-values.
  Genet Epidemiol 2009;33:79-86.
- Giambartolomei C, et al. Bayesian test for colocalisation between pairs of genetic association
  studies using summary statistics. PLoS Genet 2014;10:e1004383.
- Maller JB, et al. Bayesian refinement of association signals for 14 loci in 3 common diseases.
  Nat Genet 2012;44:1294-1301 (origin of credible sets in this form).
