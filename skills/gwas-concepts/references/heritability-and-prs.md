# Heritability, genetic correlation, and polygenic scores

## Heritability

**Heritability** is the proportion of trait variance in a population explained by genetic
variation. It describes a population at a time and in an environment, not an individual, and
not how "fixed" a trait is.

- **Twin / family heritability** (broad or narrow sense): estimated from relatives' resemblance;
  includes all genetic effects (rare variants, sometimes non-additive effects). Often higher.
- **SNP heritability (h2_SNP)**: variance explained additively by the common variants measured.
  Estimated from summary statistics by LDSC (`genor ldsc h2`) or from genotypes (GCTA-GREML, not
  in GeNor). Usually well below twin estimates (the "missing heritability" gap: rare variants,
  imperfect tagging, non-additive effects, and possibly inflated twin estimates).
- LDSC h2 uses HapMap3 SNPs and LD scores from a reference panel; it is a lower-bound-type estimate
  for common-variant heritability under the LDSC model.

### Observed vs liability scale

For binary traits, h2 on the observed 0-1 scale depends on the case fraction in the sample. The
liability-threshold model assumes an underlying continuous liability; h2 on the liability scale
is comparable across studies. Conversion needs the sample prevalence (from case/control counts)
and the **population prevalence** K, which the user must provide (`genor ldsc h2 --pop-prev`).
Liability h2 is sensitive to K; report K and consider a range.

## Genetic correlation (rg)

The correlation between the genetic effects of two traits across the genome (`genor ldsc rg`).
- Ranges -1 to 1. Estimates outside this range (`rg_out_of_bounds`) usually reflect low h2 in one
  trait.
- Reflects shared genetic architecture: pleiotropy, causal relationships, shared confounding, or
  correlated measurement. It does not say which.
- Genome-wide rg can be near 0 while local regions correlate positively and negatively (LAVA,
  not in GeNor).

## Polygenic risk scores (PRS)

A PRS sums an individual's allele counts weighted by GWAS effect sizes. GeNor does **not** compute
PRS or score genotypes; it works only with summary statistics and refuses individual-level data.

What to tell users who ask:
- Construction methods include clumping and thresholding (C+T), and Bayesian shrinkage methods
  (LDpred2, PRS-CS, SBayesR). They need an LD reference and a target genotype dataset.
- Predictive accuracy is bounded by SNP heritability and GWAS sample size; for most traits a PRS
  explains a minority of variance.
- **Portability**: PRS built from European-ancestry GWAS predict substantially worse in other
  ancestries (often several-fold lower accuracy in African ancestries) because of LD and
  frequency differences, and differing environments.
- A PRS ranks people in a population; it is not a diagnosis. Absolute risk depends on prevalence,
  age, and non-genetic factors.
- Clinical use requires validated, regulated tools. GeNor is research software.

## Explaining to a user

"An SNP heritability of 0.08 for MDD on the liability scale means common variants explain about
8% of variation in liability to depression in populations like this one, assuming a population
prevalence of 15%. It does not mean any one person's depression is 8% genetic."
