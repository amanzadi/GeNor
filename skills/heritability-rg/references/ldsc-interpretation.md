# LDSC interpretation notes

## The model

For SNP j with LD score l_j, under a polygenic model:

E[chi2_j] = N h2 l_j / M + N a + 1

The slope estimates h2 (M = number of SNPs with MAF > 5% in the reference, `M_5_50`). The
intercept (1 + N a) captures confounding from stratification and cryptic relatedness, which does
not scale with LD score. Bulik-Sullivan et al. 2015, Nat Genet 47:291.

For two traits:

E[z1_j z2_j] = sqrt(N1 N2) rho_g l_j / M + N_s rho_p / sqrt(N1 N2)

rg = rho_g / sqrt(h2_1 h2_2). The cross-trait intercept carries sample overlap (N_s) times the
phenotypic correlation in the overlap (rho_p). Bulik-Sullivan et al. 2015, Nat Genet 47:1236.

## Intercept and ratio

| Pattern | Likely explanation |
|---|---|
| intercept ~1.00-1.05, high lambda_GC | Polygenic signal in a large sample. Fine. |
| intercept > 1.1 (`ldsc_intercept_high`) | Residual stratification, relatedness, or misspecified N. In GWAS of several hundred thousand people the intercept also rises from attenuation bias (Loh et al. 2018), so read the ratio too. |
| intercept < 0.95 (`ldsc_intercept_low`) | Genomic control applied to the file, or N overstated. h2 is then biased downward. |
| ratio ~0 | Inflation is polygenic. |
| ratio 0.1-0.2 | Common in large biobank GWAS; mild. |
| ratio > 0.2 | Report as a caution; some inflation may be confounding. |

Ratio is reported only when mean chi2 > 1 and is unstable when mean chi2 is close to 1.

## Thresholds GeNor uses and where they come from

- mean chi2 < 1.02 (`ldsc_low_chi2`): the LDSC documentation warns estimates are unreliable at
  this signal level; commonly cited minimum.
- h2 z < 4 (`ldsc_h2_weak`): Bulik-Sullivan et al. 2015 (rg atlas) restricted genetic correlation
  analyses to traits with h2 z > 4.
- |rg| > 1.2 (`rg_out_of_bounds`): rg is a correlation; estimates beyond 1 come from sampling noise
  in h2 or gencov, usually with weak h2.
- |gencov intercept / SE| > 3 (`sample_overlap_signal`).
- fewer than 200,000 SNPs overlapping the LD scores (`ldsc_few_snps`): SEs grow; the standard
  HapMap3 set has about 1.2M SNPs.

## Liability scale

Observed-scale h2 for a case-control GWAS depends on the case fraction. The liability-scale
transform (Lee et al. 2011, AJHG 88:294) needs:

- K: the **population** prevalence (lifetime or point prevalence matching the case definition).
- P: the **sample** prevalence, N_cases / (N_cases + N_controls). For a meta-analysis analysed with
  effective N, P = 0.5.

Different but defensible K values can change h2_liab noticeably for rare diseases. Report the K you
used and why. rg does not depend on the scale.

## What LDSC cannot do in GeNor today

- Partitioned / stratified LDSC (functional or cell-type enrichment): not in GeNor yet.
- Non-EUR LD scores: not in the registry yet.
- Cross-ancestry genetic correlation (Popcorn): not in GeNor yet.
- Constraining the cross-trait intercept from the CLI: not exposed.
