# QC metrics in `summary.qc`

All metrics are computed after filtering, on the rows that are saved.

| Field | Definition | How to read it |
|---|---|---|
| `n_input`, `n_output` | rows before / after QC | A drop over 30% raises `qc_removed_many`. |
| `removed.bad_p` | P missing or outside (0, 1] | A few are normal; many suggest P is -log10 P (map it as `neg_log10_p`). |
| `removed.bad_se` | SE <= 0 or missing | |
| `removed.bad_beta` | non-finite beta | |
| `removed.bad_alleles` | allele missing or both alleles equal | |
| `removed.bad_eaf` | frequency outside [0, 1] | Frequency column may be a count (AC) or percent. |
| `removed.maf_lt_<x>` | minor allele frequency below `--maf` | Huge counts: frequency column mis-mapped. |
| `removed.info_lt_<x>` | INFO below `--info` | Huge counts: `R2`/`Rsq` may not be INFO in this file. |
| `removed.bad_chromosome` | chromosome not 1-26 after X/Y/XY/MT recoding | |
| `removed.duplicates` | duplicate IDs or chr:pos:alleles (smallest P kept) | |
| `removed.mhc` | only with `--drop-mhc` | |
| `flags.non_snv` | indels and multi-base alleles (kept) | LDSC and clumping use SNVs. |
| `flags.palindromic` | A/T or C/G SNPs (kept) | |
| `flags.palindromic_ambiguous_freq` | palindromic with frequency 0.4-0.6 | Likely dropped in MR harmonisation. |
| `flags.mhc` | variants in chr6:25-34 Mb | |
| `lambda_gc` | median chi2 / 0.4549 (median of chi2 with 1 df) | See below. |
| `lambda_1000` | 1 + (lambda_GC - 1) x 1000 / N (quantitative); 1 + (lambda_GC - 1) x 500 x (1/N_cases + 1/N_controls) (case-control) | Inflation rescaled to 1,000 individuals (1,000 cases and 1,000 controls for binary traits). |
| `mean_chi2` | mean of Z^2 | LDSC needs > 1.02. |
| `n_genome_wide_sig` | P < 5e-8 | |
| `min_p`, `n_median` | smallest P, median N | |

## Inflation: lambda_GC is not a confounding test

Under polygenicity, lambda_GC increases with sample size even with no confounding
(Yang et al. 2011, Eur J Hum Genet 19:807). The LDSC intercept separates confounding from
polygenic signal (Bulik-Sullivan et al. 2015): intercept near 1 with high lambda_GC means
polygenic signal. Run `genor ldsc h2` before calling a GWAS inflated.

Signals that a file is unsuitable for genome-wide methods:

- `lambda_gc` < 0.9 (`deflation`): genomic control already applied, or only significant hits.
- `n_output` < 50,000 (`few_variants`): a hit list or a subset.
- `mean_chi2` < 1.02: too little polygenic signal; LDSC and condFDR are unreliable.

## Genome build inference

GeNor compares positions of a bundled set of marker SNPs (about 300) with GRCh37 and GRCh38. It
matches by rsID when at least 50 markers are present, otherwise by chr:pos. The evidence
(`match_grch37`, `match_grch38`, `markers`) is in the `build_inferred` / `build_unknown` warning
details. A clear call needs > 50% matches on one build and 3x more than the other (or 5x for
sparse files).

## Liftover

`genor sumstats liftover` uses UCSC chain files natively. Variants outside aligned blocks, mapped to
another chromosome, or mapped to the reverse strand are dropped and counted (`n_unmapped`). More
than 5% unmapped raises `liftover_many_unmapped`. rsIDs and alleles are kept unchanged.
