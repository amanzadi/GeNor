---
name: meta-analysis
description: Variant-level meta-analysis of two or more GWAS of the same trait with GeNor (genor meta run): fixed-effect inverse-variance (IVW), DerSimonian-Laird random effects, sample-size weighted Z, Cochran's Q, I2 and tau2, cross-ancestry caveats. Use when the user wants to combine cohorts, consortia or biobank GWAS, or asks about heterogeneity between studies.
---

# Meta-analysis

## When to use

- Combining GWAS of the **same trait** from independent cohorts to gain power.
- Checking whether a hit is consistent across cohorts (heterogeneity).

Not for combining different traits (use `heritability-rg` or `pleiotropy-condfdr`) and not a
replacement for MR-MEGA style trans-ancestry modelling (not in GeNor yet).

## Inputs and prerequisites

- Two or more traits already loaded with `genor sumstats load` (skill `sumstats-qc`), all on the
  **same genome build** (lift over first).
- The same effect allele convention after loading (GeNor aligns alleles, but a mis-mapped effect
  column in one study flips its signs silently; check `columns_detected` for every study).
- N for each study if you use `--scheme samplesize`.
- No reference data or backends needed.

## Steps

```bash
genor meta run runs/ukb_t2d runs/finngen_t2d_b37 runs/biobankjapan_t2d \
    --scheme ivw --name T2D_meta --out runs/t2d_meta --json
```

MCP: `run_analysis("meta", {"traits": [...], "scheme": "ivw", "name": "T2D_meta", "out": ...})`.

The output directory is itself a loaded trait (`sumstats.parquet` + `sumstats.json`), so it can go
straight into `genor clump run`, `genor ldsc h2`, and so on. The full table with heterogeneity
columns is `meta_results.tsv.gz`.

## Defaults and why

- `--scheme ivw` (default): fixed-effect inverse-variance weighting, beta = sum(w b) / sum(w),
  w = 1/SE^2, as METAL `SCHEME STDERR` (Willer et al. 2010). Most efficient when all studies
  estimate the same effect on the same scale. GeNor also reports DerSimonian-Laird random effects
  (`beta_random`, `se_random`, `p_random`), Cochran's Q (`het_Q`, `het_df`, `het_p`), I2 (`het_I2`)
  and tau2 per variant (DerSimonian and Laird 1986; Higgins and Thompson 2002).
- `--scheme samplesize`: weighted Z with weights sqrt(N), as METAL `SCHEME SAMPLESIZE`. Use when
  studies report effects on different scales (different phenotype transformations, linear vs
  logistic models) or lack usable SEs. It needs N for every study (`n_missing` otherwise).
- Alleles are aligned to the first study that reports each variant: swapped alleles flip sign,
  strand complements are accepted, anything else is dropped and counted (`meta_allele_mismatch`).
  Palindromic SNPs are matched as written (as METAL does), so all studies must be on the same strand.
- Variants present in only some studies are kept; `n_studies` says how many contributed.
  `n` is the sum of N over contributing studies.

## Decision points

- **Ask the user** whether the studies are independent (no shared participants). Overlapping
  samples make fixed-effect SEs too small; GeNor does not correct for overlap.
- **Ask the user** whether the phenotype definitions and effect scales match. If not (for example
  inverse-normal transformed vs raw units, or different case definitions), prefer
  `--scheme samplesize` and say that the result has no interpretable beta.
- **Ask the user** before meta-analysing across ancestries. GeNor allows it but warns
  (`cross_ancestry_meta`) and labels the result's ancestry `multi`, which makes LD-based steps ask
  again which panel to use.
- **Ask the user** before reporting loci driven by one study (`n_studies` = 1 at a hit).

## Reading the output

`summary`: `scheme`, `n_studies`, `n_variants`, `n_in_all_studies`, `n_genome_wide_sig`,
`alignment` (per study: rows, `flipped`, `allele_mismatch_dropped`), and for IVW
`n_sig_with_high_I2` (genome-wide significant variants with I2 > 50%).

Heterogeneity per variant (IVW only):

- Q: weighted squared deviations from the fixed-effect estimate; chi2 with k-1 df under
  homogeneity. Low power with few studies, so a non-significant Q does not show homogeneity.
- I2 = max(0, (Q - df)/Q): share of variation beyond chance. Rough guide: 25% low, 50% moderate,
  75% high (Higgins et al. 2003). Unstable with 2-3 studies.
- tau2: between-study variance of true effects (DerSimonian-Laird moment estimator); random-effects
  weights are 1/(SE^2 + tau2).

When to prefer which estimate:

- Discovery: the fixed-effect P is the standard GWAS test (it tests "no effect in any study" most
  powerfully). Report it, and report I2 and the random-effects estimate alongside for hits.
- A hit with high I2 and a much weaker random-effects P is driven by a subset of studies; flag it.
- DerSimonian-Laird with few studies underestimates uncertainty; do not over-read `p_random` either.

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `meta_too_few` | fewer than two traits given | Pass at least two loaded trait directories. |
| `build_mismatch` | studies on different builds | `genor sumstats liftover` the odd ones out. |
| `n_missing` | sample-size scheme with a study lacking N | Reload that study with `--n`, or use `--scheme ivw`. |
| `meta_needs_position` | studies mix rsID and chr:pos keys and one lacks chromosome/position | Reload that study with chromosome/position columns, or use rsIDs everywhere. |
| `meta_no_overlap` | no variant in more than one study | Check builds and ID schemes match; liftover or reload. |
| `bad_param` | unknown `--scheme` | `ivw` or `samplesize`. |
| `not_a_trait_dir` | raw file passed | Pass the `sumstats load` output directory. |

Warnings: `cross_ancestry_meta` (studies span ancestries), `meta_allele_mismatch` (rows dropped
for incompatible alleles; many means a build or ID mismatch), `heterogeneous_hits` (significant
variants with I2 > 50%).

## Known limits

- The sample-size scheme output has Z and P but no beta or SE, so it cannot feed MR.
- The meta result is saved as a quantitative trait with total N; case and control counts are not
  carried over. For a case-control meta-analysis, liability-scale h2 and binary-trait handling
  downstream will not be automatic. Summed N is not effective N when case fractions differ.
- Genomic control per study, overlap correction, and trans-ancestry meta-regression (MR-MEGA) are
  not in GeNor yet.

## What not to conclude

- Low I2 does not prove effects are identical; with few studies the test has little power.
- A cross-ancestry fixed-effect hit does not imply the same causal variant in each ancestry.
- A meta-analysis of overlapping cohorts overstates significance.

## Related skills

`sumstats-qc`, `heritability-rg`, `fine-mapping`, `gwas-concepts`, `interpretation-and-reporting`.
