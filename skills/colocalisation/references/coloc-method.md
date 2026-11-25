# coloc ABF: method details

What `genor coloc abf` computes (source `src/genor/methods/coloc.py`, port of `coloc::coloc.abf`).

## Data flow

1. Both traits are subset to the region (from `--region`, or the `--lead` variant in trait 1
   +/- `--window-kb`).
2. Trait 2 is harmonised to trait 1: variants are matched on rsID when most trait 1 IDs are rsIDs
   and some overlap, otherwise on chromosome:position. Trait 2 betas are flipped where alleles are
   swapped or strand-flipped; palindromic (A/T, C/G) variants are resolved with allele frequencies
   when MAF <= 0.42 and dropped otherwise (or dropped if frequencies are missing). Variants with
   incompatible alleles are dropped.
3. The two tables are merged on `variant_id`. Fewer than 10 shared variants stops the run.
4. A log ABF is computed for each variant and each trait (same formula as `fine-mapping`; see
   [../../fine-mapping/references/abf-method.md](../../fine-mapping/references/abf-method.md)).

## Hypothesis posteriors

With `l1_j`, `l2_j` the log ABFs for variant j in traits 1 and 2, and logsumexp written `LSE`:

```
log BF(H0) = 0
log BF(H1) = log p1  + LSE(l1)
log BF(H2) = log p2  + LSE(l2)
log BF(H3) = log p1 + log p2 + log( exp(LSE(l1) + LSE(l2)) - exp(LSE(l1 + l2)) )
log BF(H4) = log p12 + LSE(l1 + l2)
PP.Hk      = BF(Hk) / sum_k BF(Hk)
SNP.PP.H4_j = exp(l1_j + l2_j - LSE(l1 + l2))
```

H3 sums over all pairs of different variants (j != k). H4 sums over single variants causal for
both.

## Priors

- `p1`, `p2` (1e-4): per-variant prior of association with trait 1 only, trait 2 only.
- `p12` (1e-5): per-variant prior of association with both.

The prior odds of H4 vs H3 in a region of n variants are roughly `p12 / (p1 * p2 * (n - 1))`.
With the defaults and n = 1,000: 1e-5 / (1e-8 * 999), about 1, so the data decide. Prior odds
change with region size, so very wide windows shift the balance towards H3.

`p12 = 1e-5` implies a prior probability that a trait-1 causal variant is also causal for trait 2
of about p12 / (p1 + p12) = 9%. Wallace (2020) showed this is permissive for many trait pairs and
that PP.H4 can be high with weak evidence; she suggests `p12 = 5e-6` as a more conservative value
and recommends checking sensitivity. For GWAS-QTL pairs, where QTLs are common, some analysts
accept 1e-5; for two GWAS traits, use 5e-6 or report both.

### Prior sensitivity procedure

Run at least `--p12 1e-5`, `5e-6`, and `1e-6`. Report a locus as "robustly colocalised" only if
PP.H4 > 0.8 across the range you consider plausible. If PP.H4 crosses 0.8 only at 1e-5, say the
result is prior-dependent. (coloc's R `sensitivity()` function does this graphically; GeNor does
not reproduce that plot.)

## Thresholds and conventions

- PP.H4 > 0.8: common "strong evidence of colocalisation" cut-off (also > 0.75 and > 0.9 in use).
- PP.H3 > 0.8: evidence of distinct causal variants.
- PP.H3 + PP.H4 small: at least one trait has no signal; coloc cannot assess sharing.
- Conditional ratio PP.H4 / (PP.H3 + PP.H4) is helpful when the QTL signal is weak.

## Assumptions and how they fail

| Assumption | Failure mode | Effect |
|---|---|---|
| At most one causal variant per trait | Two signals in either trait | H3 favoured wrongly, or H4 driven by a secondary signal. Use single-signal windows or coloc-SuSiE (not in GeNor). |
| Same variants measured in both | One study imputed to a different panel; QTL file restricted to significant pairs | Missing causal variant -> lower PP.H4; `coloc_too_few`. |
| Same LD in both samples | Different ancestries | Unreliable; do not run cross-ancestry coloc. |
| Correct priors | Defaults not suitable | Run sensitivity. |
| Correct N, trait type, sdY | Wrong N, case-control loaded as quantitative | Wrong prior effect size. |

## References

- Giambartolomei C, et al. PLoS Genet 2014;10:e1004383 (coloc ABF).
- Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation
  analyses. PLoS Genet 2020;16:e1008720 (p12 = 5e-6 recommendation).
- Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal
  variants. PLoS Genet 2021;17:e1009440 (coloc-SuSiE; not in GeNor).
- Hormozdiari F, et al. Colocalization of GWAS and eQTL signals detects target genes. Am J Hum
  Genet 2016;99:1245-1260 (eCAVIAR; appears in Open Targets results as CLPP).
