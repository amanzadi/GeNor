---
name: gwas-concepts
description: Teaching notes for explaining GWAS and post-GWAS concepts to the user while running GeNor - variants and alleles, LD and ancestry, association statistics, 5e-8, power, winner's curse, lambda GC and the LDSC intercept, design pitfalls (sample overlap, collider bias, ascertainment, relatedness), heritability and PRS interpretation, and an overview of advanced methods (fine-mapping, coloc, MR, TWAS, condFDR, MiXeR, LAVA) with what GeNor does and does not implement. Use when the user asks "what is...", "why does...", "explain...", or seems unfamiliar with a term in a GeNor result.
---
# GWAS concepts for explaining as you go

Use this when a user needs a concept explained, or when a GeNor result can only be understood
with some background. Explain briefly, in plain language, tied to the user's actual result. Point
to the reference file for depth; do not paste it whole.

## How to explain

- Start from the number the user is looking at (e.g. "your lambda GC is 1.18"), then the concept,
  then what it means for their analysis.
- One concept at a time. Offer more depth instead of giving it unasked.
- Say what is uncertain or debated. Do not overstate.
- Separate what GeNor computed from general background.

## Index

| Area | Covers | Reference |
|---|---|---|
| Genetic variation | SNPs, indels, alleles (effect/other, REF/ALT), strand and palindromic SNPs, allele frequency and MAF, rsIDs, genome builds, imputation INFO | [references/genetic-variation.md](references/genetic-variation.md) |
| Population genetics | Linkage disequilibrium (r2, D'), LD blocks, the MHC, ancestry and LD differences, population stratification, Hardy-Weinberg equilibrium, reference panels | [references/population-genetics.md](references/population-genetics.md) |
| Statistical principles | Association tests (linear, logistic, mixed models), beta/SE/Z/P, odds ratios, multiple testing and 5e-8, power, winner's curse, lambda GC, mean chi-square, LDSC intercept and ratio, Bayes factors and posteriors | [references/statistical-principles.md](references/statistical-principles.md) |
| Design pitfalls | Sample overlap, collider and selection bias, ascertainment, cryptic relatedness, case-control vs population prevalence, proxy phenotypes, meta-analysis heterogeneity | [references/design-pitfalls.md](references/design-pitfalls.md) |
| Heritability and PRS | Twin vs SNP heritability, observed vs liability scale, genetic correlation, PRS construction and portability, what PRS can and cannot say (GeNor does not score genotypes) | [references/heritability-and-prs.md](references/heritability-and-prs.md) |
| Advanced methods | Fine-mapping, colocalisation, MR, TWAS/SMR, condFDR/conjFDR, MiXeR, LAVA, partitioned heritability, MAGMA; which are in GeNor | [references/advanced-methods.md](references/advanced-methods.md) |

Read a reference with your skills loader, or with `genor skills show gwas-concepts/references/<file>.md`
(MCP: `get_guidance(topic="gwas-concepts/references/<file>.md")`).

## Quick answers to common questions

- **"Is my GWAS inflated?"** lambda GC and mean chi-square rise with both confounding and true
  polygenic signal in large samples. The LDSC intercept separates them: near 1 means little
  confounding. See statistical-principles.
- **"Why 5e-8?"** Bonferroni for about one million independent common-variant tests in European
  genomes. See statistical-principles.
- **"Is the top SNP the causal one?"** Usually not known: LD makes many variants look alike. That
  is what fine-mapping addresses. See population-genetics and advanced-methods.
- **"Is the nearest gene the causal gene?"** Often not. See `gene-pathway-celltype`.
- **"What does h2 = 0.2 mean?"** About 20% of trait variance in this population is explained
  additively by the common variants measured, on the stated scale. Not "20% genetic for a person".
  See heritability-and-prs.
- **"Can I use this to predict my risk?"** GeNor does not compute individual scores, and population
  summary results do not give individual predictions. Research use only.
- **"Does genetic correlation mean one trait causes the other?"** No. See design-pitfalls and
  advanced-methods (MR).
- **"Why does ancestry matter?"** LD patterns and allele frequencies differ, so effect estimates,
  LD references, and PRS do not transfer cleanly. See population-genetics.

## Where GeNor shows each concept

| Concept | Where it appears in GeNor output |
|---|---|
| lambda GC, mean chi2, lambda_1000 | `genor sumstats load` QC report; warnings `inflation_high`, `deflation` |
| LDSC intercept, ratio | `genor ldsc h2`; warnings `ldsc_intercept_high`, `ldsc_intercept_low` |
| Sample overlap | `genor ldsc rg` cross-trait intercept (`sample_overlap_signal`); `genor mr run --overlap` |
| Build | `genor sumstats load` (`build_inferred`, `build_conflict`, `build_unknown`) |
| Palindromic SNPs, strand | Harmonisation reports in MR, coloc and meta |
| Liability scale | `genor ldsc h2 --pop-prev`, `genor sumstats load --population-prevalence` |
| Winner's curse, weak instruments | `genor mr run` F statistics (`weak_instruments`) |
| Single causal variant | `single_causal_assumption` in fine-mapping and coloc |

## Boundaries

- GeNor works only with summary statistics. It does not run GWAS, score individuals (PRS), or
  handle individual-level genotypes (it refuses them: `individual_level_refused`).
- GeNor is research software, not a medical device. Do not interpret results for an individual's
  health.

## Related skills

`genor` (router), `sumstats-qc`, `heritability-rg`, `meta-analysis`, `mendelian-randomisation`,
`pleiotropy-condfdr`, `fine-mapping`, `colocalisation`, `gene-pathway-celltype`,
`open-targets-triage`, `interpretation-and-reporting`.
