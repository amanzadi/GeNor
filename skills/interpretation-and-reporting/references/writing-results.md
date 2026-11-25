# Writing results

Plain statements, numbers with uncertainty and units, and the assumptions that the number depends
on.

## Effect sizes and scales

| Result | Report | Scale notes |
|---|---|---|
| GWAS variant effect | beta (SE) or OR (95% CI), p, effect allele, EAF | Quantitative: per effect allele in trait units (or SD units if standardised). Binary: log-OR; OR = exp(beta), 95% CI = exp(beta +/- 1.96 SE). |
| Meta-analysis | Combined beta (SE), p; heterogeneity (Cochran's Q, I2); model (fixed / random) | See `meta-analysis`. |
| SNP heritability (LDSC) | h2 (SE) on observed scale; for binary traits liability-scale h2 (SE) with the population prevalence used; intercept (SE); ratio | Liability-scale h2 depends on the assumed population prevalence; state it. See `heritability-rg`. |
| Genetic correlation | rg (SE), p; cross-trait intercept | rg outside [-1, 1] is not interpretable. |
| MR | IVW estimate (95% CI), p; number of instruments; sensitivity estimators; heterogeneity; Egger intercept; F statistics; Steiger | Estimate is per unit of the exposure's GWAS scale (per SD, per log-odds for binary exposure) on the outcome's scale (log-OR -> OR for binary outcomes). See `mendelian-randomisation`. |
| Fine-mapping | Region, number of variants, PP.null, credible set size at 95%, top variant and PP, priors | Conditional on one causal variant per region. |
| Colocalisation | Region, nsnps, PP.H0-H4, priors (p1, p2, p12), p12 sensitivity | PP.H4 > 0.8 is a convention. |
| condFDR / conjFDR | Number of loci at the threshold used, threshold | See `pleiotropy-condfdr`. |
| MAGMA | Gene P with Bonferroni threshold; gene-set / property beta, SE, P with correction method | Competitive test. |
| Enrichment | Term, source, adjusted p (g:SCS), intersection / term size, background | State background and version. |

Always give 95% confidence intervals (or SEs) next to point estimates. For posterior
probabilities give the priors.

## Multiple testing

- Genome-wide significance p < 5e-8 accounts for about one million independent common-variant tests
  in European-ancestry genomes. It is not a law; it assumes that testing scope.
- Several outcomes, exposures, or traits: correct across them (Bonferroni or FDR) and state the
  number of tests, e.g. MR of one exposure on 10 outcomes: 0.05 / 10 = 0.005.
- MAGMA genes: 0.05 / number of genes (GeNor's `bonferroni_p`). Gene sets: correct over sets tested.
- Posterior probabilities (fine-mapping, coloc) are not corrected with Bonferroni; report how many
  loci, genes, and tissues were tested, so readers can judge selection.
- Exploratory analyses chosen after looking at results (e.g. picking the tissue with the best
  PP.H4) must be labelled as exploratory.

## Causal language

| Evidence | Acceptable wording | Avoid |
|---|---|---|
| GWAS association | "associated with", "the locus is associated" | "gene X causes" |
| Fine-mapping | "the credible set contains", "most likely causal variant under a single-signal model" | "the causal variant is" |
| Colocalisation PP.H4 > 0.8 | "consistent with a shared causal variant" | "eQTL proves the mechanism", "gene X mediates" |
| MR with robust sensitivity analyses | "genetic evidence consistent with a causal effect of X on Y", "under the IV assumptions" | "X causes Y", "proves" |
| MR with heterogeneity / pleiotropy warnings | "the estimate may be biased by pleiotropy" | presenting IVW alone |
| L2G / nearest gene | "prioritised gene", "candidate gene" | "the causal gene" |
| Enrichment | "genes near associated loci are enriched for ..." | "the disease is caused by pathway ..." |
| rg | "shared genetic influences" | "the same causal variants", "one causes the other" |

## Reporting warnings and human decisions

- List every `serious` warning and every `warning` that affects a reported result in the
  limitations, with what it implies.
- State every decision the human made (build, ancestry, N, population prevalence, sample overlap,
  instrument threshold, priors, tissue, loci selected) and who made it.
- Report unanswered `ask_user` questions as open limitations.
- If a result was computed with GeNor defaults the user did not review, say so.

## Methods paragraph template

"Summary statistics for <trait> (<source>, N = <n>, <ancestry>, <build>) were processed with
GeNor <version> (QC: MAF >= <maf>, INFO >= <info>; <MHC handling>). <Analyses with parameters>.
Reference data: <refdata with versions>. External tools: <backends with versions>. All inputs are
identified by SHA-256 checksums in the run manifests."

Fill from `manifest.json` (`params`, `inputs`, `refdata`, `backends`, `genor_version`), not from
memory.

## Research-use notice

Keep the notice from the report header in any derived document: GeNor is research software, not a
medical device, and must not be used for clinical decisions, diagnosis, or commercial services.
Results about populations say nothing reliable about an individual.

## References

- Skrivankova VW, et al. STROBE-MR statement. JAMA 2021;326:1614-1621.
- Little J, et al. STrengthening the REporting of Genetic Association Studies (STREGA). PLoS Med
  2009;6:e22.
