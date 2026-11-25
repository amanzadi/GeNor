# Triage of candidate targets from GWAS loci

Triage is a structured summary for humans to decide on, not a ranking GeNor computes. Present the
evidence per gene; let the user set weights.

## Suggested table

| Gene | Locus (lead, p) | Credible set size / top PP | Coding variant in set | L2G (study) | QTL coloc (tissue, PP.H4 or CLPP) | MR / direction evidence | Existing disease association score | Tractability (best modality, label) | Concerns |
|---|---|---|---|---|---|---|---|---|---|

Sources:
- Credible set: `genor finemap abf` (own data) or `genor lookup variant` (public).
- Coding variant: `genor lookup vep`.
- L2G: `genor lookup variant`.
- QTL coloc: `genor coloc abf` (own QTL data) or `genor lookup coloc`.
- MR: `genor mr run` with a gene-level exposure (e.g. cis-pQTL instruments) is possible only if the
  user supplies the QTL summary statistics; see `mendelian-randomisation`.
- Disease association and tractability: `genor lookup target`, `genor lookup disease`.

## Questions to ask per gene

1. **Is the genetic signal robust?** Genome-wide significant, replicated, not in a long-range LD
   region (MHC), consistent across ancestries if data exist.
2. **Is this the right gene at the locus?** Coding variant with high PP; QTL colocalisation in a
   relevant tissue with single-signal check; L2G >= 0.5; no stronger competitor nearby.
3. **What is the direction of effect?** Does the risk allele increase or decrease the gene's
   expression or function? Target modulation must go in the protective direction. GeNor's coloc
   does not report direction and does not save the harmonised betas; look up the top shared
   variant in both trait files, align them to the same effect allele by hand, and compare signs.
4. **Is it tractable?** Tractability modalities and labels.
5. **What is already known?** Existing association score and diseases; known drugs on the target.
6. **Safety and pleiotropy**: associations with many other diseases (`associatedDiseases`) can
   signal on-target side effects. GeNor does not fetch Open Targets safety data; say so.

## Reading genetics as drug-target support

Targets with human genetic support are more likely to succeed in clinical development (King et al.
2019; Minikel et al. 2024). This is a population-level enrichment, not a guarantee for any single
target. Support is strongest when the causal gene is confidently assigned (coding variants,
Mendelian genes), and weaker for L2G-only assignments.

## Pitfalls

- Ranking by L2G alone repeats the gold-standard bias.
- Ranking by overall association score rewards genes that are already well studied.
- Ignoring direction of effect.
- Merging genes from different credible sets of one locus as if they were independent hits.
- Presenting triage as a recommendation to develop a drug. GeNor is research software.

## References

- King EA, Davis JW, Degner JF. Are drug targets with genetic support twice as likely to be
  approved? PLoS Genet 2019;15:e1008489.
- Minikel EV, et al. Refining the impact of genetic evidence on clinical success. Nature
  2024;629:624-629.
- Buniello A, et al. Open Targets Platform. Nucleic Acids Res 2025;53:D1467-D1475.
