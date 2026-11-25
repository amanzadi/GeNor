# Population genetics for GWAS

## Linkage disequilibrium (LD)

Nearby variants are inherited together because recombination between them is rare. LD is the
non-random association of alleles at different loci.

- **r2**: squared correlation between allele dosages at two variants (0 to 1). If r2 = 1, the two
  variants carry the same information. GWAS association at one variant "tags" variants in high r2.
- **D'**: normalised measure of recombination history; D' = 1 with low r2 happens when allele
  frequencies differ. For GWAS tagging, r2 is the relevant measure.
- **LD blocks**: regions of high LD separated by recombination hotspots; length varies from a few
  kb to hundreds of kb, longer in European than African ancestries.
- **Consequence**: a GWAS signal points to a region, not a variant. Many variants in high LD with
  the causal one are also significant. Fine-mapping tries to separate them; clumping
  (`genor clump run`) collapses them into one index variant per signal.

### Long-range LD regions

The **MHC** (major histocompatibility complex, chr6 about 25-34 Mb on GRCh37) has extreme LD and
diversity. Many methods exclude it: LDSC weights exclude it, pleioFDR excludes it by default
(`--keep-mhc` overrides), and MR may exclude it (`--exclude-mhc`). GeNor flags it (`mhc_present`)
and can drop it on load (`--drop-mhc`). Other long-range LD regions (e.g. inversions on 8p23 and
17q21) cause similar trouble.

## Reference panels

Methods that need LD (clumping, LDSC, MAGMA, COJO, SuSiE-RSS) use a reference panel when the GWAS
genotypes are not available. GeNor uses 1000 Genomes Phase 3 panels (`g1000-eur`, `g1000-eas`,
GRCh37) and precomputed LD scores. The panel must match the GWAS ancestry; LD mismatch can produce
false secondary signals in conditional analysis and biased gene p-values. 1000 Genomes EUR has
about 500 individuals, so LD estimates for low-frequency variants are noisy.

## Ancestry

- Allele frequencies and LD differ between populations because of different demographic
  histories. African-ancestry genomes have more diversity and shorter LD.
- GWAS effect estimates for a tagging variant differ across ancestries when LD with the causal
  variant differs, even if the causal effect is identical.
- Most large GWAS are European ancestry. Results, LD references, and PRS do not transfer reliably.
- GeNor asks the user for ancestry on load and checks panel ancestry (`ancestry_mismatch`), and
  warns on cross-ancestry meta-analysis, rg, and MR.
- "Ancestry" here means genetic similarity to reference groups, not race or ethnicity.

## Population stratification

If allele frequencies and the trait both differ between subgroups in the sample (e.g. geography),
variants can appear associated without any causal role. GWAS control this with principal
components of genotypes and linear mixed models. Residual stratification inflates test statistics;
the LDSC intercept is the summary-statistics diagnostic (see statistical-principles).

## Hardy-Weinberg equilibrium (HWE)

In a randomly mating population without selection, genotype frequencies are p2, 2pq, q2 for allele
frequencies p and q. Strong deviation in controls often indicates genotyping error. HWE filtering
happens during GWAS QC on individual-level data; summary statistics usually do not include HWE
p-values, and GeNor does not test HWE.

## Explaining to a user

"Variants near each other travel together through generations, so when one is associated with
your trait, its neighbours are too. The lead SNP is the strongest marker of the region, not
necessarily the cause. That is why we clump to count independent signals and fine-map to rank
candidates."
