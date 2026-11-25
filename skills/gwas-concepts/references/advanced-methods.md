# Advanced post-GWAS methods: overview and GeNor status

| Method | Question | In GeNor? | Command / skill |
|---|---|---|---|
| LD clumping | How many independent signals? | Yes (PLINK 2 or distance) | `genor clump run` |
| Meta-analysis | Combine studies | Yes (IVW fixed/random, sample-size) | `genor meta run`, `meta-analysis` |
| LDSC h2 / rg | SNP heritability; genome-wide genetic correlation | Yes | `genor ldsc h2|rg`, `heritability-rg` |
| Fine-mapping, Wakefield ABF | Which variant in a region, assuming one signal | Yes | `genor finemap abf`, `fine-mapping` |
| SuSiE / SuSiE-RSS, FINEMAP | Fine-mapping with multiple signals | No (planned later) | - |
| GCTA-COJO | Conditional / joint analysis of multiple signals | No (GCTA installable, no command; planned later) | - |
| Colocalisation, coloc ABF | Do two traits share a causal variant in a region? | Yes | `genor coloc abf`, `colocalisation` |
| coloc-SuSiE, eCAVIAR | Coloc with multiple signals | No (Open Targets results via `genor lookup coloc`) | `open-targets-triage` |
| Two-sample MR | Causal effect of exposure on outcome | Yes (IVW, Egger, weighted median/mode, Steiger, STROBE-MR) | `genor mr run`, `mendelian-randomisation` |
| MR-PRESSO, CAUSE, multivariable MR | Robust / multivariable MR | No | - |
| SMR / HEIDI | MR using QTLs as instruments, with a test for linkage | No (SMR binary installable, no command) | - |
| TWAS (FUSION, S-PrediXcan) | Gene expression imputed association | No | - |
| condFDR / conjFDR (pleioFDR) | Boost discovery using a related trait; shared loci | Yes (Python port) | `genor pleio run`, `pleiotropy-condfdr` |
| MiXeR | Polygenicity, discoverability, polygenic overlap | No (container only, not wired) | - |
| LAVA | Local genetic correlation | No (R backend listed, not wired) | - |
| S-LDSC / LDSC-SEG | Partitioned heritability by annotation / tissue | No | - |
| MAGMA | Gene, gene-set, gene-property analysis | Yes, with user-installed MAGMA | `genor magma gene|sets`, `gene-pathway-celltype` |
| PRS | Individual genetic scores | No (GeNor does not score genotypes) | - |

## Short explanations

- **Fine-mapping**: ranks variants in a region by posterior probability of being causal, using
  the strength of association and (for multi-signal methods) LD. Output: credible sets.
- **Colocalisation**: asks whether the association signals of two traits in one region come from
  the same variant. Common use: GWAS + eQTL to link a locus to a gene.
- **Mendelian randomisation**: uses genetic variants as instrumental variables to estimate the
  causal effect of an exposure on an outcome. Valid only if variants affect the outcome only
  through the exposure (no horizontal pleiotropy), are associated with the exposure, and are not
  confounded. Sensitivity methods relax some of these assumptions.
- **TWAS**: predicts gene expression from genotypes using a reference panel, then tests predicted
  expression against the trait. Significant genes are often not causal: LD and co-regulation
  produce many hits per locus.
- **SMR/HEIDI**: MR with a top cis-QTL as instrument; HEIDI tests whether the GWAS and QTL signals
  reflect one variant or two in LD.
- **condFDR**: re-estimates the false discovery rate for trait 1 variants given their association
  with trait 2, boosting discovery when the traits share genetics. **conjFDR**: variants associated
  with both traits (shared loci, either direction of effect).
- **MiXeR**: models the number of causal variants and their effect size distribution; bivariate
  MiXeR estimates how many causal variants two traits share, even when rg is near zero.
- **LAVA**: genetic correlation within genomic regions; reveals local sharing hidden by opposing
  directions genome-wide.
- **Partitioned heritability (S-LDSC)**: tests whether annotations (e.g. tissue-specific
  enhancers or cell-type-specific genes) carry more heritability than expected.

When a user asks for a method GeNor lacks, say so plainly, explain what it would answer, and
suggest whether an existing GeNor method addresses part of the question.
