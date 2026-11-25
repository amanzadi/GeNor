# From a locus to candidate genes

No single method reliably names the causal gene at a GWAS locus. Collect several lines of evidence
and present them side by side.

## Evidence types and GeNor commands

| Evidence | GeNor | Strength | Main bias |
|---|---|---|---|
| Coding consequence (missense, stop, splice) of a credible-set variant | `genor lookup vep` | Strong when the variant has high PP | Few GWAS signals are coding |
| Gene overlapping or near the credible set | `genor lookup genes --region` | Weak to moderate | Gene density; long genes overlap more |
| Open Targets L2G score | `genor lookup variant` | Moderate; integrates distance, QTL coloc, chromatin, coding | Trained on known gold-standard genes; favours well-studied genes |
| QTL colocalisation in a relevant tissue | `genor coloc abf`; `genor lookup coloc` | Moderate to strong if PP.H4 high and single signals | QTLs are common; a variant can affect several genes |
| MAGMA gene association | `genor magma gene` | Gene-level aggregation | Neighbouring genes share signal through LD |
| Rare-variant burden, Mendelian disease genes, knockout phenotypes | Not in GeNor; Open Targets association evidence covers some (`genor lookup target`) | Strong when present | Coverage is limited |

## Defining the region for positional mapping

Choose and state one rule, for example:

- Lead variant +/- 250 kb (matches `--preset loci` clumping window), or
- Bounds of the credible set from `fine-mapping`, extended to the LD r2 > 0.6 boundary if known,
  or
- Fixed +/- 500 kb for a liberal candidate list.

The region must be in the trait's build. `genor lookup genes` defaults to `--build GRCh37`;
`genor lookup vep` defaults to `--build GRCh38`. Set both explicitly. (rsIDs are build-independent,
so VEP consequences do not depend on the build, but `location` does.)

## A per-locus evidence table

| Locus (lead) | Credible set size | Coding variant in set (gene) | Genes within region | L2G top gene (score) | QTL coloc PP.H4 (gene, tissue) | MAGMA gene P | Notes |
|---|---|---|---|---|---|---|---|

Rank candidate genes by the number and strength of independent lines of evidence. Keep "no data"
distinct from "evidence against".

## Common traps

- Picking the nearest gene and then citing the nearest gene in enrichment analyses as independent
  support.
- Treating MAGMA's top gene in a locus as the causal gene: MAGMA tests association aggregated over
  a window, so all genes overlapping a strong signal get small p-values.
- Ignoring non-coding genes (lncRNAs) or, conversely, listing every pseudogene. Filter by
  `biotype` and state the filter.
- Mixing builds: positions from a GRCh37 GWAS queried with `--build GRCh38` give the wrong genes.
- In the MHC (chr6:25-34 Mb on GRCh37) and other long-range LD regions, positional mapping is
  close to meaningless; flag such loci.
