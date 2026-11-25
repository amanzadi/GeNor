---
name: gene-pathway-celltype
description: Move from GWAS loci to genes, pathways, and tissues or cell types - positional mapping (genor lookup genes), variant consequences (genor lookup vep), Open Targets L2G, MAGMA gene, gene-set and gene-property analysis (genor magma gene|sets), g:Profiler enrichment (genor lookup enrich) and STRING networks (genor lookup string). Use when the user asks which genes a locus implicates, what pathways or biology are enriched, which tissues or cell types are involved, or mentions MAGMA, GO, Reactome, KEGG, gene-set or PPI analysis.
---
# Loci to genes, pathways, tissues and cell types

GWAS finds variants. Getting to genes and biology needs several lines of evidence, each with its
own biases. GeNor offers:

| Step | Command | Needs |
|---|---|---|
| Genes overlapping a region | `genor lookup genes --region chr:start-end --build GRCh37` | Internet (Ensembl REST) |
| Variant consequences | `genor lookup vep rs1 rs2 ... --build GRCh38` | Internet (Ensembl VEP) |
| Locus-to-gene (L2G) from public fine-mapping | `genor lookup variant CHR_POS_REF_ALT` | Internet (Open Targets), GRCh38 IDs |
| Gene-based association | `genor magma gene` | User-installed MAGMA + MAGMA auxiliary files |
| Gene-set / gene-property (tissue, cell type) | `genor magma sets` | MAGMA + gene-set or gene-property file |
| Functional enrichment of a gene list | `genor lookup enrich` | Internet (g:Profiler) |
| Protein interaction network | `genor lookup string` | Internet (STRING v12) |

Not in GeNor yet: S-LDSC / LDSC-SEG partitioned heritability, TWAS (FUSION, S-PrediXcan), SMR,
PoPS, FUMA-style eQTL / chromatin-interaction mapping, DEPICT. Say so when asked.

## When to use

- After clumping (and ideally `fine-mapping`) to name candidate genes per locus.
- Genome-wide: MAGMA gene and gene-set analysis on the full summary statistics.
- To ask which tissues or cell types are enriched for trait heritability-related genes (MAGMA
  gene-property with an expression specificity file).

## Inputs and prerequisites

- A loaded trait (`genor sumstats load`), and for loci, `genor clump run` output.
- For MAGMA (see [references/magma.md](references/magma.md)):
  - MAGMA itself, downloaded by the user (licence forbids redistribution) and registered:
    `genor backend add magma /path/to/magma`. Check with `genor backend list`.
  - Gene location file matching the trait's build (e.g. `NCBI37.3.gene.loc` for GRCh37) and a
    reference panel `g1000_eur` (`.bed/.bim/.fam`) matching the GWAS ancestry, both from the MAGMA
    website (https://cncr.nl/research/magma/).
  - Trait `variant_id` as rsIDs (the MAGMA reference panels are keyed by rsID) and N.
- For enrichment and networks: gene symbols (or Ensembl IDs).

## Steps

1. Loci: `genor clump run runs/bmi --preset loci --out runs/clump_bmi --json`.
2. Positional mapping per locus (region up to 5 Mb; use the trait's build):
   ```bash
   genor lookup genes --region 16:53700000-53900000 --build GRCh37 --json
   ```
   MCP: `lookup_genes_in_region(region="16:53700000-53900000", build="GRCh37")`.
3. Variant consequences for credible-set variants (up to 200 rsIDs per call):
   ```bash
   genor lookup vep rs1421085 rs9939609 --build GRCh38 --json
   ```
   MCP: `lookup_variants(rsids=[...], build="GRCh38")`. Output per variant: `id`,
   `most_severe_consequence`, `location`, `genes` (up to 10 overlapping gene symbols).
4. L2G and existing credible sets (GRCh38 variant ID; see `open-targets-triage`):
   ```bash
   genor lookup variant 16_53767042_T_C --n 10 --json
   ```
5. **Ask the user whether to install MAGMA** (licence acceptance is theirs). If yes, after
   `genor backend add magma /path/to/magma`:
   ```bash
   genor magma gene runs/bmi --gene-loc aux/NCBI37.3.gene.loc --ref-bfile aux/g1000_eur --out runs/magma_bmi --json
   genor magma sets runs/magma_bmi/magma/genes.genes.raw --set-annot aux/gene_sets.txt --out runs/magma_sets_bmi --json
   genor magma sets runs/magma_bmi/magma/genes.genes.raw --covar aux/tissue_expression.txt --out runs/magma_tissue_bmi --json
   ```
   MCP: `run_analysis(method="magma.gene", params={"trait": ..., "gene_loc": ..., "ref_bfile": ...,
   "up_kb": 35, "down_kb": 10, "out": ...})` and `run_analysis(method="magma.sets",
   params={"genes_raw": ..., "set_annot": ... or "covar": ..., "out": ...})`.
6. Enrichment of a prioritised gene list, with an explicit background:
   ```bash
   genor lookup enrich FTO MC4R TMEM18 ... --background background_genes.txt --sources GO:BP,REAC --json
   ```
   MCP: `lookup_enrichment(genes=[...], background=[...], sources=["GO:BP","REAC"])`.
7. Network context:
   ```bash
   genor lookup string FTO MC4R TMEM18 ... --required-score 700 --json
   ```
   MCP: `lookup_string(genes=[...], required_score=700)`.
8. Combine evidence per locus in a table (see [references/locus-to-gene.md](references/locus-to-gene.md))
   and report.

## Defaults and why

| Setting | Default | Why |
|---|---|---|
| MAGMA gene window | 35 kb upstream, 10 kb downstream (`--up-kb`, `--down-kb`) | Captures nearby regulatory variants; widely used (e.g. in FUMA-style analyses and many MAGMA papers). |
| MAGMA gene significance | Bonferroni 0.05 / number of genes (`bonferroni_p` in output) | About 18-19k protein-coding genes gives ~2.7e-6. |
| g:Profiler correction | g:SCS, threshold 0.05 | g:Profiler's method accounts for the hierarchical dependence of GO terms. |
| g:Profiler sources | GO:BP, GO:MF, GO:CC, REAC, KEGG, HPA | Broad coverage; restrict with `--sources` to cut multiple testing. |
| g:Profiler background | All annotated genes unless `--background` given | Usually too broad for GWAS-derived lists; see below. |
| STRING score | 700 | STRING's "high confidence" cut-off (scores run 0-1000). |

Citations: MAGMA (de Leeuw et al. 2015) is recorded in MAGMA run manifests. Lookups do not create
run directories, so cite Ensembl VEP (McLaren et al. 2016), Open Targets (Buniello et al. 2025),
STRING (Szklarczyk et al. 2023) and g:Profiler (Kolberg et al. 2023) by hand, with the date of the
query and the `gprofiler_version` / `string_version` from the output.

## Decision points

- **Ask the user which loci to follow up** (all genome-wide significant ones, or a subset).
- **Ask the user whether to install MAGMA** and confirm they accept its licence. Do not download
  it for them.
- **Ask the user which gene-set or gene-property file** to use (MSigDB collections, GO, tissue
  expression such as GTEx, single-cell atlases). They download it; check its licence and gene ID type.
- **Ask the user which enrichment background** is appropriate: all genes tested by MAGMA, all
  protein-coding genes, or all genes in GWAS loci. The choice changes results.
- **Ask the user how genes were prioritised** before running enrichment, to avoid circularity.

## Reading the output

- `lookup genes`: list of `gene_id`, `symbol`, `biotype`, `start`, `end`, `strand`. Includes
  non-coding genes and pseudogenes; filter by `biotype` deliberately and say so.
- `lookup vep`: `most_severe_consequence` across transcripts (e.g. `missense_variant`,
  `intron_variant`). A coding consequence in the credible set is strong evidence for that gene; an
  intronic one is not evidence for the host gene.
- `magma gene` summary: `n_genes`, `bonferroni_p`, `n_sig`, `top_genes` (GENE ID and P),
  `genes_out`, `genes_raw` (input for `magma sets`). Gene IDs follow the gene-loc file (Entrez IDs
  in the NCBI files); map to symbols before enrichment.
- `magma sets` summary: `n_sets`, `top` (10 lowest P rows from MAGMA's `.gsa.out`), `gsa_out`.
  GeNor does not correct across sets; apply Bonferroni or FDR over `n_sets` yourself.
- `lookup enrich`: `n_terms`, `terms` (up to 100: `source`, `term_id`, `name`, `p_value`
  (g:SCS adjusted), `term_size`, `intersection_size`), `gprofiler_version`, `background`.
- `lookup string`: `string_version`, `n_edges`, `edges` (up to 300), `ppi_enrichment`
  (STRING's test of more interactions than expected for a random set of the same size and degree).

Details: [references/magma.md](references/magma.md),
[references/enrichment-and-networks.md](references/enrichment-and-networks.md).

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `backend_missing` | MAGMA not registered | User downloads MAGMA, then `genor backend add magma /path/to/magma`. |
| `backend_not_found` | Path given to `backend add` does not exist or lacks the executable | Point to the binary or its folder. |
| `backend_failed` | MAGMA exited non-zero; hint has the log tail | Read `magma_annotate.log` / `magma_gene.log` / `magma_sets.log` in the run's `magma/` folder. Common: IDs not matching the reference (non-rsID IDs), build mismatch between gene-loc and positions. |
| `file_not_found` | Gene-loc file or `<ref-bfile>.bed` missing | Download MAGMA auxiliary files; pass the prefix without `.bed`. |
| `n_missing` | Trait has no N | Reload with `--n`. |
| `bad_param` | `magma sets` without `--set-annot` or `--covar` | Give one. |
| `region_too_large` | Region > 5 Mb for `lookup genes` | Split it. |
| `too_many` | > 200 rsIDs to VEP, or > 500 genes to STRING | Batch the calls. |
| `remote_unavailable` | Service down / 5xx after 3 retries | Retry later; report it. The GRCh37 Ensembl mirror is less reliable than GRCh38. |
| `remote_bad_request` | API rejected the query | Check IDs and parameters. |

## What not to conclude

- The nearest gene is not necessarily the causal gene. In benchmark sets of loci with known
  causal genes it is right in only part of cases (estimates vary by study, often around half to
  two-thirds). The FTO obesity variants act on IRX3/IRX5, not FTO.
- A gene in a credible-set region is a candidate, not a finding.
- MAGMA gene p-values are not independent between neighbouring genes; adjacent significant genes
  often reflect one signal.
- Enrichment of a gene list built from GWAS hits is biased towards long, well-annotated,
  well-studied genes. Significance does not show that the pathway is causal.
- STRING networks reflect what has been studied as much as biology. Literature-derived edges can
  reflect the trait itself (circularity).
- Tissue or cell-type enrichment shows that trait-associated genes are specifically expressed
  there, not that the disease process starts there.

## Related skills

- `fine-mapping`, `colocalisation` (QTL link from variant to gene), `open-targets-triage` (L2G,
  tractability), `interpretation-and-reporting`, `gwas-concepts`.
- References: [locus-to-gene.md](references/locus-to-gene.md), [magma.md](references/magma.md),
  [enrichment-and-networks.md](references/enrichment-and-networks.md).
