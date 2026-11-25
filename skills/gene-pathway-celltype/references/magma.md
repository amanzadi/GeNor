# MAGMA through GeNor

MAGMA (de Leeuw et al. 2015) aggregates variant associations into gene-level statistics while
accounting for LD from a reference panel, then tests gene sets or continuous gene properties.
GeNor drives a MAGMA binary the user installed; it records the binary path and version in the
manifest.

## Setup (user actions)

1. Download MAGMA (v1.10 in GeNor's registry) from https://cncr.nl/research/magma/ and accept its
   licence (free for academic use, not redistributable). GeNor will not download it.
2. Register it: `genor backend add magma /path/to/magma` (a folder containing `magma` also works).
   Check with `genor backend list` or `genor doctor`.
3. From the same site, download:
   - a gene location file in the GWAS build: `NCBI37.3.gene.loc` for GRCh37, `NCBI38.gene.loc` for
     GRCh38 (Entrez gene IDs; Ensembl-based versions also exist);
   - a reference panel matching the GWAS ancestry, e.g. `g1000_eur.bed/.bim/.fam`
     (1000 Genomes; EUR, EAS, AFR, SAS, AMR available). Pass the prefix without extension.

## What GeNor runs

`genor magma gene TRAIT --gene-loc FILE --ref-bfile PREFIX [--up-kb 35 --down-kb 10]`:

1. Writes `snp.loc` (variant_id, chromosome, position) and `pvals.txt` (SNP, P, N) from the trait.
2. `magma --annotate window=35,10 --snp-loc snp.loc --gene-loc FILE --out annot`
3. `magma --bfile PREFIX --pval pvals.txt use=SNP,P ncol=N --gene-annot annot.genes.annot --out genes`
4. Reads `genes.genes.out` and reports `n_genes`, `bonferroni_p` (0.05 / n_genes), `n_sig`,
   `top_genes` (10 smallest P), and paths to `genes.genes.out` and `genes.genes.raw`.

All files are in `<out>/magma/` with logs `magma_annotate.log` and `magma_gene.log`.

`genor magma sets GENES_RAW --set-annot FILE` runs `magma --gene-results GENES_RAW --set-annot FILE`
(competitive gene-set analysis). `--covar FILE` instead runs `magma --gene-results GENES_RAW
--gene-covar FILE` (gene-property analysis). Output: `n_sets`, `top` (10 smallest P), `gsa_out`
(`<out>/magma/geneset.gsa.out` or `geneprop.gsa.out`). GeNor passes no `--model` options, so
MAGMA's defaults apply.

## Requirements that commonly break runs

- **Variant IDs must be rsIDs** that match the reference `.bim`. Files keyed by `chr:pos` will
  annotate but then match almost nothing in the reference. Check the number of SNPs used in
  `magma_gene.log`.
- **Build**: positions in `snp.loc` come from the GWAS; the gene-loc file must be in the same build.
- **N**: per-variant N if the trait has an `n` column, otherwise the trait-level N. Missing N stops
  with `n_missing`.
- **Ancestry**: the reference panel should match the GWAS ancestry; LD mismatch distorts gene
  p-values.

## Models and their meaning

- Gene analysis with p-value input uses MAGMA's default SNP-wise mean model: it combines the
  variant statistics in the gene window, adjusting for LD. It has most power when many variants
  carry signal.
- Gene-set analysis is **competitive**: it tests whether genes in the set are more associated
  than other genes, not merely associated. By default MAGMA conditions on gene size, gene density
  and minor allele count (and their logs), which addresses the main gene-size bias.
- Gene-property analysis regresses gene association on a continuous covariate per gene (e.g. a
  tissue's expression or specificity). Each column is tested separately by default.
  FUMA-style tissue analyses condition on average expression across tissues
  (`--model condition-hide=Average` in MAGMA); GeNor cannot pass that option yet, so results are
  not directly comparable to FUMA's. Not in GeNor yet. Check the MAGMA manual and the header of
  `.gsa.out` for the test direction your MAGMA version used.

## Gene-set and gene-property files

The user supplies them. Examples:
- Gene sets: MSigDB collections (check licence; MSigDB requires registration), GO, Reactome.
  MAGMA `--set-annot` format is one set per row: set name followed by gene IDs. **Gene IDs must
  match the gene-loc file** (Entrez IDs for the NCBI files). A mismatch yields zero-size sets.
- Gene properties: GTEx tissue expression or specificity tables, single-cell cell-type specificity
  matrices. First column gene ID, one column per tissue or cell type.

## Multiple testing

- Genes: Bonferroni at 0.05 / `n_genes` (GeNor reports it). With ~18-19k genes this is ~2.7e-6.
- Sets and properties: GeNor reports raw P only. Apply Bonferroni over `n_sets` (or over the
  number of tissues / cell types tested) or Benjamini-Hochberg FDR, and state which.

## Not in GeNor yet

- MAGMA gene-property conditioning options and interaction models.
- S-LDSC / LDSC-SEG partitioned heritability for tissues (only `genor ldsc h2` / `rg`).
- Automatic download of auxiliary files or gene sets.

## Reference

de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set analysis of GWAS data.
PLoS Comput Biol 2015;11:e1004219.
