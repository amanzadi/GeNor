# Functional enrichment (g:Profiler) and networks (STRING)

## g:Profiler: `genor lookup enrich`

```bash
genor lookup enrich GENE1 GENE2 ... [--background bg.txt] [--sources GO:BP,REAC] --json
```

GeNor sends an over-representation query to g:Profiler's g:GOSt with:
- organism `hsapiens`;
- significance correction `g_SCS`, threshold 0.05 (results above threshold are not returned);
- sources default `GO:BP, GO:MF, GO:CC, REAC, KEGG, HPA`;
- custom background (`domain_scope=custom`) only if `--background` is given (a text file, one
  gene per line). Otherwise g:Profiler uses all genes annotated in each source.

Output: `n_terms`, up to 100 `terms` (`source`, `term_id`, `name`, `p_value` adjusted, `term_size`,
`intersection_size`), `gprofiler_version` (record it), `background`.

### Gene identifiers

Use gene symbols or Ensembl gene IDs. GeNor does not set g:Profiler's namespace for numeric IDs,
so bare Entrez IDs (as in MAGMA's NCBI output) can be misread; convert them first.

### Choosing the background

The background must be the set of genes that *could* have been selected. Options:

| Gene list came from | Background |
|---|---|
| MAGMA significant genes | All genes in `genes.genes.out` (genes MAGMA tested) |
| Genes in GWAS loci (positional) | All genes with SNP coverage in the GWAS, or all protein-coding genes if you filtered to protein-coding |
| An expression-filtered list | Genes expressed in that tissue |

Using all annotated genes as background for a GWAS-derived list inflates enrichment for
brain-expressed, long, and well-annotated genes. **Ask the user** which background fits, and
report it.

### Biases

- Gene length and SNP density: long genes contain more variants and are more likely to be hit.
  MAGMA corrects for this inside its own gene-set analysis; g:Profiler on a hit list does not.
- Annotation bias: well-studied genes have more GO terms.
- Locus redundancy: several genes from one locus (e.g. a cluster of related genes) can make one
  locus look like a pathway. Collapse to one gene per locus, or report per-locus counts.
- Circularity: if genes were prioritised because of their known function (e.g. "picked the
  neuronal gene in each locus"), enrichment for that function is guaranteed.

Prefer MAGMA's competitive gene-set analysis for genome-wide pathway questions; use g:Profiler
for describing a short prioritised list.

## STRING: `genor lookup string`

```bash
genor lookup string GENE1 GENE2 ... --required-score 700 --json
```

Queries STRING v12 (`version-12-0.string-db.org`) for human (taxon 9606):
- `edges` (up to 300): pairs of genes with combined score >= `required_score` (scale 0-1000;
  150 low, 400 medium, 700 high, 900 highest confidence; GeNor's `score` field is 0-1).
- `ppi_enrichment`: STRING's test of whether the list has more interactions than expected for a
  random gene set of the same size and degree distribution (`number_of_nodes`, `number_of_edges`,
  `expected_number_of_edges`, `p_value`).

Limit: 500 genes per call (`too_many`).

### Interpretation

- STRING scores combine experiments, co-expression, curated databases, and text mining. Text-mined
  edges can link genes simply because they appear in papers about the same disease, which is
  circular for GWAS interpretation.
- A significant PPI enrichment says the genes are more connected than random genes, which is
  expected for well-studied genes. It does not show a shared mechanism.
- Hub genes are central partly because they are studied more.

## Citations (add by hand; lookups do not create manifests)

- Kolberg L, et al. g:Profiler (2023 update). Nucleic Acids Res 2023;51:W207-W212. Report the
  `gprofiler_version`.
- Szklarczyk D, et al. The STRING database in 2023. Nucleic Acids Res 2023;51:D638-D646. Report
  STRING v12.0 and the score cut-off.
- McLaren W, et al. The Ensembl Variant Effect Predictor. Genome Biol 2016;17:122.
