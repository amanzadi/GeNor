# What Open Targets fields mean

GeNor requests a small, fixed set of fields (see `src/genor/lookup/services.py`). Anything not
listed here is not fetched; tell the user to use the Open Targets web interface or API directly for
more.

## `lookup variant` (credible sets and L2G)

| Field | Meaning |
|---|---|
| `id`, `rsIds` | Open Targets variant ID and rsIDs |
| `credibleSets.count` | Number of credible sets (across all studies) containing the variant |
| `studyLocusId` | ID of one credible set (one signal in one study). Input to `lookup coloc` |
| `pValueMantissa`, `pValueExponent` | Lead p-value of the credible set = mantissa x 10^exponent |
| `study.id` | GWAS Catalog accession (GCST...), FinnGen (FINNGEN_R...), UK Biobank or QTL study IDs |
| `study.traitFromSource` | Trait label as reported by the study |
| `study.studyType` | `gwas` or a molecular QTL type (e.g. `eqtl`, `pqtl`, `sqtl`) |
| `l2GPredictions` | Top 3 genes by L2G score for this credible set |

The variant's own posterior inclusion probability within each set is not requested by GeNor.
Credible sets come from SuSiE fine-mapping (or summary-statistic based methods for studies without
full data); a variant can be in a credible set with low PIP.

### L2G (locus-to-gene)

A gradient-boosting model that scores every gene within a window of a credible set using features
such as distance to the credible-set variants (weighted by PIP), molecular QTL colocalisation,
chromatin interaction and enhancer-gene links, variant consequence (VEP), and gene density. It is
trained on gold-standard GWAS loci with a known causal gene (often drug targets or Mendelian
genes). Interpretation:

- 0.5 and above: commonly treated as a confident gene assignment.
- Scores are per credible set; compare within a locus, not across loci.
- Gold-standard bias: well-studied, coding-variant-driven loci are over-represented in training.

## `lookup coloc`

| Field | Meaning |
|---|---|
| `colocalisationMethod` | Method that produced the row |
| `h4` | coloc posterior for a shared causal variant (for coloc rows) |
| `clpp` | eCAVIAR colocalisation posterior probability: sum over variants of the product of PIPs (for eCAVIAR rows) |
| `otherStudyLocus.study` | The other study's ID, `studyType`, and `traitFromSource` |

`h4` and `clpp` are on different scales: CLPP values are much smaller than PP.H4 for the same
locus. Do not threshold CLPP at 0.8.

## `lookup target`

| Field | Meaning |
|---|---|
| `approvedSymbol`, `approvedName`, `biotype` | HGNC symbol, name, Ensembl biotype |
| `functionDescriptions` | First two UniProt function descriptions |
| `tractability` | Entries with value true: `modality` (`SM`, `AB`, `PR`, `OC`) and `label` |
| `associatedDiseases.count` / `rows` | Diseases ranked by overall association score |

Tractability buckets come from the Open Targets tractability pipeline: e.g. for small molecules,
"Approved Drug", "Advanced Clinical", "Phase 1 Clinical", "Structure with Ligand",
"High-Quality Ligand", "High-Quality Pocket", "Med-Quality Pocket", "Druggable Family"; for
antibodies, clinical labels and subcellular location labels (e.g. "UniProt loc high conf",
"GO CC high conf"). Clinical labels mean a drug exists for the target, for some indication.

## `lookup disease`

`associatedTargets.rows`: targets with an overall `score` (0-1) for the disease.

## Association scores

The overall target-disease score aggregates evidence across data sources (GWAS credible sets,
gene burden, ClinVar/Gene2Phenotype and other rare-disease sources, somatic mutations, ChEMBL
known drugs, pathways, literature co-mention, animal models, expression changes). Each source is
scored, weighted, and combined with a harmonic sum, so many weak items add up less than one
strong item.

- Not a probability, not calibrated across diseases.
- GeNor fetches only the overall score, not the per-datatype breakdown. A high overall score can
  come entirely from literature mining or known drugs, with no genetic evidence. For genetic
  support, use `lookup variant` evidence instead.

## Data release

Open Targets releases change content and sometimes the schema. GeNor does not record the release.
Note the date of the query in the report.
