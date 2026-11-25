# Using QTL data with genor coloc abf

GeNor does not download QTL data or call eQTL APIs. The user downloads a QTL summary statistics
file and GeNor loads it like any GWAS. Sources listed by `genor sumstats sources`:

| Key | Source | Build | Notes |
|---|---|---|---|
| `gtex` | GTEx portal (eQTL/sQTL), 49 tissues | GRCh38 | Significant pairs are open; full (all-pairs) statistics via GTEx/AnVIL buckets. Coloc needs all pairs, not significant pairs only. |
| `eqtl-catalogue` | eQTL Catalogue (uniformly processed eQTL/sQTL/... datasets) | GRCh38 | Open FTP, tabix-indexed per dataset. The old REST API is retired. |
| `eqtlgen` | eQTLGen (blood cis/trans-eQTL, N ~ 31k) | GRCh37 | Open; some releases report Z-scores and N rather than beta/SE. |

Show details with `genor sumstats sources <key>`.

Open Targets also hosts precomputed GWAS-QTL colocalisations; see "Precomputed colocalisation"
below.

## Preparing one molecular trait

coloc compares two association signals. A QTL file holds thousands of molecular traits (genes,
transcripts, splice junctions). Extract **one** molecular trait for the region before loading:

1. Filter the QTL file to the gene of interest (e.g. eQTL Catalogue `molecular_trait_id` or
   `gene_id`; GTEx `gene_id`), and to the region.
2. Keep variant ID or chromosome/position, effect allele, other allele, beta, SE (or p-value),
   allele frequency, and N.
3. Load it as a quantitative trait, giving N and ancestry:
   ```bash
   genor sumstats load eqtl_GENE_adipose.tsv --trait GENE_adipose_eQTL --ancestry EUR --n 581 --out runs/eqtl_GENE_adipose --no-plots --json
   ```
   Expect warnings such as `few_variants` (a single-region file is not genome-wide; that is fine
   here) and possibly `build_inferred`. Check the build GeNor reports; a region-only file can make
   build inference unreliable, so pass `--build` when you know it.
4. Lift to the GWAS build if needed (`genor sumstats liftover ... --to GRCh37`), or lift the GWAS.

The filtering in step 1 is a plain table operation; GeNor has no command for it. Do it in the
shell or Python and tell the user what you filtered.

## Scale and prior

eQTL effect sizes are usually on a normalised expression scale (sdY about 1). If the file gives
allele frequency and N, GeNor estimates sdY and reports `sdY_estimated`; values far from 1 for
normalised expression suggest wrong N.

## Matching variants

- GTEx and eQTL Catalogue use GRCh38 positions and their own variant IDs (e.g. `chr1_123_A_G_b38`).
  If the GWAS uses rsIDs, GeNor will fall back to chromosome:position matching only when IDs do
  not match; make sure both files have chromosome and position columns.
- `coloc_too_few` most often means a build mismatch or a QTL file that only holds significant
  pairs.
- `allele_mismatch` (serious) above 5% means the alleles do not line up; check build and strand.

## Tissue choice

Ask the user. Guidance:
- Prefer tissues with a plausible link to the trait.
- Larger QTL samples have more power; blood (eQTLGen) is powerful but may miss tissue-specific
  regulation.
- QTLs are often shared across tissues, so a colocalisation in one tissue does not localise the
  effect to that tissue.
- Testing many tissues and genes raises the chance of a spurious high PP.H4. Report the number of
  tests.

## Precomputed colocalisation (Open Targets)

Open Targets runs coloc and eCAVIAR between its fine-mapped GWAS credible sets and QTL credible
sets (eQTL Catalogue, UKB-PPP pQTLs and others) on GRCh38.

```bash
genor lookup variant 16_53767042_T_C --n 10 --json     # credible sets containing the variant -> studyLocusId
genor lookup coloc 007688a660da8613b177ce1ad1db36bb --n 20 --json
```

MCP: `lookup_open_targets(kind="variant", query="16_53767042_T_C")`, then
`lookup_open_targets(kind="coloc", query="<studyLocusId>")`.

Each row has `colocalisationMethod`, `h4` (coloc PP.H4, where the method is coloc), `clpp`
(eCAVIAR colocalisation posterior), and `otherStudyLocus.study` with `studyType` (e.g. `gwas`,
`eqtl`, `pqtl`, `sqtl`) and `traitFromSource`. The page returns the first `--n` rows of a possibly
long list (1,852 for the example above), not necessarily the strongest; raise `--n` when needed.
Open Targets fine-mapping uses SuSiE, so its colocalisations can handle multiple signals that
GeNor's coloc ABF cannot.

## Not in GeNor yet

- SMR/HEIDI (summary-data MR with QTLs): `genor backend install smr` downloads the binary, but no
  GeNor command runs it.
- TWAS / PrediXcan / FUSION.
- coloc-SuSiE.
- Direct eQTL Catalogue or GTEx API queries.
