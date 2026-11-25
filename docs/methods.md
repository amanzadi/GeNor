# Methods, defaults and validation status

Each method lists its implementation, key defaults (with the reason), and how it is validated.
The skills in `skills/` explain interpretation and the decisions left to the user.

## Summary statistics

**`genor sumstats load`**: reads TSV/CSV/space-delimited text (optionally gzip/zip), Parquet and
single-trait GWAS-VCF. Columns are detected from ~150 aliases (`src/genor/sumstats/columns.py`),
including PLINK 2 `.glm` (A1 is the effect allele), REGENIE, SAIGE (Allele2 is the effect allele),
BOLT-LMM, METAL, FinnGen and GWAS Catalog harmonised files. `--col canonical=FILECOLUMN` overrides.

- Odds ratios become log-odds; a `beta` column that looks like odds ratios stops with a question.
- Z from beta/SE; SE from beta and P when SE is missing; P from Z when P is missing.
- Build inferred from ~4,000 marker SNPs with known GRCh37/GRCh38 positions.
- Individual-level files (`.bed/.pgen/.bgen/...`, VCFs with genotypes) are refused.

QC defaults: MAF ≥ 0.01 (LD panels of ~500 people do not tag rarer variants), INFO ≥ 0.8,
valid P/SE/beta/alleles, duplicates resolved by smallest P. Palindromic SNPs and the MHC are
flagged, not removed. Reports λGC, λ1000, mean χ², genome-wide significant count.

**Liftover**: native UCSC chain-file liftover for SNP positions (no external tool). Variants
mapping to the reverse strand or another chromosome are dropped and counted.

Validation: unit tests for each format, the refusal path, build inference, QC counts, allele
harmonisation (aligned, swapped, strand-flipped, palindromic resolved by frequency, dropped).

## Meta-analysis (`genor meta run`)

Fixed-effect inverse-variance weighting (METAL `SCHEME STDERR`) with Cochran's Q, I² and
DerSimonian-Laird random effects; or sample-size-weighted Z (METAL `SCHEME SAMPLESIZE`).
Alleles are aligned to the first study reporting each variant.

Validation: identical to R `meta::metagen` (fixed, random, Q, I², τ²) to 1e-8 relative on a
60-variant, 3-study fixture.

## LD score regression (`genor ldsc h2`, `genor ldsc rg`)

Regression code vendored unmodified from [CBIIT/ldsc](https://github.com/CBIIT/ldsc) (Python 3
port of bulik/ldsc, GPL-3). Data handling reproduces `munge_sumstats.py` (SNVs only, no strand-
ambiguous SNPs, INFO > 0.9, MAF > 0.01, N ≥ N90/1.5) and `ldsc.py` (200 blocks; two-step
estimator at χ² 30 for h2; reference-panel SNP order). Reference: 1000G Phase 3 EUR HapMap3 LD
scores and weights from Zenodo 10515792 (CC BY 4.0).

Liability-scale h2 uses the sample prevalence from case/control counts and the user's
population prevalence. Warnings: mean χ² < 1.02, h2 z < 4, intercept > 1.1 or < 0.95,
rg outside [-1.2, 1.2], cross-trait intercept ≠ 0 (sample overlap).

Validation: identical to the original `ldsc.py --h2/--rg` (h2, rg, SEs, intercepts to 1e-6
relative) on the same munged input (`tests/test_ldsc_golden.py`, marked slow: downloads 25 MB).

## Clumping (`genor clump run`)

PLINK 2 `--clump` per chromosome against 1000G Phase 3 EUR or EAS. Presets: `loci` (p 5e-8,
r² 0.1, 250 kb), `mr` (p 5e-8, r² 0.001, 10 Mb). `--mode distance` keeps the strongest variant
per window without LD, always with a `distance_clumping` warning.

## Mendelian randomisation (`genor mr run`)

Instrument selection (p < 5e-8, `mr` clumping), harmonisation (TwoSampleMR action 2:
palindromes with MAF > 0.42 dropped), then: IVW (multiplicative random effects, main), IVW fixed
effects, MR-Egger with intercept and I²_GX, weighted median (bootstrap SE), weighted and simple
mode (bootstrap SE), Wald ratios, leave-one-out IVW, Cochran's Q (IVW and Egger), Steiger
directionality, F statistics. Plots: scatter, forest, funnel, leave-one-out. STROBE-MR checklist
with each of the 20 items marked filled / partial / human.

The analysis asks the user about sample overlap unless `--overlap` is given.

Validation against TwoSampleMR 0.7.14 on a 141-instrument fixture: IVW, IVW-FE, MR-Egger
(slope, intercept, SEs, p, Q), Wald ratio and leave-one-out identical to 1e-8; weighted median
and weighted/simple mode point estimates identical to 1e-6 (the mode uses a port of R's
`density()`); bootstrap SEs agree within Monte Carlo error.

## condFDR / conjFDR (`genor pleio run`)

GeNor aligns each trait to pleioFDR's 9,545,380-SNP template (z signed to the template A1),
writes `logpvec`/`zvec` MAT files and `config.txt`, and runs the
[pleioFDR Python port](https://github.com/espenhgn/pleiofdr/tree/python-port) (commit `45e03b9`),
which is itself validated against the MATLAB version. Defaults: condFDR threshold 0.01, conjFDR
0.05, 500 random prunings, MHC excluded.

Validation: a round trip on the chr21 demo data (GeNor alignment of re-encoded traits, half with
swapped alleles, vs pleioFDR on the original MAT files) gives identical result tables
(`scripts/golden/pleio_roundtrip.py`).

## Fine-mapping and colocalisation

`genor finemap abf`: Wakefield approximate Bayes factors, one causal variant, prior 1e-4 per
variant, prior effect SD 0.15·sdY (quantitative, sdY estimated from MAF and N as in coloc) or 0.2
(case-control); credible set conditional on a causal variant being present.

`genor coloc abf`: coloc ABF with p1 = p2 = 1e-4, p12 = 1e-5; reports PP.H0–H4 and per-SNP
PP.H4.

Validation: identical to `coloc::coloc.abf` and `coloc::finemap.abf` (coloc 5.2.3) to 1e-6,
both with supplied and estimated sdY.

## Gene-level and lookups

- MAGMA gene / gene-set / gene-property: drives a user-installed MAGMA (annotation window
  35 kb up / 10 kb down by default).
- Open Targets Platform GraphQL: search, target, disease, variant credible sets with L2G,
  colocalisations of a credible set.
- Ensembl: VEP consequences (GRCh38 or GRCh37 server), genes overlapping a region.
- STRING v12 (version-pinned URL, ≥ 1 s between requests): network and PPI enrichment.
- g:Profiler: GO/Reactome/KEGG/HPA enrichment with g:SCS correction and optional custom background.

Validation: offline tests with recorded responses (`tests/fixtures/http`), and a nightly live job
that checks the services still answer in the expected shape.
