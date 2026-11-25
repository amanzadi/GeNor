# Where to find GWAS summary statistics

Generated from `src/genor/sumstats/sources.toml` (also `genor sumstats sources`). GeNor never downloads GWAS: you download the file, accept the provider's terms, and give GeNor the local path. Checked 2026-11-25; services change.

## NHGRI-EBI GWAS Catalog (summary statistics) (`gwas-catalog`)

- URL: https://www.ebi.ac.uk/gwas/downloads/summary-statistics
- FTP: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/
- Access: Open. Files by study accession (GCST...). The old summary-statistics REST API and REST API v1 are retired; browse the website or FTP.
- Format: GWAS-SSF (harmonised *.h.tsv.gz files are GRCh38) and submitter formats
- Build: Harmonised files: GRCh38. Check each study.
- Good for: Any trait; first place to look.
- Cite: Cite the original study (PubMed ID on the study page) and the GWAS Catalog (Cerezo et al. 2025, NAR).

## Psychiatric Genomics Consortium (`pgc`)

- URL: https://pgc.unc.edu/for-researchers/download-results/
- Access: Free after agreeing to terms of use on the website (per dataset). Some datasets exclude 23andMe or require a separate request.
- Format: Mostly PLINK/ricopili style text, GRCh37
- Build: GRCh37
- Good for: Psychiatric disorders (MDD, SCZ, BIP, ADHD, ASD, PTSD, ...).
- Cite: Cite the paper listed for each download and follow PGC's acknowledgement rules.

## IEU OpenGWAS (`opengwas`)

- URL: https://opengwas.io
- Access: Needs a free account and a personal JWT token for API use; files can be downloaded as GWAS-VCF per study ID. Token handling stays with the user, never in GeNor.
- Format: GWAS-VCF (GeNor reads single-trait GWAS-VCF directly)
- Build: GRCh37
- Good for: Large harmonised collection, popular for MR.
- Cite: Cite the original study and Elsworth et al. 2020 (bioRxiv 10.1101/2020.08.10.244293).

## FinnGen public results (`finngen`)

- URL: https://www.finngen.fi/en/access_results
- Access: Open after filling a short form; per-release bucket of per-endpoint files.
- Format: Tab-separated, GRCh38, columns #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...
- Build: GRCh38
- Good for: Disease endpoints from Finnish registries; Finnish (EUR, isolate) ancestry.
- Cite: Kurki MI et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature 2023.

## Pan-UK Biobank (`pan-ukb`)

- URL: https://pan.ukbb.broadinstitute.org/downloads
- Access: Open (AWS open data).
- Format: Tab-separated per phenotype with per-ancestry columns (EUR, AFR, AMR, CSA, EAS, MID). Pick one ancestry's columns with --col.
- Build: GRCh37
- Good for: Thousands of UK Biobank phenotypes across six ancestry groups.
- Cite: Pan-UKB team, https://pan.ukbb.broadinstitute.org (2020); cite UK Biobank.

## Neale lab UK Biobank GWAS (`neale-ukb`)

- URL: http://www.nealelab.is/uk-biobank
- Access: Open.
- Format: Tab-separated, GRCh37, variant as chr:pos:ref:alt
- Build: GRCh37
- Good for: Many UK Biobank traits, EUR.
- Cite: Neale lab UK Biobank GWAS round 2 (2018).

## GIANT consortium (anthropometric traits) (`giant`)

- URL: https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
- Access: Open.
- Format: Text, GRCh37
- Build: GRCh37
- Good for: BMI, height, waist-hip ratio.
- Cite: Cite the specific GIANT paper (e.g. Yengo et al. 2018 for BMI/height).

## eQTLGen (blood cis/trans-eQTL) (`eqtlgen`)

- URL: https://www.eqtlgen.org/phase1.html
- Access: Open.
- Format: Text; SMR .besd versions available
- Build: GRCh37
- Good for: Exposure data for SMR/MR on gene expression.
- Cite: Vosa U et al. Nat Genet 2021.

## GTEx portal (eQTL/sQTL) (`gtex`)

- URL: https://gtexportal.org/home/downloads/adult-gtex/qtl
- Access: Significant pairs open; full summary statistics via the GTEx/AnVIL buckets.
- Format: Parquet / text, GRCh38
- Build: GRCh38
- Good for: Tissue QTLs for colocalisation.
- Cite: GTEx Consortium. Science 2020;369:1318-1330.

## eQTL Catalogue (`eqtl-catalogue`)

- URL: https://www.ebi.ac.uk/eqtl/Data_access/
- Access: Open via FTP (tabix-indexed per dataset). The old REST API is retired.
- Format: Tab-separated, tabix-indexed, GRCh38
- Build: GRCh38
- Good for: Uniformly processed QTLs from many studies.
- Cite: Kerimov N et al. Nat Genet 2021.

## PGS Catalog (`pgs-catalog`)

- URL: https://www.pgscatalog.org
- Access: Open.
- Format: Scoring files
- Build: GRCh37/GRCh38 harmonised
- Good for: Published polygenic score weights (not summary statistics).
- Cite: Lambert SA et al. Nat Genet 2021.
