# Formats and column mapping

## What GeNor reads

- Delimited text: tab, comma, whitespace, `;` or `|` (sniffed from the first lines). Leading
  `##` comment lines are skipped. Missing-value tokens: `NA`, `NaN`, `.`, empty, `-`.
- Compression: `.gz`, `.bgz`, `.zip` (exactly one file inside).
- Parquet (`.parquet`).
- Single-trait GWAS-VCF (IEU OpenGWAS format): FORMAT fields `ES` (beta), `SE`, `LP` (-log10 P),
  `AF`, `SS` (N), `ID`. ALT is the effect allele. Multi-allelic rows are skipped. Multi-trait
  GWAS-VCFs are refused (`gwas_vcf_multi`); sites-only VCFs are allowed; VCFs with `GT`/`DS`/`GP`
  sample columns are refused as individual-level data.

## Canonical names (use these on the left of `--col`)

`variant_id`, `chromosome`, `base_pair_location`, `effect_allele`, `other_allele`, `beta`,
`odds_ratio`, `standard_error`, `p_value`, `neg_log10_p`, `z`, `n`, `n_cases`, `n_controls`,
`effect_allele_frequency`, `info`.

`--col canonical=FILECOLUMN` uses the exact header text (case and punctuation as in the file).
Overrides win over auto-detection. MCP: `columns={"canonical": "FILECOLUMN"}`.

## Auto-detection

Headers are lower-cased and stripped to letters and digits, then matched against aliases
(`src/genor/sumstats/columns.py`). The first header column that matches a canonical name wins;
later columns that would map to the same name are ignored. Examples of recognised aliases:

| Canonical | Some aliases |
|---|---|
| variant_id | SNP, rsid, MarkerName, ID, variant, hm_rsid |
| chromosome | CHR, chrom, #CHROM, hm_chrom |
| base_pair_location | BP, POS, position, bp_hg19, hm_pos |
| effect_allele | A1, EA, Allele1, ALT, tested_allele, effect_allele |
| other_allele | A2, OA, Allele2, ALLELE0, NEA, REF |
| beta | BETA, b, Effect, logOR, est |
| odds_ratio | OR, odds_ratio |
| standard_error | SE, StdErr, sebeta, LOG(OR)_SE |
| p_value | P, pval, P-value, P_BOLT_LMM, P_BOLT_LMM_INF |
| neg_log10_p | LOG10P, mlogp, neg_log10_pval, LP |
| effect_allele_frequency | EAF, FRQ, Freq1, A1FREQ, A1_FREQ, af_alt, AF |
| n | N, N_total, OBS_CT, Neff, Weight (METAL), NMISS |
| n_cases / n_controls | Ncases, N_CAS / Ncontrols, N_CON |
| z | Z, zscore, Z_STAT, T_STAT |
| info | INFO, Rsq, R2 |

Notes GeNor emits as `column_note`: A1 is assumed to be the effect allele (true for PLINK, LDSC,
BOLT-LMM, METAL `Allele1`); `Rsq` is treated as imputation quality; `Neff` is used as N;
METAL `Weight` is used as N.

## Format recipes (check the header; versions differ)

- **PLINK 1.9 `.assoc.linear/.logistic`**: A1 is the effect allele, but there is no A2 column in
  older outputs. Add the other allele from the `.bim` before loading, or GeNor stops with
  `columns_not_detected`.
- **PLINK 2 `.glm.*`**: effects are for `A1`, which can be REF or ALT. GeNor detects the format
  (A1 + REF + ALT columns), uses `A1` as the effect allele, takes the other allele as whichever of
  REF/ALT is not A1, and uses `A1_FREQ` as its frequency (a `column_note` says so).
  `.glm.logistic` reports `OR`; `.glm.logistic.hybrid` may report `BETA` on the log-odds scale.
- **REGENIE**: `ALLELE1` is the effect allele, `ALLELE0` the other; `LOG10P` is -log10 P.
  Auto-detection handles this.
- **SAIGE**: `BETA` refers to `Allele2`, not `Allele1`. GeNor detects SAIGE output (an
  `AF_Allele2` or `AC_Allele2` column) and maps `Allele2` as the effect allele and `AF_Allele2` as
  its frequency (a `column_note` says so). Check the note if the file was reformatted.
- **BOLT-LMM**: `ALLELE1` is the effect allele; P columns `P_BOLT_LMM_INF` / `P_BOLT_LMM`. When
  both exist, the first in the header is used; name the one you want with `--col p_value=...`.
- **METAL**: `Allele1` is the effect allele, alleles often lower-case (GeNor upper-cases them);
  `Weight` is N in the sample-size scheme.
- **GWAS Catalog harmonised (`*.h.tsv.gz`)**: GWAS-SSF names, GRCh38. Lift to GRCh37 for GeNor's
  LD references unless rsIDs are present for LDSC/condFDR.
- **FinnGen**: `#chrom pos ref alt rsids ... pval mlogp beta sebeta af_alt`; ALT is the effect
  allele; GRCh38; no N column, so pass `--n-cases` and `--n-controls` from the endpoint page.
- **Pan-UKB**: one column per ancestry (for example `beta_EUR`, `se_EUR`, `af_EUR`, and a P or
  -log10 P column per ancestry). Pick one ancestry with `--col beta=beta_EUR --col
  standard_error=se_EUR ...` and set `--ancestry` to match. Check whether the P column is
  -log10 P (`neg_log10_p`) or P (`p_value`).
- **Variant IDs like `1:12345:A:G`** with no CHR/POS columns: positions are parsed from the ID
  (`position_from_id`). If there are no allele columns either, alleles are taken as REF:ALT with ALT
  as the effect allele (`alleles_from_id`, serious). Confirm against the README.
