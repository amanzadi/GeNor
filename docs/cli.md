# CLI reference

Generated from the code by `scripts/gen_cli_docs.py`. Every command accepts `--json` (also the default when stdout is not a terminal). Exit codes: 0 ok, 2 actionable error, 3 a decision for the human is needed, 1 unexpected failure.

### `genor doctor`

Check platform, backends and reference data; say what each missing piece blocks.



### `genor version`

Print the GeNor version.



### `genor example`

Copy the bundled simulated GWAS pair (BMI_sim, MDD_sim) for tests and walkthroughs.

- `--out`, `-o`:  (default: `genor_example`)

### `genor runs`

List run directories under the current directory.



### `genor mcp`

Start the MCP server on stdio (for Claude Desktop, Claude Code, Codex, opencode).



## `genor sumstats`

Load, QC, lift over summary statistics; where to find GWAS.

### `genor sumstats sources`

Where to find GWAS summary statistics (GeNor never downloads them itself).

- `NAME`

### `genor sumstats load`

Parse a summary statistics file (any common format, GWAS-VCF), standardise, QC, plot.

- `FILE`
- `--trait`: Short trait name, e.g. BMI.
- `--out`, `-o`: Run directory to write (created).
- `--build`: GRCh37 or GRCh38 (inferred if omitted).
- `--ancestry`: EUR, EAS, AFR, AMR, SAS, MID, CSA, or multi.
- `--n`: Total sample size if not in the file.
- `--n-cases`: Number of cases (binary traits).
- `--n-controls`: Number of controls (binary traits).
- `--population-prevalence`: For liability-scale h2.
- `--effect-type`: beta or or (odds ratio).
- `--col`: Column override canonical=FILECOL, repeatable.
- `--maf`: Minimum MAF. (default: `0.01`)
- `--info`: Minimum imputation INFO. (default: `0.8`)
- `--drop-mhc`, `--no-drop-mhc`: Remove chr6:25-34 Mb (default: keep and flag).
- `--plots`, `--no-plots`: Write Manhattan and QQ plots. (default: `True`)

### `genor sumstats liftover`

Lift a loaded trait to another genome build (UCSC chain files).

- `TRAIT`
- `--to`: GRCh37 or GRCh38
- `--out`, `-o`: Run directory to write (created).

## `genor meta`

Meta-analysis (IVW fixed/random effects, sample-size weighted).

### `genor meta run`

Meta-analyse two or more loaded traits.

- `TRAITS` (one or more)
- `--out`, `-o`: Run directory to write (created).
- `--scheme`: ivw or samplesize (default: `ivw`)
- `--name`: Trait name for the result. (default: `meta`)

## `genor ldsc`

LD score regression: SNP heritability and genetic correlation.

### `genor ldsc h2`

SNP heritability with LD score regression.

- `TRAIT`
- `--out`, `-o`: Run directory to write (created).
- `--ref`: refdata name of LD scores (default: `ldsc-ldscores-eur`)
- `--w`: refdata name of regression weights (default: `ldsc-weights-eur`)
- `--pop-prev`: Population prevalence (binary traits).
- `--intercept`: Constrain the intercept.

### `genor ldsc rg`

Genetic correlation between two traits with LD score regression.

- `TRAIT1`
- `TRAIT2`
- `--out`, `-o`: Run directory to write (created).
- `--ref`:  (default: `ldsc-ldscores-eur`)
- `--w`:  (default: `ldsc-weights-eur`)
- `--pop-prev1`: Population prevalence of trait1 (binary)
- `--pop-prev2`: Population prevalence of trait2 (binary)

## `genor clump`

LD clumping into independent loci.

### `genor clump run`

Clump a trait into independent index variants.

- `TRAIT`
- `--out`, `-o`: Run directory to write (created).
- `--mode`: plink (LD panel) or distance (no LD; approximate) (default: `plink`)
- `--preset`: loci (r2 0.1, 250 kb) or mr (r2 0.001, 10 Mb) (default: `loci`)
- `--p1`: Index p-value threshold (preset default 5e-8)
- `--r2`: LD r2 threshold (preset default)
- `--kb`: Window in kb (preset default)
- `--panel`: refdata name, e.g. g1000-eur

## `genor mr`

Two-sample Mendelian randomisation with STROBE-MR report.

### `genor mr run`

Instrument selection, harmonisation, IVW/Egger/median/mode, diagnostics, plots, STROBE-MR.

- `--exposure`: Loaded exposure trait directory.
- `--outcome`: Loaded outcome trait directory.
- `--out`, `-o`: Run directory to write (created).
- `--p-threshold`: Instrument p-value threshold. (default: `5e-08`)
- `--clump-mode`: plink or distance (default: `plink`)
- `--panel`: LD panel refdata name (default from ancestry, e.g. g1000-eur)
- `--clump-r2`: Clumping r2 for instruments (default: `0.001`)
- `--clump-kb`: Clumping window (kb) (default: `10000`)
- `--instruments`: File whose first column lists instrument rsIDs.
- `--exclude-mhc`, `--no-exclude-mhc`: Drop instruments in chr6:25-34 Mb (default: keep)
- `--overlap`: Sample overlap stated by the user: none|partial|substantial|unknown
- `--nboot`: Bootstrap replicates (weighted median/mode SE) (default: `1000`)
- `--seed`: Random seed (default: `1`)

## `genor pleio`

condFDR / conjFDR with pleioFDR (Python port).

### `genor pleio run`

condFDR (trait1 conditional on trait2) or conjFDR (shared loci).

- `TRAIT1`
- `TRAIT2`
- `--out`, `-o`: Run directory to write (created).
- `--stat`: condfdr or conjfdr (default: `condfdr`)
- `--fdr`: Threshold (default 0.01 condFDR, 0.05 conjFDR)
- `--randprune-n`:  (default: `500`)
- `--seed`:  (default: `1`)
- `--keep-mhc`, `--no-keep-mhc`: Do not exclude the MHC
- `--prepare-only`, `--no-prepare-only`: Write MAT files and config, don't run
- `--ref-mat`: Existing pleioFDR reference .mat (instead of refdata)
- `--template`: Existing 9545380.ref SNP template (instead of refdata)

## `genor finemap`

Fine-mapping (ABF, single causal variant).

### `genor finemap abf`

Wakefield ABF fine-mapping with a credible set.

- `TRAIT`
- `--out`, `-o`: Run directory to write (created).
- `--region`: chr:start-end
- `--lead`: Lead variant ID; region = lead +/- window
- `--window-kb`:  (default: `500`)
- `--p1`: Prior per variant (default: `0.0001`)
- `--level`: Credible set level (default: `0.95`)

## `genor coloc`

Colocalisation (coloc ABF).

### `genor coloc abf`

Colocalisation of two traits in a region (coloc.abf).

- `TRAIT1`
- `TRAIT2`
- `--out`, `-o`: Run directory to write (created).
- `--region`: 
- `--lead`: 
- `--window-kb`:  (default: `500`)
- `--p1`:  (default: `0.0001`)
- `--p2`:  (default: `0.0001`)
- `--p12`:  (default: `1e-05`)

## `genor magma`

MAGMA gene / gene-set / gene-property analysis (user-installed MAGMA).

### `genor magma gene`

MAGMA gene-based analysis.

- `TRAIT`
- `--gene-loc`: MAGMA gene location file
- `--ref-bfile`: MAGMA reference prefix, e.g. g1000_eur
- `--out`, `-o`: Run directory to write (created).
- `--up-kb`:  (default: `35`)
- `--down-kb`:  (default: `10`)

### `genor magma sets`

MAGMA gene-set or gene-property (cell type / tissue) analysis.

- `GENES_RAW`
- `--out`, `-o`: Run directory to write (created).
- `--set-annot`: Gene-set file (e.g. MSigDB)
- `--covar`: Gene-property file (e.g. tissue/cell-type expression)

## `genor lookup`

Open Targets, Ensembl, STRING and g:Profiler lookups.

### `genor lookup search`

Search Open Targets for a gene, disease, variant or study ID.

- `TERM`
- `--entity`: target | disease | variant | study
- `--n`: Number of hits (default: `5`)

### `genor lookup target`

Open Targets target summary with top associated diseases and tractability.

- `ENSEMBL_ID`
- `--n`:  (default: `10`)

### `genor lookup disease`

Open Targets top targets for a disease (EFO/MONDO ID).

- `EFO_ID`
- `--n`:  (default: `25`)

### `genor lookup variant`

Open Targets credible sets and L2G genes that include a variant.

- `VARIANT_ID`
- `--n`:  (default: `10`)

### `genor lookup coloc`

Open Targets colocalisations for a credible set (studyLocusId).

- `STUDY_LOCUS_ID`
- `--n`:  (default: `20`)

### `genor lookup vep`

Ensembl VEP consequences and genes for up to 200 rsIDs.

- `RSIDS` (one or more)
- `--build`:  (default: `GRCh38`)

### `genor lookup genes`

Genes overlapping a region (Ensembl).

- `--region`: chr:start-end
- `--build`:  (default: `GRCh37`)

### `genor lookup string`

STRING v12 protein interaction network and PPI enrichment for a gene list.

- `GENES` (one or more)
- `--required-score`:  (default: `700`)

### `genor lookup enrich`

g:Profiler enrichment (GO, Reactome, KEGG, HPA).

- `GENES` (one or more)
- `--background`: File with background genes, one per line
- `--sources`: Comma-separated, e.g. GO:BP,REAC

## `genor report`

Build reports from run directories.

### `genor report build`

Combine run directories into report.md + report.html with all citations and warnings.

- `RUN_DIRS` (one or more)
- `--out`, `-o`:  (default: `genor_report.md`)
- `--title`:  (default: `GeNor analysis report`)

## `genor refdata`

Reference data (LD scores, LD panels, chain files).

### `genor refdata list`

Reference datasets, sizes and whether they are downloaded.



### `genor refdata fetch`

Download reference datasets (checksum-verified).

- `NAMES` (one or more)
- `--yes`, `-y`: Confirm large downloads
- `--allow-unpinned`, `--no-allow-unpinned`: Allow datasets without a pinned checksum

### `genor refdata remove`

Delete a downloaded dataset.

- `NAME`

## `genor backend`

External tools (PLINK 2, GCTA, MAGMA, pleioFDR, R).

### `genor backend list`

External tools and how to enable the missing ones.



### `genor backend install`

Download a pinned, redistributable binary (plink2, gcta, smr).

- `NAME`

### `genor backend add`

Register a tool you installed yourself (e.g. MAGMA).

- `NAME`
- `PATH`

### `genor backend remove`

Forget a registered tool.

- `NAME`

## `genor job`

Background jobs for long analyses.

### `genor job submit`

Run any genor command in the background: genor job submit -- ldsc rg A B --out R



### `genor job status`

Job state (running / done / failed / cancelled).

- `JOB_ID`

### `genor job result`

The finished job's result JSON.

- `JOB_ID`

### `genor job logs`

Last lines of the job's stdout/stderr.

- `JOB_ID`
- `--n`:  (default: `50`)

### `genor job cancel`

Stop a running job.

- `JOB_ID`

### `genor job list`

Recent jobs.

- `--limit`:  (default: `20`)

## `genor skills`

Agent skills: list, show, install for Claude/Codex/opencode.

### `genor skills list`

Skills shipped with GeNor.



### `genor skills show`

Print a skill (or a reference file: <skill>/references/<file>.md).

- `TOPIC`

### `genor skills install`

Copy GeNor's skills to .claude/skills or .agents/skills.

- `--agent`: claude | codex | opencode (default: `claude`)
- `--scope`: project | user (default: `project`)
- `--dest`: Explicit target directory
- `--force`, `--no-force`: 
