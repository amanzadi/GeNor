# Genetic variation

## Variant types

- **SNP / SNV**: a single base differs between chromosomes (e.g. C vs T). Most GWAS variants.
  "SNP" traditionally implies a common variant; "SNV" any frequency.
- **Indel**: insertion or deletion of one or more bases. Written in VCF style with an anchor base
  (REF `A`, ALT `AT`), or in some GWAS files as `I`/`D`, which cannot be mapped without context.
- **Structural variants and copy-number variants**: larger changes, mostly not in standard GWAS
  summary statistics.
- **Multi-allelic sites**: more than two alleles at one position; GWAS files usually split them
  into separate biallelic rows.

## Alleles in summary statistics

- **Effect allele** (also "tested allele", A1 in some tools): the allele the beta refers to. Beta
  is the change in trait per copy of this allele (additive model).
- **Other allele** (non-effect, A2).
- **REF / ALT**: reference genome allele and alternative allele. Not the same as effect / other.
  Open Targets IDs use REF/ALT on GRCh38.
- Effect allele is not necessarily the minor allele or the risk allele. Flipping the effect allele
  flips the sign of beta (and frequency f becomes 1 - f).

## Strand and palindromic SNPs

DNA has two strands. A/G on one strand is T/C on the other. Files may report alleles on different
strands. Harmonisation (in GeNor's MR, coloc and meta) flips strands by complementing alleles.

**Palindromic SNPs** (A/T or C/G) look the same on both strands, so strand cannot be inferred
from alleles. GeNor uses allele frequencies to resolve them when MAF <= 0.42, and drops them
otherwise (or if frequencies are missing: warning `palindromic_dropped_no_freq`).

## Allele frequency and MAF

- **Effect allele frequency (EAF)**: frequency of the effect allele in the study (or a reference).
- **Minor allele frequency (MAF)**: min(EAF, 1 - EAF). Common: MAF >= 5%; low-frequency 1-5%;
  rare < 1%.
- GeNor's QC removes MAF < 0.01 by default (`--maf`). GWAS power for rare variants is low and
  their statistics less stable.
- Frequencies differ between ancestries, sometimes a lot.

## Identifiers and builds

- **rsID** (dbSNP, e.g. rs1421085): stable-ish name independent of build; some rsIDs merge or
  split between dbSNP releases.
- **Positional IDs** (`chr:pos:ref:alt`): depend on genome build.
- **Genome builds**: GRCh37 (hg19) and GRCh38 (hg38) differ in coordinates. Mixing builds silently
  breaks matching. GeNor infers the build from positions on load (`build_inferred`) and flags
  conflicts (`build_conflict`). `genor sumstats liftover` converts between them using UCSC chain
  files; a few percent of variants may not map.

## Imputation

Most GWAS genotype a few hundred thousand variants on an array and impute millions more using a
reference panel (1000 Genomes, HRC, TOPMed). **INFO** (or R2) measures imputation quality from 0
to 1; GeNor's QC removes INFO < 0.8 by default. Poorly imputed variants have attenuated effects
and unreliable statistics. A variant absent from the reference panel cannot be imputed, so it
cannot appear in fine-mapping credible sets.

## Explaining to a user

"Each row of your file is one variant. Beta tells you how much the trait changes per copy of the
effect allele. If another study used the other allele as the effect allele, the sign flips; that
is why GeNor harmonises alleles before comparing studies."
