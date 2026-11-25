# condFDR / conjFDR: method and interpretation

## Definitions

- Unconditional (tail-area) FDR for trait 1 at p-value p1: FDR(p1) = pi0 x p1 / F(p1), with F the
  empirical cumulative distribution of trait 1 p-values.
- Conditional FDR: condFDR(p1 | p2) = pi0(p2) x p1 / F(p1 | p2), the FDR for trait 1 among SNPs
  whose trait 2 p-value is at or below p2. If SNPs that are associated with trait 2 are enriched for
  trait 1 association, condFDR is smaller than FDR for those SNPs and more trait 1 loci pass the
  threshold (Andreassen et al. 2013, PLoS Genet 9:e1003455).
- Conjunctional FDR: conjFDR(p1, p2) = max(condFDR(p1 | p2), condFDR(p2 | p1)), a conservative
  estimate of the FDR for association with both traits (Smeland et al. 2020, Hum Genet 139:85).

## What GeNor does and what pleioFDR does

GeNor aligns each trait to the pleioFDR 9,545,380-SNP template by rsID, signs z-scores to the
template's A1 allele (swapped alleles flip sign, strand complements accepted, other mismatches set
to missing), writes `logpvec` / `zvec` MAT files and a `config.txt`, and calls
`pleiofdr.cli.run_from_config`. The FDR estimation, random pruning, locus definition, tables and
plots are pleioFDR's own code. The reference (`pleiofdr-ref-eur`) holds the LD structure
(r2 > 0.1, 1000G Phase 3 EUR) used for random pruning and locus clumping.

Config written by GeNor (see `<out>/pleiofdr/config.txt`): `stattype`, `fdrthresh`,
`randprune_n`, `exclude_chr_pos` (MHC unless `--keep-mhc`), `onscreen=false`, trait files and
names, reference paths, output directory. Other pleioFDR settings stay at the port's defaults;
GeNor does not expose them.

## Random pruning

LD makes neighbouring SNPs' statistics dependent, which inflates enrichment estimates. pleioFDR
repeatedly selects one random SNP per LD block (r2 > 0.1) and averages the FDR estimates over
`randprune_n` iterations. 500 is the usual published setting; fewer iterations give noisier
FDR values (`pleio_few_prunes` below 100).

## Excluded regions

GeNor excludes the MHC by default. Some published analyses also exclude other long-range LD
regions (for example the 8p23.1 inversion or APOE for neurodegenerative traits). GeNor does not
do this; if the user's field expects it, say it was not done.

## Interpreting results

- Read the conditional Q-Q plots in `results/`: leftward deflection of trait 1 curves with
  stricter trait 2 thresholds shows cross-trait enrichment. Without enrichment, condFDR adds little.
- condFDR loci for trait 1 that were not genome-wide significant in trait 1 alone are "novel"
  only relative to that GWAS; check the GWAS Catalog / Open Targets before calling them new.
- conjFDR loci: report the lead SNP, both z-score signs (concordant / discordant), and whether
  each trait reached genome-wide significance on its own.
- Mixed effect directions can coexist with rg near zero. conjFDR is often used for exactly that
  situation; MiXeR quantifies polygenic overlap formally but is not wired into GeNor yet.

## Pitfalls

- **Sample overlap / shared controls** inflate cross-trait enrichment and shared loci. pleioFDR
  as called by GeNor applies no correction. Ask the user and report it.
- **Low template overlap**: non-rsID IDs or a filtered file reduce power; check `matched`.
- **Different ancestries**: the reference is EUR only.
- **Memory**: the reference needs about 12 GB RAM; a killed job usually means out of memory.
- **`--prepare-only`** needs the template (262 MB) but not the 2.4 GB reference: it writes the
  aligned MAT files and `config.txt` (pointing where the reference will be) and skips the run.
- Users who already have the pleioFDR files (e.g. on a cluster) can pass `--template` and
  `--ref-mat` instead of downloading them through GeNor.
