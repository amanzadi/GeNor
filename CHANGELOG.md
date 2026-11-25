# Changelog

## 0.1.0 (unreleased)

First usable release (plan milestones 1–3 plus agent surfaces).

- Summary statistics: format detection (incl. GWAS-VCF, PLINK 2, REGENIE, SAIGE, METAL, FinnGen),
  GWAS-SSF standardisation, build inference, QC, native liftover, plots, sources guide,
  bundled simulated example.
- Meta-analysis (IVW fixed/random, sample-size), LDSC h2/rg (vendored LDSC), PLINK 2 clumping,
  two-sample MR with STROBE-MR checklist, condFDR/conjFDR via the pleioFDR Python port,
  ABF fine-mapping, coloc ABF, MAGMA driver, Open Targets / Ensembl / STRING / g:Profiler lookups,
  Markdown/HTML reports.
- Run manifests, background jobs, refdata and backend registries with pinned checksums.
- CLI with JSON output, MCP server, 12 skills, Claude Code plugin and marketplace.
- Golden tests against TwoSampleMR, coloc, R meta and original LDSC.
