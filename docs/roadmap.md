# Roadmap

Release 0.1 covers milestones 1–3 of the plan (skeleton, sumstats, core statistics) plus the
agent surfaces (skills, MCP server, Claude Code plugin) and parts of the locus and biology
milestones (ABF fine-mapping, coloc ABF, MAGMA driver, Open Targets / Ensembl / STRING /
g:Profiler lookups).

## Next

- **Locus level**: native SuSiE-RSS with LD from the 1000G panel (golden test vs susieR);
  coloc-SuSiE; GCTA-COJO driver (binary already pinned); SMR/HEIDI driver (binary pinned);
  eQTL Catalogue / GTEx v10 lookups over tabix (their REST APIs are retired).
- **Biology**: FUMA-style SNP2GENE (positional + eQTL + L2G) in one command; MAGMA auxiliary
  files via refdata; gene-property (tissue / cell type) inputs.
- **LDSC**: partitioned h2 (baselineLD v2.2) and LDSC-SEG; non-EUR LD scores.
- **Opt-in R backends** (`Rscript` subprocess): susieR cross-check, MR-PRESSO, MVMR, CAUSE, LAVA.
- **MiXeR** through Apptainer/Docker on x86-64 Linux.
- **PRS weights** (PRS-CS, vendored MIT) and PGS Catalog lookup; no genotype scoring.
- **Release**: PyPI (trusted publishing; needs the maintainer's PyPI project), docs site,
  recorded end-to-end agent demo, validated `.mcpb` Desktop extension.

## Release gate

1. `uv run pytest` and `uv run pytest -m slow` pass on Linux x86-64 and macOS arm64.
2. Clean-machine run: `uvx --from . genor doctor` and the example pipeline with no R/conda/Docker.
3. End-to-end on the example plus two public traits downloaded by hand (e.g. PGC MDD and a GWAS
   Catalog BMI study); report inspected by hand.
4. Agent test in Claude Code (plugin) and Claude Desktop (MCP only): "estimate the genetic
   correlation and causal effect of BMI on depression from these two files" completes and stops
   to ask about ancestry, prevalence, sample overlap and large downloads. Codex and opencode:
   smoke run of the same prompt.
