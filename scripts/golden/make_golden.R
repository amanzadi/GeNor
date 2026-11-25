# Regenerate golden reference values (maintainers only).
# Rscript scripts/golden/make_golden.R   (needs TwoSampleMR, coloc, meta, jsonlite)
suppressPackageStartupMessages({library(TwoSampleMR); library(coloc); library(meta); library(jsonlite)})
g <- "tests/golden"
out <- list(versions = list(R = R.version.string, TwoSampleMR = as.character(packageVersion("TwoSampleMR")),
            coloc = as.character(packageVersion("coloc")), meta = as.character(packageVersion("meta"))))

# --- MR ---------------------------------------------------------------
h <- read.delim(file.path(g, "mr_harmonised.tsv"))
bx <- h$beta_exp; by <- h$beta_out; sx <- h$se_exp; sy <- h$se_out
p <- default_parameters()
ivw <- mr_ivw(bx, by, sx, sy, p); ivwfe <- mr_ivw_fe(bx, by, sx, sy, p)
eg <- mr_egger_regression(bx, by, sx, sy, p)
set.seed(1); wm <- mr_weighted_median(bx, by, sx, sy, p)
set.seed(1); wmo <- mr_weighted_mode(bx, by, sx, sy, p)
set.seed(1); smo <- mr_simple_mode(bx, by, sx, sy, p)
wr <- mr_wald_ratio(bx[1], by[1], sx[1], sy[1], p)
loo1 <- mr_ivw(bx[-1], by[-1], sx[-1], sy[-1], p)
out$mr <- list(
  ivw = ivw[c("b","se","pval","Q","Q_df","Q_pval")], ivw_fe = ivwfe[c("b","se","pval")],
  egger = eg[c("b","se","pval","b_i","se_i","pval_i","Q","Q_df","Q_pval")],
  weighted_median = wm[c("b","se")], weighted_mode = wmo[c("b","se")], simple_mode = smo[c("b","se")],
  wald_first = wr[c("b","se","pval")], loo_drop_first = loo1[c("b","se")])

# --- coloc ------------------------------------------------------------
r <- read.delim(file.path(g, "coloc_region.tsv"))
d1 <- list(beta = r$beta1, varbeta = r$se1^2, snp = r$variant_id, MAF = r$maf, N = r$n1[1], type = "quant", sdY = 1)
d2 <- list(beta = r$beta2, varbeta = r$se2^2, snp = r$variant_id, MAF = r$maf, N = r$n2[1], type = "cc", s = 0.4)
cr <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
d1e <- d1; d1e$sdY <- NULL
cre <- coloc.abf(d1e, d2)
fm <- finemap.abf(d1, p1 = 1e-4)
out$coloc <- list(summary = as.list(cr$summary), summary_sdy_estimated = as.list(cre$summary),
                  snp_pp_h4_top = cr$results$SNP.PP.H4[order(-cr$results$SNP.PP.H4)][1:5],
                  finemap_top_pp = fm$SNP.PP[order(-fm$SNP.PP)][1:5],
                  finemap_null_pp = fm$SNP.PP[fm$snp == "null"])

# --- meta -------------------------------------------------------------
m <- read.delim(file.path(g, "meta_studies.tsv"))
res <- lapply(split(m, m$variant_id), function(x) {
  mg <- metagen(TE = x$beta, seTE = x$se, method.tau = "DL", common = TRUE, random = TRUE)
  data.frame(variant_id = x$variant_id[1], beta = mg$TE.common, se = mg$seTE.common,
             beta_random = mg$TE.random, se_random = mg$seTE.random, Q = mg$Q, I2 = mg$I2, tau2 = mg$tau2)
})
out$meta <- do.call(rbind, res)
write(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = TRUE), file.path(g, "expected_r.json"))
cat("ok\n")
# LDSC reference values (tests/golden/ldsc_expected.json) come from the original
# LDSC (CBIIT/ldsc ldsc313 @ 128f996) via scripts/golden/ldsc_reference.py on
# tests/golden/ldsc_{bmi,mdd}.sumstats.gz, which were made by munge_sumstats.py
# from `genor example` files.
