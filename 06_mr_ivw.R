# ============================================================
# 06_mr_ivw.R
# Two-sample Mendelian randomization (IVW) of 5 trait exposures
# on Male Pattern Baldness (outcome).
# ============================================================
 
if (!require("data.table")) install.packages("data.table")
library(data.table)
 
# Per-exposure P-value thresholds 
THRESH <- list(
  Insomnia   = 5e-8,
  Cannabis   = 5e-6,    # relaxed: only 3 SNPs reach 5e-8
  Income     = 5e-8,
  Chronotype = 5e-8,
  ADHD       = 5e-8
)
THIN_WINDOW <- 100       
 
# Helper to load one munged sumstats file
load_munged <- function(path, label) {
  d <- fread(path)
  d <- d[!is.na(Z) & !is.na(N) & A1 != "" & A2 != ""]
  d[, orig_idx := .I]
  d[, P    := 2 * pnorm(-abs(Z))]
  d[, BETA := Z / sqrt(N)]
  d[, SE   := 1 / sqrt(N)]
  setnames(d, c("A1","A2","Z","N","P","BETA","SE"),
              paste0(c("A1","A2","Z","N","P","BETA","SE"), "_", label))
  d
}
 
# Greedy LD-thinning 
greedy_thin <- function(d, label, window = THIN_WINDOW) {
  z_col <- paste0("Z_", label)
  d <- d[order(-abs(get(z_col)))]
  d[!duplicated(orig_idx %/% window)]
}
harmonize <- function(exp_d, out_d, exp_label, out_label) {
  m <- merge(exp_d, out_d, by = "SNP")
  a1e <- m[[paste0("A1_", exp_label)]]; a2e <- m[[paste0("A2_", exp_label)]]
  a1o <- m[[paste0("A1_", out_label)]]; a2o <- m[[paste0("A2_", out_label)]]
  same <- a1e == a1o & a2e == a2o
  swap <- a1e == a2o & a2e == a1o
  pal  <- (a1e == "A" & a2e == "T") | (a1e == "T" & a2e == "A") |
          (a1e == "G" & a2e == "C") | (a1e == "C" & a2e == "G")
  keep <- (same | swap) & !pal
  m <- m[keep]
  swap_kept <- swap[keep]
  beta_o <- paste0("BETA_", out_label); z_o <- paste0("Z_", out_label)
  m[swap_kept, (beta_o) := -get(beta_o)]
  m[swap_kept, (z_o)    := -get(z_o)]
  m
}
 
# IVW estimator 
ivw <- function(bx, by, sy) {
  w     <- 1 / sy^2
  bhat  <- sum(w * bx * by) / sum(w * bx^2)
  se    <- sqrt(1 / sum(w * bx^2))
  Q     <- sum(w * (by - bhat * bx)^2)
  df    <- length(bx) - 1
  pQ    <- pchisq(Q, df, lower.tail = FALSE)
  se_re <- if (df > 0) se * max(1, sqrt(Q / df)) else se
  p     <- 2 * pnorm(-abs(bhat / se_re))
  list(beta = bhat, se = se_re, p = p, Q = Q, df = df, p_Q = pQ)
}
 
#Run MR for each exposure
out_d <- load_munged("Baldness_sumstats.gz", "Baldness")
cat(sprintf("Outcome Baldness: %d SNPs\n", nrow(out_d)))
 
results <- list()
for (exp in names(THRESH)) {
  p_thresh <- THRESH[[exp]]
  cat(sprintf("\n=== %s -> Baldness (P < %.0e) ===\n", exp, p_thresh))
  exp_d <- load_munged(paste0(exp, "_sumstats.gz"), exp)
  
  instr <- exp_d[get(paste0("P_", exp)) < p_thresh]
  cat(sprintf("  %d instruments at this threshold\n", nrow(instr)))
  if (nrow(instr) < 2) { cat("  too few instruments, skipping\n"); next }
  
  instr <- greedy_thin(instr, exp)
  cat(sprintf("  %d after LD-thinning\n", nrow(instr)))
  
  m <- harmonize(instr, out_d, exp, "Baldness")
  cat(sprintf("  %d after allele harmonization\n", nrow(m)))
  if (nrow(m) < 2) { cat("  too few SNPs after harmonization\n"); next }
  
  bx <- m[[paste0("BETA_", exp)]]
  by <- m[["BETA_Baldness"]]
  sy <- m[["SE_Baldness"]]
  r  <- ivw(bx, by, sy)
  
  cat(sprintf("  IVW beta = %+.4f   SE = %.4f   P = %.3g\n", r$beta, r$se, r$p))
  cat(sprintf("  Cochran Q = %.1f (df=%d), P_Q = %.3g\n", r$Q, r$df, r$p_Q))
  
  results[[exp]] <- data.table(
    exposure = exp, outcome = "Baldness",
    P_threshold = p_thresh,
    n_instruments = nrow(m),
    beta_IVW = r$beta, SE_IVW = r$se, P_IVW = r$p,
    Q = r$Q, df = r$df, P_Q = r$p_Q
  )
}
 
# Save the results table
res <- rbindlist(results)
fwrite(res, "MR_IVW_results.tsv", sep = "\t")
cat("\n=== Summary table ===\n")
print(res)
 
# Forest plot
res <- res[order(beta_IVW)]
png("MR_forest.png", width = 1600, height = 700, res = 180)
par(mar = c(5, 9, 4, 5))
y <- seq_len(nrow(res))
xrng <- range(res$beta_IVW - 2 * res$SE_IVW,
              res$beta_IVW + 2 * res$SE_IVW)
plot(res$beta_IVW, y, xlim = xrng,
     ylim = c(0.5, nrow(res) + 0.5), pch = 19, cex = 1.5,
     col = ifelse(res$P_IVW < 0.05, "red", "gray40"),
     yaxt = "n", xlab = "IVW causal estimate (per-SD effect on Baldness)",
     ylab = "", main = "Mendelian randomization (IVW): trait \u2192 MPB")
abline(v = 0, lty = 2)
arrows(res$beta_IVW - 1.96 * res$SE_IVW, y,
       res$beta_IVW + 1.96 * res$SE_IVW, y,
       length = 0.05, angle = 90, code = 3,
       col = ifelse(res$P_IVW < 0.05, "red", "gray40"), lwd = 2)
# Y-axis labels include relaxed-threshold flag for Cannabis
ylabs <- paste0(res$exposure, "  (k=", res$n_instruments,
                ifelse(res$P_threshold < 5e-8, "", ""),
                ifelse(res$P_threshold > 5e-8, "*", ""), ")")
axis(2, at = y, labels = ylabs, las = 1, cex.axis = 0.9)
text(par("usr")[2], y, labels = sprintf("P=%.3g%s", res$P_IVW,
                                        ifelse(res$P_IVW < 0.05, "*", "")),
     pos = 2, cex = 0.85,
     col = ifelse(res$P_IVW < 0.05, "red", "gray30"))
mtext("* Cannabis used relaxed P<5e-6 threshold (only 3 SNPs at 5e-8)",
      side = 1, line = 4, cex = 0.75, adj = 0)
dev.off()
cat("\nSaved MR_IVW_results.tsv and MR_forest.png\n")
