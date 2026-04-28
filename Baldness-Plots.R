# ============================================================
# baldness_plots.R
# Manhattan + QQ plot for the MPB GWAS (Baldness.tab).
# Run this in RStudio — open this file and click "Source"
# (or press Cmd+Shift+Enter on Mac).
# ============================================================

if (!require("qqman"))      install.packages("qqman")
if (!require("data.table")) install.packages("data.table")
library(qqman)
library(data.table)

# Load 
gwas <- fread("Baldness.tab")
gwas[, CHR := suppressWarnings(as.numeric(CHR))]
gwas[, P_BOLT_LMM := suppressWarnings(as.numeric(P_BOLT_LMM))]
gwas <- gwas[!is.na(CHR) & CHR %in% 1:23 &
               !is.na(P_BOLT_LMM) & P_BOLT_LMM > 0 & P_BOLT_LMM <= 1]
plot_data <- gwas[, .(SNP = SNP, CHR = CHR, BP = BP, P = P_BOLT_LMM)]
cat("After cleaning:", nrow(plot_data), "SNPs\n")

# Manhattan plot
png("Baldness_manhattan.png", width = 1400, height = 550, res = 130)
manhattan(plot_data,
          main = "Male Pattern Baldness GWAS (N = 205,327)",
          ylim = c(0, 50),
          chrlabs = c(as.character(1:22), "X"),
          col = c("#1f77b4", "#7f7f7f"),
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8))
dev.off()

# QQ plot
png("Baldness_qq.png", width = 600, height = 600, res = 130)
qq(plot_data$P, main = "QQ plot: MPB GWAS")
dev.off()

# Lambda + top hits
chisq <- qchisq(1 - plot_data$P, df = 1)
cat(sprintf("\nlambda_GC = %.3f\nmean chi^2 = %.3f\n",
            median(chisq) / qchisq(0.5, df = 1), mean(chisq)))
print(plot_data[order(P)][1:10])
fwrite(plot_data[order(P)][1:20], "Baldness_top20_SNPs.tsv", sep = "\t")
cat("\nDone. PNGs are in your working directory.\n")

