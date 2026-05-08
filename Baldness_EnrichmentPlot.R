library(ggplot2)
library(dplyr)
library(scales)

# =============================================================================
# Baldness Partitioned Heritability Enrichment
# GenoSkylinePlus Tier 3 + Baseline annotations
# Bonferroni threshold: 0.05 / 120 = 4.17e-4
# =============================================================================

results_file <- "/Users/jenniferbrant/Desktop/stat620/results/enrichment/Baldness_enrichment.results"
outdir       <- "/Users/jenniferbrant/Desktop/stat620/results/enrichment"

df_raw <- read.table(results_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# --------------------------------------------------------------------------
# 1. Clean up: remove baseline row, rows with NA p-value, non-cell-type rows
# --------------------------------------------------------------------------
df_raw <- df_raw %>%
  filter(Category != "baseL2_0", !is.na(Enrichment_p), Enrichment_p != "NA")

df_raw$Enrichment   <- as.numeric(df_raw$Enrichment)
df_raw$Enrichment_std_error <- as.numeric(df_raw$Enrichment_std_error)
df_raw$Enrichment_p <- as.numeric(df_raw$Enrichment_p)
df_raw$Prop._h2     <- as.numeric(df_raw$Prop._h2)
df_raw$Prop._SNPs   <- as.numeric(df_raw$Prop._SNPs)

n_tests <- nrow(df_raw)
bonf_threshold <- 0.05 / n_tests

df_raw <- df_raw %>%
  mutate(
    P_bonf = pmin(Enrichment_p * n_tests, 1),
    Sig    = case_when(
      P_bonf < 0.001 ~ "***",
      P_bonf < 0.01  ~ "**",
      P_bonf < 0.05  ~ "*",
      TRUE           ~ ""
    ),
    CI_lo = Enrichment - 1.96 * Enrichment_std_error,
    CI_hi = Enrichment + 1.96 * Enrichment_std_error,
    # Clean up category labels
    Label = gsub("L2_[01]$", "", Category),
    Label = gsub("\\.bed$", "", Label),
    Label = gsub("_", " ", Label),
    # Flag significant hits
    Significant = P_bonf < 0.05
  )

# --------------------------------------------------------------------------
# 2. Split into GenoSkylinePlus cell types vs baseline functional annotations
# --------------------------------------------------------------------------
df_celltypes <- df_raw %>%
  filter(grepl("L2_1$", Category)) %>%   # GenoSkylinePlus Tier3
  arrange(desc(Enrichment))

df_baseline <- df_raw %>%
  filter(grepl("L2_0$", Category)) %>%   # Baseline functional annotations
  arrange(desc(Enrichment))

# --------------------------------------------------------------------------
# 3. TOP CELL TYPES BAR PLOT — top 20 by enrichment
# --------------------------------------------------------------------------
top20 <- df_celltypes %>%
  slice_max(Enrichment, n=20) %>%
  mutate(Label = factor(Label, levels=rev(Label)))   # descending on coord_flip

p_bar <- ggplot(top20, aes(x=Label, y=Enrichment, fill=Significant)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50", linewidth=0.6) +
  geom_col(alpha=0.88, width=0.72) +
  geom_errorbar(aes(ymin=pmax(CI_lo,0), ymax=CI_hi),
                width=0.35, linewidth=0.55, colour="grey30") +
  geom_text(aes(y=CI_hi + 0.4, label=Sig), size=3.5, fontface="bold", colour="grey20") +
  coord_flip() +
  scale_fill_manual(values=c("TRUE"="#C0392B", "FALSE"="#7F8C8D"),
                    labels=c("TRUE"="Bonferroni p<0.05", "FALSE"="Not significant"),
                    name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0.02, 0.12))) +
  labs(
    title    = "Partitioned Heritability Enrichment — Male Pattern Baldness",
    subtitle = paste0("Top 20 GenoSkylinePlus Tier 3 cell types  |  ",
                      "Bonferroni threshold: p < ", formatC(bonf_threshold, format="e", digits=2)),
    x = NULL,
    y = "Fold Enrichment (Prop. h² / Prop. SNPs)",
    caption = "Significance: Bonferroni-corrected across all annotation categories\n(* p<0.05  ** p<0.01  *** p<0.001)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title    = element_text(face="bold", size=14),
    plot.subtitle = element_text(size=10, colour="grey40"),
    axis.text.y   = element_text(size=10),
    legend.position = "bottom",
    panel.grid.major.x = element_line(colour="grey92")
  )

pdf(file.path(outdir, "Baldness_celltypes_barplot.pdf"), width=12, height=9)
print(p_bar)
dev.off()
cat("Saved: Baldness_celltypes_barplot.pdf\n")

# --------------------------------------------------------------------------
# 4. ALL CELL TYPES — full dotplot sorted by enrichment
# --------------------------------------------------------------------------
df_celltypes_plot <- df_celltypes %>%
  mutate(Label = factor(Label, levels=rev(Label)))

p_dot <- ggplot(df_celltypes_plot, aes(x=Enrichment, y=Label, colour=Significant)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey50", linewidth=0.5) +
  geom_errorbarh(aes(xmin=pmax(CI_lo,0), xmax=CI_hi),
                 height=0.4, linewidth=0.4, colour="grey60") +
  geom_point(size=2.5) +
  geom_text(aes(x=CI_hi + 0.3, label=Sig), size=2.8, fontface="bold", colour="grey20") +
  scale_colour_manual(values=c("TRUE"="#C0392B", "FALSE"="#95A5A6"),
                      labels=c("TRUE"="Bonferroni p<0.05", "FALSE"="Not significant"),
                      name=NULL) +
  labs(
    title    = "Partitioned Heritability Enrichment — All Cell Types",
    subtitle = "Male Pattern Baldness | GenoSkylinePlus Tier 3",
    x = "Fold Enrichment (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size=9) +
  theme(
    plot.title    = element_text(face="bold", size=12),
    axis.text.y   = element_text(size=7),
    legend.position = "bottom",
    panel.grid.major.x = element_line(colour="grey93")
  )

pdf(file.path(outdir, "Baldness_celltypes_dotplot.pdf"), width=10, height=16)
print(p_dot)
dev.off()
cat("Saved: Baldness_celltypes_dotplot.pdf\n")

# --------------------------------------------------------------------------
# 5. BASELINE FUNCTIONAL ANNOTATIONS BAR PLOT
# --------------------------------------------------------------------------
df_base_sig <- df_baseline %>%
  filter(Significant) %>%
  mutate(Label = factor(Label, levels=rev(Label)))

p_base <- ggplot(df_base_sig, aes(x=Label, y=Enrichment)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50", linewidth=0.6) +
  geom_col(fill="#2980B9", alpha=0.85, width=0.7) +
  geom_errorbar(aes(ymin=pmax(CI_lo,0), ymax=CI_hi),
                width=0.35, linewidth=0.55, colour="grey30") +
  geom_text(aes(y=CI_hi + 0.1, label=Sig), size=3.5, fontface="bold") +
  coord_flip() +
  scale_y_continuous(expand=expansion(mult=c(0.02, 0.15))) +
  labs(
    title    = "Enrichment in Baseline Functional Annotations — Male Pattern Baldness",
    subtitle = "Bonferroni-significant categories only",
    x = NULL,
    y = "Fold Enrichment (Prop. h² / Prop. SNPs)"
  ) +
  theme_classic(base_size=11) +
  theme(
    plot.title = element_text(face="bold", size=13),
    axis.text.y = element_text(size=9),
    panel.grid.major.x = element_line(colour="grey92")
  )

pdf(file.path(outdir, "Baldness_baseline_barplot.pdf"), width=11, height=7)
print(p_base)
dev.off()
cat("Saved: Baldness_baseline_barplot.pdf\n")

# --------------------------------------------------------------------------
# 6. Print summary table of significant cell types
# --------------------------------------------------------------------------
cat("\n=== Significant Cell Types (Bonferroni p < 0.05) ===\n")
df_celltypes %>%
  filter(Significant) %>%
  select(Label, Prop._SNPs, Prop._h2, Enrichment, Enrichment_std_error, Enrichment_p, P_bonf, Sig) %>%
  arrange(Enrichment_p) %>%
  print(n=Inf)

cat("\n=== Significant Baseline Annotations (Bonferroni p < 0.05) ===\n")
df_baseline %>%
  filter(Significant) %>%
  select(Label, Prop._SNPs, Prop._h2, Enrichment, Enrichment_std_error, Enrichment_p, P_bonf, Sig) %>%
  arrange(Enrichment_p) %>%
  print(n=Inf)

cat("\nAll plots saved to:", outdir, "\n")
cat("Bonferroni threshold used: p <", bonf_threshold, "\n")
