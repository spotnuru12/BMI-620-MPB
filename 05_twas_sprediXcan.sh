#!/bin/bash
# ==============================================================
# 05_twas_sprediXcan.sh
# TWAS using S-PrediXcan with GTEx v7 Elastic Net skin models
#
# Input:  results/Baldness.sumstats.gz
# Models: GTEx v7 Elastic Net Skin (Sun Exposed & Not Sun Exposed)
# Output: results/Baldness_TWAS_skin_sun_exposed.csv
#         results/Baldness_TWAS_skin_not_sun_exposed.csv
#         results/Baldness_TWAS_summary.txt
#
# Note: GTEx v7 models use rs IDs compatible with munged sumstats.
#       GTEx v8 MASHR models use hg38 variant IDs and require
#       harmonization, so we use v7 for simplicity.
# ==============================================================
set -e

METAXCAN=~/MetaXcan/software
MODELS=~/bmi620_final/models
OUT=~/bmi620_final/results
GWAS=$OUT/Baldness.sumstats.gz

# Clean munged sumstats: remove rows with missing data
# (some HM3 SNPs had no match in the original GWAS)
echo "=== Cleaning munged sumstats ==="
zcat $GWAS | awk 'NR==1 || NF==5' | gzip > $OUT/Baldness_clean.sumstats.gz
GWAS_CLEAN=$OUT/Baldness_clean.sumstats.gz

# 1. Skin - Sun Exposed (Lower leg)
echo "=== Running S-PrediXcan: Skin Sun Exposed ==="
python $METAXCAN/SPrediXcan.py \
    --model_db_path $MODELS/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db \
    --covariance $MODELS/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_eur_covariances.txt.gz \
    --gwas_file $GWAS_CLEAN \
    --snp_column SNP \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --zscore_column Z \
    --output_file $OUT/Baldness_TWAS_skin_sun_exposed.csv

# 2. Skin - Not Sun Exposed (Suprapubic)
echo "=== Running S-PrediXcan: Skin Not Sun Exposed ==="
python $METAXCAN/SPrediXcan.py \
    --model_db_path $MODELS/gtex_v7_Skin_Not_Sun_Exposed_Suprapubic_imputed_europeans_tw_0.5_signif.db \
    --covariance $MODELS/gtex_v7_Skin_Not_Sun_Exposed_Suprapubic_imputed_eur_covariances.txt.gz \
    --gwas_file $GWAS_CLEAN \
    --snp_column SNP \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --zscore_column Z \
    --output_file $OUT/Baldness_TWAS_skin_not_sun_exposed.csv

# Clean up intermediate file
rm -f $OUT/Baldness_clean.sumstats.gz

# Summary: top 10 genes by p-value from each tissue
{
    echo "============================================="
    echo "TWAS Results: S-PrediXcan with GTEx v7"
    echo "Elastic Net Skin Tissue Models"
    echo "============================================="
    echo ""
    echo "--- Skin Sun Exposed (Lower Leg) - Top 10 by p-value ---"
    echo "gene_name | zscore | pvalue | n_snps_used"
    tail -n +2 $OUT/Baldness_TWAS_skin_sun_exposed.csv | \
        sort -t',' -k5,5g | head -10 | \
        awk -F',' '{printf "%s | %.2f | %s | %s\n", $2, $3, $5, $10}'
    echo ""
    echo "--- Skin Not Sun Exposed (Suprapubic) - Top 10 by p-value ---"
    echo "gene_name | zscore | pvalue | n_snps_used"
    tail -n +2 $OUT/Baldness_TWAS_skin_not_sun_exposed.csv | \
        sort -t',' -k5,5g | head -10 | \
        awk -F',' '{printf "%s | %.2f | %s | %s\n", $2, $3, $5, $10}'
    echo ""
    echo "Total genes tested (sun exposed): $(tail -n +2 $OUT/Baldness_TWAS_skin_sun_exposed.csv | grep -v 'NA,NA,NA,NA' | wc -l)"
    echo "Total genes tested (not sun exposed): $(tail -n +2 $OUT/Baldness_TWAS_skin_not_sun_exposed.csv | grep -v 'NA,NA,NA,NA' | wc -l)"
} > $OUT/Baldness_TWAS_summary.txt

echo ""
echo "=== Top TWAS Genes ==="
cat $OUT/Baldness_TWAS_summary.txt
echo ""
echo "Full results: $OUT/Baldness_TWAS_skin_sun_exposed.csv"
echo "              $OUT/Baldness_TWAS_skin_not_sun_exposed.csv"
echo "Summary:      $OUT/Baldness_TWAS_summary.txt"
