#!/bin/bash
# ==============================================================
# 02_ldsc_h2.sh
# LDSC heritability and intercept for Baldness GWAS
#
# Input:  results/Baldness.sumstats.gz
# Output: results/Baldness_h2.log (full LDSC log)
#         results/Baldness_h2_summary.txt (key results)
#
# Note: The @ in the home directory path confuses LDSC's
# chromosome file parsing, so we symlink LD scores to /tmp/ld_chr
# ==============================================================
set -e

LDSC=~/ldsc
OUT=~/bmi620_final/results

# Create symlink to avoid path issues with @ in home directory
ln -sf ~/ldsc/ldsc_inputs/for_h2/eur_w_ld_chr /tmp/ld_chr
LDCHR=/tmp/ld_chr/

python $LDSC/ldsc.py \
    --h2 $OUT/Baldness.sumstats.gz \
    --ref-ld-chr $LDCHR \
    --w-ld-chr $LDCHR \
    --out $OUT/Baldness_h2

# Save key results
{
    echo "========================================="
    echo "LDSC Heritability Results: Baldness (MPB)"
    echo "========================================="
    echo ""
    grep "Total Observed" $OUT/Baldness_h2.log
    grep "Lambda GC" $OUT/Baldness_h2.log
    grep "Mean Chi" $OUT/Baldness_h2.log
    grep "Intercept" $OUT/Baldness_h2.log
    grep "Ratio" $OUT/Baldness_h2.log
} > $OUT/Baldness_h2_summary.txt

echo ""
echo "=== Key Results ==="
cat $OUT/Baldness_h2_summary.txt
echo ""
echo "Full log: $OUT/Baldness_h2.log"
echo "Summary:  $OUT/Baldness_h2_summary.txt"
