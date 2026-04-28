#!/bin/bash
# ==============================================================
# 03_ldsc_rg.sh
# Pairwise genetic correlations between Baldness and 5 traits
#
# Output: results/Baldness_rg.log (full LDSC log)
#         results/Baldness_rg_summary.txt (key results)
# ==============================================================
set -e

LDSC=~/ldsc
ln -sf ~/ldsc/ldsc_inputs/for_h2/eur_w_ld_chr /tmp/ld_chr
LDCHR=/tmp/ld_chr/
OUT=~/bmi620_final/results

python $LDSC/ldsc.py \
    --rg $OUT/Baldness.sumstats.gz,$OUT/Chronotype.sumstats.gz,$OUT/ADHD.sumstats.gz,$OUT/Cannabis.sumstats.gz,$OUT/Insomnia.sumstats.gz,$OUT/Income.sumstats.gz \
    --ref-ld-chr $LDCHR \
    --w-ld-chr $LDCHR \
    --out $OUT/Baldness_rg

# Save key results
{
    echo "============================================="
    echo "Genetic Correlations with Baldness (MPB)"
    echo "============================================="
    echo ""
    echo "p1: Baldness"
    echo ""
    # Extract the results table
    sed -n '/^p1/,/^$/p' $OUT/Baldness_rg.log
} > $OUT/Baldness_rg_summary.txt

echo ""
echo "=== Genetic Correlations with Baldness ==="
cat $OUT/Baldness_rg_summary.txt
echo ""
echo "Full log: $OUT/Baldness_rg.log"
echo "Summary:  $OUT/Baldness_rg_summary.txt"
