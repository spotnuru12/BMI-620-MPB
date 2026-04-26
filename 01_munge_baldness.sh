#!/bin/bash
# ==============================================================
# 01_munge_baldness.sh
# Munge Baldness GWAS summary statistics for LDSC
# 
# Input:  Baldness.tab (BOLT-LMM output, ~19.1M SNPs)
# Output: results/Baldness.sumstats.gz (munged, HM3 SNPs only)
#
# Notes:
# - Filter to autosomes (chr 1-22) since LDSC expects autosomes
# - Remove rows where P_BOLT_LMM is not numeric (extreme values)
# - Use P_BOLT_LMM (non-infinitesimal) as primary p-value
# - ALLELE1 = effect allele, ALLELE0 = reference allele
# - N = 205,327 (Hagenaars et al. 2017, PMID: 30573740)
# ==============================================================
set -e

LDSC=~/ldsc
HM3=~/ldsc/ldsc_inputs/w_hm3.snplist
WORKDIR=~/bmi620_final
OUT=$WORKDIR/results
mkdir -p $OUT

# Step 1: Filter to autosomes and ensure P is numeric
echo "Filtering to autosomes and cleaning P values..."
awk -F'\t' 'NR==1 {print; next} 
     $2+0 >= 1 && $2+0 <= 22 && $12+0 > 0 && $12+0 <= 1 {print}' \
     $WORKDIR/Baldness.tab > $WORKDIR/Baldness_autosomes.tab

echo "Original SNP count:"
wc -l $WORKDIR/Baldness.tab
echo "Autosome SNP count:"
wc -l $WORKDIR/Baldness_autosomes.tab

# Step 2: Munge
echo "Running munge_sumstats..."
python $LDSC/munge_sumstats.py \
    --sumstats $WORKDIR/Baldness_autosomes.tab \
    --N 205327 \
    --signed-sumstats BETA,0 \
    --snp SNP --a1 ALLELE1 --a2 ALLELE0 --p P_BOLT_LMM \
    --merge-alleles $HM3 --chunksize 500000 \
    --out $OUT/Baldness

echo "Munging complete. Output: $OUT/Baldness.sumstats.gz"
