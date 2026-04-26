#!/bin/bash
# ==============================================================
# 01b_munge_other_gwas.sh
# Munge 5 GWAS summary statistics for genetic correlations,
# SUPERGNOVA, and other downstream analyses.
#
# Traits:
#   1. Chronotype (quantitative, N=449,734)
#   2. ADHD (binary, Nca=38,691 Nco=186,843)
#   3. Cannabis Dependence (binary, per-SNP N in file)
#   4. Insomnia (binary, per-SNP N in NMISS column)
#   5. Income (quantitative, N=286,301, Hill et al. 2019)
# ==============================================================
set -e

LDSC=~/ldsc
HM3=~/ldsc/ldsc_inputs/w_hm3.snplist
WORKDIR=~/bmi620_final
OUT=$WORKDIR/results
mkdir -p $OUT

# 1. CHRONOTYPE - quantitative, BOLT-LMM
echo "=== Munging Chronotype ==="
python $LDSC/munge_sumstats.py \
    --sumstats $WORKDIR/Chronotype.txt.gz \
    --N 449734 \
    --signed-sumstats BETA,0 \
    --snp SNP --a1 ALLELE1 --a2 ALLELE0 --p P_BOLT_LMM \
    --merge-alleles $HM3 --chunksize 500000 \
    --out $OUT/Chronotype

# 2. ADHD - binary, case-control
echo "=== Munging ADHD ==="
python $LDSC/munge_sumstats.py \
    --sumstats $WORKDIR/Attention-deficit_hyperactivity_disorder_\(ADHD\).gz \
    --N-cas 38691 --N-con 186843 \
    --signed-sumstats OR,1 \
    --snp SNP --a1 A1 --a2 A2 --p P \
    --merge-alleles $HM3 --chunksize 500000 \
    --out $OUT/ADHD

# 3. CANNABIS DEPENDENCE - binary, per-SNP N
echo "=== Munging Cannabis Dependence ==="
python $LDSC/munge_sumstats.py \
    --sumstats $WORKDIR/Cannabis_Dependence \
    --signed-sumstats Z,0 \
    --snp SNP --a1 A1 --a2 A2 --p P --N-col N \
    --merge-alleles $HM3 --chunksize 500000 \
    --out $OUT/Cannabis

# 4. INSOMNIA - binary, per-SNP N
# Preprocessing: Insomnia file has two columns LDSC recognizes as
# SNP identifiers (SNP in chr:pos format, RSID_UKB with rs IDs).
# We drop the original SNP column and use RSID_UKB instead.
echo "=== Munging Insomnia ==="
zcat $WORKDIR/Insomnia.txt.gz | \
    awk -F'\t' '{
        for(i=2;i<=NF;i++) printf "%s%s", $i, (i<NF?"\t":"\n")
    }' | gzip > $WORKDIR/Insomnia_clean.txt.gz

python $LDSC/munge_sumstats.py \
    --sumstats $WORKDIR/Insomnia_clean.txt.gz \
    --signed-sumstats OR,1 \
    --snp RSID_UKB --a1 A1 --a2 A2 --p P --N-col NMISS \
    --merge-alleles $HM3 --chunksize 500000 \
    --out $OUT/Insomnia

rm -f $WORKDIR/Insomnia_clean.txt.gz

# 5. INCOME - quantitative
echo "=== Munging Income ==="
python $LDSC/munge_sumstats.py \
    --sumstats $WORKDIR/Income.txt.gz \
    --N 286301 \
    --signed-sumstats Beta,0 \
    --snp SNP --a1 Effect_Allele --a2 Non_effect_Allele --p P \
    --merge-alleles $HM3 --chunksize 500000 \
    --out $OUT/Income

# Summary
echo ""
echo "=== Munging Summary ==="
for TRAIT in Chronotype ADHD Cannabis Insomnia Income; do
    if [ -f $OUT/${TRAIT}.sumstats.gz ]; then
        NSNP=$(zcat $OUT/${TRAIT}.sumstats.gz | wc -l)
        echo "$TRAIT: $((NSNP - 1)) SNPs"
    else
        echo "$TRAIT: FAILED"
    fi
done
