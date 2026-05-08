#!/usr/bin/env bash
# =============================================================================
# 03_partitioned_heritability_baldness.sh
# Partitioned SNP-Heritability Enrichment for Male Pattern Baldness
#
# Method: Stratified LD Score Regression (S-LDSC; Finucane et al. 2015)
#
# Input:  Baldness.sumstats.gz (munged HapMap3 summary statistics)
# Output: results/enrichment/Baldness_enrichment.results
#         results/enrichment/Baldness_enrichment.log
#
# Annotations: GenoSkylinePlus Tier 3 (67 cell/tissue types, Lu et al. 2016)
#              jointly modeled with Baseline v1.1 functional annotations (52 categories)
#
# Enrichment = Prop. h2 / Prop. SNPs
# Bonferroni correction: p < 0.05 / 120 = 4.17e-4
#
# Author: Jenny Brant
# Course: BMI 620 Final Project
# =============================================================================
set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE=~/Desktop/stat620
LDSC="${BASE}/ldsc"
INPUTS="${BASE}/ldsc_inputs"
ENR="${INPUTS}/for_enrichment"

# Annotation file prefixes (chromosomes 1-22 appended automatically by ldsc.py)
BASELINE="${ENR}/Baseline/baseline."
TIER3="${ENR}/GenoSkylinePlus/GSplus_Tier3_1KGphase3."
WEIGHTS="${ENR}/weights/weights.hm3_noMHC."
FRQ="${ENR}/genotype/1000G.EUR.QC."

# Input summary statistics
SUMSTATS="${BASE}/Baldness.sumstats.gz"

# Output directory
ENR_OUT="${BASE}/results/enrichment"
mkdir -p "${ENR_OUT}"

# Python from ldsc conda environment
PYTHON=$(conda run -n ldsc which python)
LDSC_PY="${LDSC}/ldsc.py"

# =============================================================================
# Run S-LDSC partitioned heritability
#
# Flags:
#   --h2              munged summary statistics
#   --ref-ld-chr      comma-separated annotation prefixes (Baseline + Tier3)
#   --w-ld-chr        regression weights (HM3 SNPs, no MHC)
#   --overlap-annot   required when annotations overlap (which Baseline does)
#   --frqfile-chr     allele frequencies for filtering to common SNPs (MAF > 5%)
#   --out             output prefix
# =============================================================================
echo "================================================================"
echo " Partitioned Heritability Enrichment: Male Pattern Baldness"
echo " Annotations: Baseline v1.1 + GenoSkylinePlus Tier 3"
echo " Bonferroni threshold: p < 4.17e-4 (0.05 / 120 categories)"
echo "================================================================"

${PYTHON} "${LDSC_PY}" \
    --h2          "${SUMSTATS}" \
    --ref-ld-chr  "${BASELINE},${TIER3}" \
    --w-ld-chr    "${WEIGHTS}" \
    --overlap-annot \
    --frqfile-chr "${FRQ}" \
    --out         "${ENR_OUT}/Baldness_enrichment" \
    2>&1 | tee "${ENR_OUT}/Baldness_enrichment.log"

echo ""
echo "================================================================"
echo " Output files:"
echo "   ${ENR_OUT}/Baldness_enrichment.results"
echo "   ${ENR_OUT}/Baldness_enrichment.log"
echo ""
echo " Next step: run Baldness_EnrichmentPlot.R to generate figures"
echo "================================================================"
