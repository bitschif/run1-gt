#!/bin/bash
#===============================================================================
# STEP 06: Variant Calling - FreeBayes
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/conf/pipeline.config.sh"
source "${ROOT_DIR}/lib/helper_functions.sh"

CALLER="freebayes"

# Input
source "${PREPROC_DIR}/bam_path.sh"

OUT_DIR="${VARIANT_DIR}/${CALLER}"

#-------------------------------------------------------------------------------
# 1. Run FreeBayes
#-------------------------------------------------------------------------------

RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf"

freebayes \
    -f "${REF_FASTA}" \
    -b "${FINAL_BAM}" \
    --min-alternate-count "${FB_MIN_ALT_COUNT}" \
    --min-alternate-fraction "${FB_MIN_ALT_FRACTION}" \
    --min-mapping-quality "${MIN_MAPPING_QUALITY}" \
    --min-base-quality "${MIN_BASE_QUALITY}" \
    --genotype-qualities \
    > "${RAW_VCF}" \
    2> "${LOG_DIR}/${CALLER}.log"


# Compress and index
bgzip -f "${RAW_VCF}"
tabix -f -p vcf "${RAW_VCF}.gz"

#-------------------------------------------------------------------------------
# 2. Filter and normalize
#-------------------------------------------------------------------------------

FILTERED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz"

# Quality filter
bcftools filter \
    -i 'QUAL>30 && INFO/DP>10' \
    "${RAW_VCF}.gz" | \
bcftools norm \
    -f "${REF_FASTA}" \
    -m -both \
    -Oz -o "${FILTERED_VCF}"

tabix -f -p vcf "${FILTERED_VCF}"

#-------------------------------------------------------------------------------
# 3. Extract high-quality variants
#-------------------------------------------------------------------------------

PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"

# Add PASS tag to high-quality variants
bcftools filter \
    -i 'QUAL>30 && INFO/DP>10' \
    -s LowQual \
    "${FILTERED_VCF}" | \
bcftools view -f "PASS,." -Oz -o "${PASS_VCF}"

tabix -f -p vcf "${PASS_VCF}"

# Split by type
bcftools view -v snps "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

#-------------------------------------------------------------------------------
# 4. Normalize for benchmarking
#-------------------------------------------------------------------------------

NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
"${ROOT_DIR}/lib/normalize_vcf.sh" "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
if [[ ! -f "${TRUTH_NORM}" ]]; then
    ensure_dir "$(dirname "${TRUTH_NORM}")"
    "${ROOT_DIR}/lib/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

#-------------------------------------------------------------------------------
# 5. Stats
#-------------------------------------------------------------------------------
bcftools stats "${PASS_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${PASS_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${PASS_VCF}" | wc -l)


