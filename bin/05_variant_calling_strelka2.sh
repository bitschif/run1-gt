#!/bin/bash
#===============================================================================
# STEP 05: Variant Calling - Strelka2 Germline (via Docker)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/conf/pipeline.config.sh"
source "${ROOT_DIR}/lib/helper_functions.sh"

CALLER="strelka2"

# Check Docker

# Input
source "${PREPROC_DIR}/bam_path.sh"

OUT_DIR="${VARIANT_DIR}/${CALLER}"
STRELKA_RUNDIR="${OUT_DIR}/strelka_run"
rm -rf "${STRELKA_RUNDIR}"
ensure_dir "${STRELKA_RUNDIR}"

#-------------------------------------------------------------------------------
# 1. Prepare paths for Docker
#-------------------------------------------------------------------------------
ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}" && pwd)
ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)

BAM_BASENAME=$(basename "${FINAL_BAM}")
REF_BASENAME=$(basename "${REF_FASTA}")

#-------------------------------------------------------------------------------
# 2. Configure Strelka2 via Docker
#-------------------------------------------------------------------------------

docker run \
    --rm \
    -v "${ABS_REF_DIR}:/ref:ro" \
    -v "${ABS_PREPROC_DIR}:/input:ro" \
    -v "${ABS_OUT_DIR}:/output" \
    ${STRELKA2_IMAGE} \
    configureStrelkaGermlineWorkflow.py \
    --bam "/input/${BAM_BASENAME}" \
    --referenceFasta "/ref/${REF_BASENAME}" \
    --runDir "/output/strelka_run" \
    2>&1 | tee "${LOG_DIR}/${CALLER}_config.log"


#-------------------------------------------------------------------------------
# 3. Run Strelka2 via Docker
#-------------------------------------------------------------------------------

docker run \
    --rm \
    -v "${ABS_REF_DIR}:/ref:ro" \
    -v "${ABS_PREPROC_DIR}:/input:ro" \
    -v "${ABS_OUT_DIR}:/output" \
    ${STRELKA2_IMAGE} \
    /output/strelka_run/runWorkflow.py \
    -m local \
    -j "${THREADS}" \
    2>&1 | tee "${LOG_DIR}/${CALLER}_run.log"


#-------------------------------------------------------------------------------
# 4. Process output
#-------------------------------------------------------------------------------

STRELKA_VCF="${STRELKA_RUNDIR}/results/variants/variants.vcf.gz"

if [[ !  -f "${STRELKA_VCF}" ]]; then
    log_error "Strelka2 output not found: ${STRELKA_VCF}"
    exit 1
fi

# Copy to standard location
RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf.gz"
cp "${STRELKA_VCF}" "${RAW_VCF}"
tabix -f -p vcf "${RAW_VCF}"

# Extract PASS variants
PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"
bcftools view -f PASS "${RAW_VCF}" -Oz -o "${PASS_VCF}"
tabix -f -p vcf "${PASS_VCF}"

# Split by type
bcftools view -v snps "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

#-------------------------------------------------------------------------------
# 5. Normalize for benchmarking
#-------------------------------------------------------------------------------

NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
"${ROOT_DIR}/lib/normalize_vcf.sh" "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
if [[ ! -f "${TRUTH_NORM}" ]]; then
    ensure_dir "$(dirname "${TRUTH_NORM}")"
    "${ROOT_DIR}/lib/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

#-------------------------------------------------------------------------------
# 6. Stats
#-------------------------------------------------------------------------------
bcftools stats "${PASS_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${PASS_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${PASS_VCF}" | wc -l)


