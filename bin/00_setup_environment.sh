#!/bin/bash
#===============================================================================
# STEP 00: Setup Environment
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/conf/pipeline.config.sh"
source "${ROOT_DIR}/lib/helper_functions.sh"

log_info "===== STEP 00: Setup Environment ====="

#-------------------------------------------------------------------------------
# 1. Create directories
#-------------------------------------------------------------------------------
for dir in "${DATA_DIR}" "${REF_DIR}" "${SIM_DIR}" "${RESULTS_DIR}" \
           "${LOG_DIR}" "${PREPROC_DIR}" "${VARIANT_DIR}" "${BENCH_DIR}" \
           "${FIGURE_DIR}" "${METRICS_DIR}"; do
    ensure_dir "$dir"
done

for caller in gatk deepvariant strelka2 freebayes; do
    ensure_dir "${VARIANT_DIR}/${caller}"
    ensure_dir "${BENCH_DIR}/${caller}"
done

echo "step,duration_seconds" > "${LOG_DIR}/runtime.csv"

#-------------------------------------------------------------------------------
# 2. Check required tools
#-------------------------------------------------------------------------------
log_info "Checking required tools..."

for tool in bwa samtools bcftools gatk fastp fastqc bgzip tabix art_illumina freebayes simutator docker; do
    if command -v "$tool" &> /dev/null; then
        log_info "  ✓ $tool"
    else
        log_warn "  ✗ $tool (missing)"
    fi
done

#-------------------------------------------------------------------------------
# 3. Download reference genome
#-------------------------------------------------------------------------------
if [[ ! -f "${REF_FASTA}" ]]; then
    log_info "Downloading ${CHR_TO_USE} from UCSC..."
    cd "${REF_DIR}"
    
    UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${CHR_TO_USE}.fa.gz"
    ENSEMBL_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
    
    if wget -q --show-progress --timeout=60 -O "${CHR_TO_USE}.fa.gz" "${UCSC_URL}" 2>/dev/null || \
       wget -q --show-progress --timeout=60 -O "${CHR_TO_USE}.fa.gz" "${ENSEMBL_URL}" 2>/dev/null || \
       curl -L -f -o "${CHR_TO_USE}.fa.gz" "${UCSC_URL}" 2>/dev/null; then
        log_info "  Download successful"
    else
        log_error "Failed to download reference genome!"
        log_info "Please download manually from: ${UCSC_URL}"
        log_info "Place the file in: ${REF_DIR}/"
        exit 1
    fi
    
    # Check file size (should be > 1MB for chr22)
    FILE_SIZE=$(stat -c%s "${CHR_TO_USE}.fa.gz" 2>/dev/null || stat -f%z "${CHR_TO_USE}.fa.gz" 2>/dev/null || echo "0")
    if [[ "${FILE_SIZE}" -lt 1000000 ]]; then
        log_error "Downloaded file is too small (${FILE_SIZE} bytes). Download may have failed."
        rm -f "${CHR_TO_USE}.fa.gz"
        exit 1
    fi
    
    log_info "  Decompressing..."
    gunzip -f "${CHR_TO_USE}.fa.gz"
    
    if [[ ! -f "${REF_FASTA}" ]]; then
        log_error "Decompression failed!"
        exit 1
    fi
    
    log_info "  Indexing reference..."
    samtools faidx "${REF_FASTA}"
    bwa index "${REF_FASTA}"
    gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}" 2>/dev/null || \
        samtools dict "${REF_FASTA}" > "${REF_DICT}"
    
    log_info "  Reference preparation complete!"
else
    log_info "Reference already exists: ${REF_FASTA}"
    
    [[ ! -f "${REF_FAI}" ]] && samtools faidx "${REF_FASTA}"
    [[ ! -f "${REF_FASTA}.bwt" ]] && bwa index "${REF_FASTA}"
    [[ ! -f "${REF_DICT}" ]] && (gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}" 2>/dev/null || samtools dict "${REF_FASTA}" > "${REF_DICT}")
fi

#-------------------------------------------------------------------------------
# 4. Verify reference files
#-------------------------------------------------------------------------------
log_info "Verifying reference files..."

MISSING_FILES=false
for file in "${REF_FASTA}" "${REF_FAI}"; do
    if [[ -f "${file}" ]]; then
        log_info "  ✓ $(basename ${file})"
    else
        log_error "  ✗ $(basename ${file}) - MISSING"
        MISSING_FILES=true
    fi
done

if [[ "${MISSING_FILES}" == true ]]; then
    log_error "Some reference files are missing!"
    exit 1
fi

REF_SIZE=$(ls -lh "${REF_FASTA}" | awk '{print $5}')
REF_SEQS=$(grep -c "^>" "${REF_FASTA}")
log_info "  Reference size: ${REF_SIZE}, Sequences: ${REF_SEQS}"

log_info "===== Setup Complete ====="
