#!/bin/bash
#===============================================================================
# STEP 00: Setup Environment
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 00: Setup Environment ====="

# Create directories
for dir in "${DATA_DIR}" "${REF_DIR}" "${SIM_DIR}" "${RESULTS_DIR}" \
           "${LOG_DIR}" "${PREPROC_DIR}" "${VARIANT_DIR}" "${BENCH_DIR}" \
           "${FIGURE_DIR}" "${METRICS_DIR}"; do
    ensure_dir "$dir"
done

for caller in gatk deepvariant strelka2 freebayes; do
    ensure_dir "${VARIANT_DIR}/${caller}"
    ensure_dir "${BENCH_DIR}/${caller}"
done

# Runtime log header
echo "step,duration_seconds" > "${LOG_DIR}/runtime.csv"

# Check required tools
log_info "Checking required tools..."
TOOLS_OK=true

for tool in bwa samtools bcftools gatk fastp fastqc bgzip tabix art_illumina freebayes simutator; do
    if check_tool "$tool"; then
        log_info "  ✓ $tool"
    else
        log_warn "  ✗ $tool (missing)"
        TOOLS_OK=false
    fi
done

# Check Docker
if check_tool docker; then
    log_info "  ✓ docker (for DeepVariant, Strelka2)"
else
    log_warn "  ✗ docker (DeepVariant, Strelka2 sẽ không chạy được)"
fi

# Download reference if needed
if [[ ! -f "${REF_FASTA}" ]]; then
    log_info "Downloading ${CHR_TO_USE} from UCSC..."
    cd "${REF_DIR}"
    wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${CHR_TO_USE}. fa.gz"
    gunzip "${CHR_TO_USE}.fa.gz"
    
    log_info "Indexing reference..."
    samtools faidx "${REF_FASTA}"
    bwa index "${REF_FASTA}"
    gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}"
fi

log_info "===== Setup Complete ====="