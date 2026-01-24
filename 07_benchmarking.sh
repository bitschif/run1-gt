#!/bin/bash
#===============================================================================
# STEP 07: Benchmarking
# Compare variant calls against truth set
# Tools: hap.py (Illumina)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 07: Benchmarking ====="
start_timer

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")

#-------------------------------------------------------------------------------
# Main benchmarking summary
#-------------------------------------------------------------------------------
SUMMARY="${BENCH_DIR}/benchmark_summary.tsv"
echo -e "Caller\tVariantType\tTP\tFP\tFN\tPrecision\tRecall\tF1\tROC_AUC" > "${SUMMARY}"

for caller in "${CALLERS[@]}"; do
    METRICS_TSV="${BENCH_DIR}/${caller}/${PREFIX}_${caller}_metrics.tsv"
    if [[ -f "${METRICS_TSV}" ]]; then
        update_benchmark_summary "${caller}" "${METRICS_TSV}"
    else
        log_warn "  Metrics not found for ${caller}. Run variant calling script first."
    fi
done

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
log_info "===== Benchmarking Summary ====="
if command -v column &>/dev/null; then
    column -t -s$'\t' "${SUMMARY}"
else
    cat "${SUMMARY}"
fi

end_timer "07_benchmarking"
log_info "===== Benchmarking Complete ====="
