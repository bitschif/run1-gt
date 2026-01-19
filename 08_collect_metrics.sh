#!/bin/bash
#===============================================================================
# STEP 08: Collect All Metrics
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 08: Collect Metrics ====="
start_timer

ensure_dir "${METRICS_DIR}"

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")

#-------------------------------------------------------------------------------
# 1. Copy benchmark summary
#-------------------------------------------------------------------------------
cp "${BENCH_DIR}/benchmark_summary.tsv" "${METRICS_DIR}/"

#-------------------------------------------------------------------------------
# 2. Collect variant statistics
#-------------------------------------------------------------------------------
log_info "Collecting variant statistics..."

STATS_FILE="${METRICS_DIR}/variant_statistics.csv"
echo "caller,total,snps,indels,titv_ratio" > "${STATS_FILE}"

for caller in "${CALLERS[@]}"; do
    vcf="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_pass.vcf.gz"
    stats="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_stats.txt"
    
    if [[ -f "${vcf}" ]]; then
        total=$(bcftools view -H "${vcf}" | wc -l)
        snps=$(bcftools view -H -v snps "${vcf}" | wc -l)
        indels=$(bcftools view -H -v indels "${vcf}" | wc -l)
        titv=$(grep "TSTV" "${stats}" 2>/dev/null | head -1 | cut -f5 || echo "NA")
        echo "${caller},${total},${snps},${indels},${titv}" >> "${STATS_FILE}"
    fi
done

#-------------------------------------------------------------------------------
# 3. Collect runtime metrics
#-------------------------------------------------------------------------------
cp "${LOG_DIR}/runtime.csv" "${METRICS_DIR}/"

#-------------------------------------------------------------------------------
# 4. Generate summary report
#-------------------------------------------------------------------------------
log_info "Generating summary report..."

REPORT="${METRICS_DIR}/summary_report.txt"

cat > "${REPORT}" << EOF
================================================================================
VARIANT CALLING BENCHMARKING REPORT
================================================================================
Generated: $(date)
Sample: ${SAMPLE_NAME}
Chromosome: ${CHR_TO_USE}
Coverage: ${COVERAGE}x

================================================================================
TRUTH SET
================================================================================
EOF

cat "${SIM_DIR}/simulation_info.txt" >> "${REPORT}"

cat >> "${REPORT}" << EOF

================================================================================
VARIANT STATISTICS
================================================================================
EOF

column -t -s',' "${STATS_FILE}" >> "${REPORT}"

cat >> "${REPORT}" << EOF

================================================================================
BENCHMARKING RESULTS
================================================================================
EOF

column -t -s$'\t' "${BENCH_DIR}/benchmark_summary.tsv" >> "${REPORT}"

cat >> "${REPORT}" << EOF

================================================================================
RUNTIME
================================================================================
EOF

column -t -s',' "${LOG_DIR}/runtime.csv" >> "${REPORT}"

cat "${REPORT}"

end_timer "08_collect_metrics"
log_info "===== Metrics Collection Complete ====="
log_info "Report:  ${REPORT}"
