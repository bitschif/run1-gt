#!/bin/bash
#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $1" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        return 1
    fi
    return 0
}

check_tool() {
    if ! command -v "$1" &> /dev/null; then
        log_error "Tool not found: $1"
        return 1
    fi
    return 0
}

ensure_dir() {
    [[ ! -d "$1" ]] && mkdir -p "$1"
}

start_timer() {
    export STEP_START_TIME=$(date +%s)
}

end_timer() {
    local step_name="$1"
    local end_time=$(date +%s)
    local duration=$((end_time - STEP_START_TIME))
    log_info "${step_name} completed in ${duration} seconds"
    echo "${step_name},${duration}" >> "${LOG_DIR}/runtime.csv"
}

check_exit() {
    if [[ $? -ne 0 ]]; then
        log_error "$1 failed"
        exit 1
    fi
    log_info "$1 completed"
}

get_python_bin() {
    if command -v python3 &>/dev/null; then
        echo "python3"
        return 0
    fi
    if command -v python &>/dev/null; then
        echo "python"
        return 0
    fi
    return 1
}

normalize_vcf() {
    local input_vcf="$1"
    local output_vcf="$2"
    local reference_fasta="$3"

    check_tool bcftools || return 1

    local tmp_vcf="${output_vcf}.tmp"
    bcftools norm -f "${reference_fasta}" -m -both "${input_vcf}" -Oz -o "${tmp_vcf}"
    bcftools sort "${tmp_vcf}" -Oz -o "${output_vcf}"
    rm -f "${tmp_vcf}"
    tabix -f -p vcf "${output_vcf}"
}

run_happy() {
    local truth_vcf="$1"
    local query_vcf="$2"
    local out_prefix="$3"

    check_tool "hap.py" || return 1
    hap.py "${truth_vcf}" "${query_vcf}" \
        -f "${HIGH_CONF_BED}" \
        -r "${REF_FASTA}" \
        -o "${out_prefix}"
}

write_happy_metrics() {
    local summary_csv="$1"
    local caller="$2"
    local hap_prefix="$3"
    local output_tsv="$4"

    local python_bin
    if ! python_bin=$(get_python_bin); then
        log_error "Python not found for parsing hap.py summary"
        return 1
    fi

    "${python_bin}" - "${summary_csv}" "${caller}" "${hap_prefix}" "${output_tsv}" <<'PY'
import csv
import os
import re
import sys

summary_csv, caller, hap_prefix, output_tsv = sys.argv[1:5]

def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())

def detect_dialect(path: str) -> csv.Dialect:
    with open(path, "r", newline="") as handle:
        sample = handle.read(2048)
        handle.seek(0)
        try:
            return csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            return csv.excel

def read_table(path: str):
    dialect = detect_dialect(path)
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, dialect=dialect)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return rows, fieldnames

def auc_from_roc(path: str) -> str:
    if not path or not os.path.exists(path):
        return "NA"
    rows, fieldnames = read_table(path)
    if not fieldnames:
        return "NA"
    normalized = {norm(name): name for name in fieldnames}
    fpr_key = normalized.get("fpr")
    tpr_key = normalized.get("tpr")
    if not fpr_key or not tpr_key:
        return "NA"
    points = []
    for row in rows:
        try:
            fpr = float(row.get(fpr_key, ""))
            tpr = float(row.get(tpr_key, ""))
        except ValueError:
            continue
        points.append((fpr, tpr))
    if not points:
        return "NA"
    points.sort()
    auc = 0.0
    for (x1, y1), (x2, y2) in zip(points, points[1:]):
        auc += (x2 - x1) * (y1 + y2) / 2.0
    return f"{auc:.6f}"

def find_roc(prefix: str, variant: str) -> str:
    candidates = [
        f"{prefix}.roc.{variant}.csv",
        f"{prefix}.roc.{variant}.tsv",
        f"{prefix}.roc.{variant}.txt",
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return ""

rows, fieldnames = read_table(summary_csv)
if not fieldnames:
    sys.exit(0)

normalized = {norm(name): name for name in fieldnames}

def pick(*names):
    for name in names:
        key = norm(name)
        if key in normalized:
            return normalized[key]
    return None

type_key = pick("type")
tp_key = pick("tp")
fp_key = pick("fp")
fn_key = pick("fn")
precision_key = pick("precision", "prec")
recall_key = pick("recall")
f1_key = pick("f1", "f1score")

if not all([type_key, tp_key, fp_key, fn_key, precision_key, recall_key, f1_key]):
    sys.exit(0)

wanted = {
    "ALL": {"all", "total", "overall"},
    "SNP": {"snp"},
    "INDEL": {"indel"},
}

header = ["Caller", "VariantType", "TP", "FP", "FN", "Precision", "Recall", "F1", "ROC_AUC"]
lines = [header]

for row in rows:
    row_type = row.get(type_key, "").strip().lower()
    for out_type, labels in wanted.items():
        if row_type in labels:
            roc_key = out_type.lower()
            roc_file = find_roc(hap_prefix, roc_key)
            roc_auc = auc_from_roc(roc_file)
            lines.append([
                caller,
                out_type,
                row.get(tp_key, ""),
                row.get(fp_key, ""),
                row.get(fn_key, ""),
                row.get(precision_key, ""),
                row.get(recall_key, ""),
                row.get(f1_key, ""),
                roc_auc,
            ])
            break

with open(output_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerows(lines)
PY
}

update_benchmark_summary() {
    local caller="$1"
    local metrics_tsv="$2"
    local runtime_seconds="${3:-NA}"
    local summary_file="${BENCH_DIR}/benchmark_summary.tsv"
    local header="Caller\tVariantType\tTP\tFP\tFN\tPrecision\tRecall\tF1\tROC_AUC\tRuntimeSeconds"

    if [[ -f "${summary_file}" ]]; then
        awk -v caller="${caller}" 'NR==1 {print; next} $1!=caller {print}' "${summary_file}" > "${summary_file}.tmp"
        mv "${summary_file}.tmp" "${summary_file}"
    else
        echo -e "${header}" > "${summary_file}"
    fi

    if [[ -f "${metrics_tsv}" ]]; then
        awk -v runtime="${runtime_seconds}" 'NR>1 {print $0 "\t" runtime}' "${metrics_tsv}" >> "${summary_file}"
    fi
}
