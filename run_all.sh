#!/bin/bash
#===============================================================================
# MASTER SCRIPT:  Run entire variant calling benchmarking pipeline
# Usage: ./run_all.sh [OPTIONS]
#
# Options:
#   --skip-setup         Skip environment setup
#   --skip-simulation    Skip data simulation
#   --skip-preprocessing Skip preprocessing
#   --skip-calling       Skip variant calling
#   --callers LIST       Comma-separated callers (gatk,deepvariant,strelka2,freebayes)
#   -h, --help           Show help
#===============================================================================

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR="${SCRIPT_DIR}"

# Default options
SKIP_SETUP=false
SKIP_SIMULATION=false
SKIP_PREPROCESSING=false
SKIP_CALLING=false
CALLERS=""

#-------------------------------------------------------------------------------
# Parse arguments
#-------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-setup)
            SKIP_SETUP=true
            shift
            ;;
        --skip-simulation)
            SKIP_SIMULATION=true
            shift
            ;;
        --skip-preprocessing)
            SKIP_PREPROCESSING=true
            shift
            ;;
        --skip-calling)
            SKIP_CALLING=true
            shift
            ;;
        --callers)
            CALLERS="$2"
            shift 2
            ;;
        -h|--help)
            head -20 "$0" | tail -15
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

#-------------------------------------------------------------------------------
# Print banner
#-------------------------------------------------------------------------------
echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║     VARIANT CALLING BENCHMARKING PIPELINE (Nextflow)              ║"
echo "║     GATK | DeepVariant | Strelka2 | FreeBayes                      ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""

echo "Start time: $(date)"
echo "Project directory: ${PROJECT_DIR}"
echo ""

#-------------------------------------------------------------------------------
# Run Nextflow pipeline
#-------------------------------------------------------------------------------
NF_ARGS=(run "${SCRIPT_DIR}/main.nf")

if [[ "${SKIP_SETUP}" == true ]]; then
    NF_ARGS+=(--run_setup false)
fi
if [[ "${SKIP_SIMULATION}" == true ]]; then
    NF_ARGS+=(--run_simulation false)
fi
if [[ "${SKIP_PREPROCESSING}" == true ]]; then
    NF_ARGS+=(--run_preprocessing false)
fi
if [[ "${SKIP_CALLING}" == true ]]; then
    NF_ARGS+=(--run_calling false)
fi
if [[ -n "${CALLERS}" ]]; then
    NF_ARGS+=(--callers "${CALLERS}")
fi

nextflow "${NF_ARGS[@]}"
