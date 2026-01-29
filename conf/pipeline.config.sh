#!/bin/bash
#===============================================================================
# CONFIG: Cấu hình chung cho pipeline variant calling benchmarking
# Theo GATK Best Practices (nf-core/sarek)
#===============================================================================

#-------------------------------------------------------------------------------
# SYSTEM RESOURCES (RAM 16GB, 4 CPU)
#-------------------------------------------------------------------------------
export THREADS="${THREADS:-4}"
export MAX_MEMORY="${MAX_MEMORY:-14G}"
export JAVA_OPTS="${JAVA_OPTS:--Xmx12G -XX:ParallelGCThreads=2}"

#-------------------------------------------------------------------------------
# REFERENCE GENOME
#-------------------------------------------------------------------------------
export GENOME_VERSION="${GENOME_VERSION:-hg38}"
export CHR_TO_USE="${CHR_TO_USE:-chr22}"

#-------------------------------------------------------------------------------
# DIRECTORY PATHS
#-------------------------------------------------------------------------------
export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
export DATA_DIR="${DATA_DIR:-${PROJECT_DIR}/data}"
export RESULTS_DIR="${OUTDIR:-${PROJECT_DIR}/results}"
export LOG_DIR="${LOG_DIR:-${PROJECT_DIR}/logs}"

export REF_FASTA="${FASTA:-}"
if [[ -n "${REF_FASTA}" ]]; then
    export REF_DIR="${REF_DIR:-$(cd "$(dirname "${REF_FASTA}")" && pwd)}"
else
    export REF_DIR="${REF_DIR:-${DATA_DIR}/reference}"
    export REF_FASTA="${REF_DIR}/${CHR_TO_USE}.fa"
fi

export SIM_DIR="${SIM_DIR:-${DATA_DIR}/simulated}"

export PREPROC_DIR="${RESULTS_DIR}/preprocessing"
export VARIANT_DIR="${RESULTS_DIR}/variants"
export BENCH_DIR="${RESULTS_DIR}/benchmarks"
export FIGURE_DIR="${RESULTS_DIR}/plots"
export METRICS_DIR="${RESULTS_DIR}/final_metrics"

#-------------------------------------------------------------------------------
# REFERENCE GENOME PATHS
#-------------------------------------------------------------------------------
export REF_DICT="${DICT:-${REF_DIR}/${CHR_TO_USE}.dict}"
export REF_FAI="${FASTA_FAI:-${REF_FASTA}.fai}"
export BWA_INDEX="${BWA_INDEX:-}"
export BWA_REFERENCE="${BWA_INDEX:-${REF_FASTA}}"

# Known sites for BQSR (GATK Bundle)
export DBSNP="${KNOWN_SNPS:-${REF_DIR}/dbsnp_146.hg38.${CHR_TO_USE}.vcf.gz}"
export KNOWN_INDELS="${KNOWN_INDELS:-${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.${CHR_TO_USE}.vcf.gz}"
export KNOWN_SNPS="${KNOWN_SNPS:-${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.${CHR_TO_USE}.vcf.gz}"

#-------------------------------------------------------------------------------
# SIMULATION PARAMETERS (simutator + ART)
#-------------------------------------------------------------------------------
# SNP_DIST=7000 -> ~7,000 SNPs on chr22 (~50Mb)
export SNP_DIST=7000
export DEL_DIST=2000
export DEL_LEN=3
export INS_DIST=2000
export INS_LEN=2

# ART Illumina parameters
export COVERAGE=60
export READ_LENGTH=150
export FRAGMENT_MEAN=350
export FRAGMENT_SD=50
export ART_PLATFORM="HS25"
export SEED=42

#-------------------------------------------------------------------------------
# QUALITY THRESHOLDS
#-------------------------------------------------------------------------------
export MIN_BASE_QUALITY=20
export MIN_MAPPING_QUALITY=20
export MIN_READ_LENGTH=50

#-------------------------------------------------------------------------------
# SAMPLE INFO
#-------------------------------------------------------------------------------
export SAMPLE_NAME="${SAMPLE_NAME:-SIMULATED_SAMPLE}"
export PREFIX="${SAMPLE_NAME}_${CHR_TO_USE}"
export READ_GROUP_ID="${READGROUP_ID:-${SAMPLE_NAME}}"
export SAMPLE_LANE="${SAMPLE_LANE:-${LANE:-1}}"
export PLATFORM="${PLATFORM:-ILLUMINA}"
export LIBRARY="${LIBRARY:-lib1}"
export CENTER="${CENTER:-}"
export INSTRUMENT="${INSTRUMENT:-}"
export READ_GROUP="${READ_GROUP:-@RG\\tID:${READ_GROUP_ID}\\tSM:${SAMPLE_NAME}\\tPL:${PLATFORM}\\tLB:${LIBRARY}\\tPU:${SAMPLE_LANE}}"

#-------------------------------------------------------------------------------
# TOOL VERSIONS (align with nf-core/sarek latest)
#-------------------------------------------------------------------------------
export DEEPVARIANT_VERSION="1.7.0"
export GATK_VERSION="4.6.0.0"
export FREEBAYES_VERSION="1.3.6"
export STRELKA2_VERSION="2.9.10"

#-------------------------------------------------------------------------------
# DOCKER IMAGES
#-------------------------------------------------------------------------------
export DEEPVARIANT_IMAGE="google/deepvariant:${DEEPVARIANT_VERSION}"
export STRELKA2_IMAGE="quay.io/biocontainers/strelka:${STRELKA2_VERSION}--h9ee0642_1"

#-------------------------------------------------------------------------------
# VARIANT CALLER PARAMETERS
#-------------------------------------------------------------------------------
# GATK HaplotypeCaller
export GATK_STAND_CALL_CONF=30

# FreeBayes
export FB_MIN_ALT_COUNT=3
export FB_MIN_ALT_FRACTION=0.2

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------
export TRUTH_VCF="${TRUTH_VCF:-${SIM_DIR}/${PREFIX}_truth.vcf.gz}"
export HIGH_CONF_BED="${INTERVALS:-${TRUTH_BED:-${SIM_DIR}/callable_regions.bed}}"

echo "[CONFIG] Loaded successfully"
echo "[CONFIG] Project: ${PROJECT_DIR}"
echo "[CONFIG] Chromosome: ${CHR_TO_USE}, Coverage: ${COVERAGE}x"
