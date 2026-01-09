#!/bin/bash
set -euo pipefail

# ============================================
# Batch Worker Script
# ============================================
# This script processes all samples in a given sample list file.
# It is intended to be submitted as an independent SLURM job.
#
# Usage:
#   sbatch process_batch_worker.sh <sample_list_file>
# ============================================

# SLURM directives (can be overridden at submission time)
#SBATCH --job-name=TE_batch
#SBATCH --output=logs/TE_batch_%j.out
#SBATCH --error=logs/TE_batch_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# Validate argument
if [[ $# -lt 1 ]]; then
    echo "ERROR: Sample list file required"
    echo "Usage: $0 <sample_list_file>"
    exit 1
fi

SAMPLE_LIST="$1"

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST}"
    exit 1
fi

echo "========================================="
echo "Batch Worker - Job ${SLURM_JOB_ID:-UNKNOWN}"
echo "========================================="
echo "Sample list: ${SAMPLE_LIST}"
echo "Started: $(date)"
echo "========================================="
echo ""

# Determine script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load configuration
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
    echo "Loading configuration from config.sh..."
    source "${SCRIPT_DIR}/config.sh"
elif [[ -f "${SCRIPT_DIR}/config_template.sh" ]]; then
    echo "WARNING: config.sh not found, using config_template.sh as fallback"
    source "${SCRIPT_DIR}/config_template.sh"
else
    echo "ERROR: Neither config.sh nor config_template.sh found"
    exit 1
fi

# Source batch configuration (optional)
if [[ -f "${SCRIPT_DIR}/batch-config.sh" ]]; then
    echo "Loading batch configuration from batch-config.sh..."
    source "${SCRIPT_DIR}/batch-config.sh"
fi

# Validate required configuration
required_vars=(
    "DATA_HOME"
    "OUTPUT_BASE"
    "JET2"
    "TEProf2"
    "refDir"
)

missing_vars=()
for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        missing_vars+=("$var")
    fi
done

if [[ ${#missing_vars[@]} -gt 0 ]]; then
    echo "ERROR: Required configuration variables not set:"
    printf '  - %s\n' "${missing_vars[@]}"
    exit 1
fi

# Load modules if needed
if command -v module &>/dev/null; then
    module load singularity 2>/dev/null || true
fi

# Create logs directory
mkdir -p logs

# Source shared functions
source "${SCRIPT_DIR}/function.sh"

# Track overall status
OVERALL_SUCCESS=0
FAILED_SAMPLES=()
PROCESSED_COUNT=0

# Count total samples
TOTAL_SAMPLES=$(grep -v "^#" "${SAMPLE_LIST}" | wc -l)

echo "Processing ${TOTAL_SAMPLES} samples..."
echo ""

# Process each sample in the list
while IFS=' ' read -r SAMPLE_NAME FQ1 FQ2 REL_PATH; do
    # Skip header or comments
    if [[ "${SAMPLE_NAME}" == "#"* ]] || [[ -z "${SAMPLE_NAME}" ]]; then
        continue
    fi
    
    PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
    
    # Setup sample output directory
    dataDir="${OUTPUT_BASE}/${REL_PATH}/${SAMPLE_NAME}"
    mkdir -p "${dataDir}"
    mkdir -p "${dataDir}/log"
    mkdir -p "${dataDir}/err"
    
    # Set outputDir for JET functions
    outputDir="${dataDir}"
    
    # Export required variables for functions
    export FQ1 FQ2 SAMPLE_NAME outputDir dataDir
    
    echo "========================================="
    echo "Processing sample ${PROCESSED_COUNT}/${TOTAL_SAMPLES}: ${SAMPLE_NAME}"
    echo "Input files:"
    echo "  R1: ${FQ1}"
    echo "  R2: ${FQ2}"
    echo "Output directory: ${dataDir}"
    echo "========================================="
    
    # Run JET Step 1
    run_jet_step1
    JET_STEP1_STATUS=$?
    
    # Run JET Step 2 (only if Step 1 succeeded)
    if [[ ${JET_STEP1_STATUS} -eq 0 ]]; then
        run_jet_step2
        JET_STEP2_STATUS=$?
    else
        JET_STEP2_STATUS=1
    fi
    
    # Run TEProf2
    run_teprof2
    TEPROF2_STATUS=$?
    
    # Sample summary
    echo ""
    echo "========================================="
    echo "Processing Summary for ${SAMPLE_NAME}"
    echo "========================================="
    echo "JET Step 1 Status: $([[ ${JET_STEP1_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "JET Step 2 Status: $([[ ${JET_STEP2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "TEProf2 Status: $([[ ${TEPROF2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "========================================="
    echo ""
    
    # Track failures
    if [[ ${JET_STEP1_STATUS} -ne 0 ]] || [[ ${JET_STEP2_STATUS} -ne 0 ]] || [[ ${TEPROF2_STATUS} -ne 0 ]]; then
        echo "ERROR: One or more pipeline steps failed for ${SAMPLE_NAME}"
        FAILED_SAMPLES+=("${SAMPLE_NAME}")
        OVERALL_SUCCESS=1
    else
        echo "SUCCESS: All pipeline steps completed for ${SAMPLE_NAME}"
    fi
    
done < "${SAMPLE_LIST}"

# Job summary
echo ""
echo "========================================="
echo "Batch Worker Job Summary"
echo "========================================="
echo "Job ID: ${SLURM_JOB_ID:-UNKNOWN}"
echo "Sample list: ${SAMPLE_LIST}"
echo "Total samples processed: ${PROCESSED_COUNT}"
echo "Completed: $(date)"

if [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
    echo "Failed samples: ${#FAILED_SAMPLES[@]}"
    printf '  - %s\n' "${FAILED_SAMPLES[@]}"
    echo "========================================="
    exit 1
else
    echo "All samples completed successfully!"
    echo "========================================="
    exit 0
fi
