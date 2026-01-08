#!/bin/bash

# ============================================
# TE Analysis Pipeline - Array Job Processing
# ============================================
# Process multiple samples in parallel using SLURM job arrays
#
# Usage:
#   1. Create sample_list.txt with format: sample_name fq1_path fq2_path rel_path
#   2. Submit job: sbatch process_array.sh
#   3. Override SBATCH settings via environment variables:
#      sbatch --export=ALL,BATCH_TIME=48:00:00,BATCH_MEM=64G process_array.sh
#
# Configuration:
#   - Sources config.sh (or config_template.sh as fallback) for pipeline settings
#   - Optionally sources batch-config.sh for SBATCH job settings
#   - Can override SBATCH settings via environment variables at submission time
#
# SBATCH Settings (can be overridden via env vars or batch-config.sh):
#   BATCH_TIME          - Job time limit (default: 24:00:00)
#   BATCH_MEM           - Memory per job (default: 32G)
#   BATCH_CPUS          - CPUs per task (default: 8)
#   BATCH_ARRAY         - Array job range (default: 1-400)
#   BATCH_ARRAY_THROTTLE - Max concurrent jobs (default: 20)
#   BATCH_LOG_DIR       - Log directory (default: logs)
# ============================================

#SBATCH --job-name=TE_array
#SBATCH --output=logs/TE_array_%A_%a.out
#SBATCH --error=logs/TE_array_%A_%a.err
#SBATCH --array=1-400%20
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# Enable strict error handling
set -euo pipefail

# ============================================
# Load Configuration
# ============================================

# Determine script directory (robust to symlinks)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source main configuration (required)
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
    echo "Loading configuration from config.sh..."
    source "${SCRIPT_DIR}/config.sh"
elif [[ -f "${SCRIPT_DIR}/config_template.sh" ]]; then
    echo "WARNING: config.sh not found, using config_template.sh as fallback"
    echo "Please copy config_template.sh to config.sh and customize it for your setup"
    source "${SCRIPT_DIR}/config_template.sh"
else
    echo "ERROR: Neither config.sh nor config_template.sh found in ${SCRIPT_DIR}"
    echo "Please create config.sh from config_template.sh with your settings"
    exit 1
fi

# Source batch configuration (optional)
if [[ -f "${SCRIPT_DIR}/batch-config.sh" ]]; then
    echo "Loading batch configuration from batch-config.sh..."
    source "${SCRIPT_DIR}/batch-config.sh"
else
    echo "INFO: batch-config.sh not found, using default SBATCH settings"
    echo "To customize SBATCH settings, copy batch-config.sh.template to batch-config.sh"
fi

# ============================================
# Validate Required Configuration
# ============================================

# Check required variables from config
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
    echo "Please check your config.sh file"
    exit 1
fi

# ============================================
# Setup
# ============================================

# Load modules if needed
if command -v module &>/dev/null; then
    module load singularity 2>/dev/null || true
fi

# Sample list file
SAMPLE_LIST="${SAMPLE_LIST:-sample_list.txt}"

# Create logs directory
mkdir -p logs

# Validate sample list exists
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST}"
    exit 1
fi

# ============================================
# Parse Sample Information from Array Task ID
# ============================================

# Get sample info from array task ID (skip header line)
SAMPLE_INFO=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_LIST}")

if [[ -z "${SAMPLE_INFO}" ]]; then
    echo "ERROR: No sample information found for array task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Parse sample information
SAMPLE_NAME=$(echo "${SAMPLE_INFO}" | awk '{print $1}')
FQ1=$(echo "${SAMPLE_INFO}" | awk '{print $2}')
FQ2=$(echo "${SAMPLE_INFO}" | awk '{print $3}')
REL_PATH=$(echo "${SAMPLE_INFO}" | awk '{print $4}')

# Validate parsed information
if [[ -z "${SAMPLE_NAME}" || -z "${FQ1}" || -z "${FQ2}" ]]; then
    echo "ERROR: Failed to parse sample information from: ${SAMPLE_INFO}"
    exit 1
fi

# ============================================
# Source Shared Functions
# ============================================

source "${SCRIPT_DIR}/function.sh"

# ============================================
# Setup Sample Output Directory
# ============================================

# Create output directory for this sample
dataDir="${OUTPUT_BASE}/${REL_PATH}/${SAMPLE_NAME}"
mkdir -p "${dataDir}"
mkdir -p "${dataDir}/log"
mkdir -p "${dataDir}/err"

# Set outputDir for JET functions
outputDir="${dataDir}"

# ============================================
# Process Sample
# ============================================

echo "========================================="
echo "Array Job ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing sample: ${SAMPLE_NAME}"
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

# ============================================
# Summary
# ============================================

echo ""
echo "========================================="
echo "Processing Summary for ${SAMPLE_NAME}"
echo "========================================="
echo "JET Step 1 Status: $([[ ${JET_STEP1_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
echo "JET Step 2 Status: $([[ ${JET_STEP2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
echo "TEProf2 Status: $([[ ${TEPROF2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
echo "========================================="

# Exit with error if any tool failed
if [[ ${JET_STEP1_STATUS} -ne 0 ]] || [[ ${JET_STEP2_STATUS} -ne 0 ]] || [[ ${TEPROF2_STATUS} -ne 0 ]]; then
    echo "ERROR: One or more pipeline steps failed for ${SAMPLE_NAME}"
    exit 1
fi

echo "SUCCESS: All pipeline steps completed for ${SAMPLE_NAME}"
exit 0
