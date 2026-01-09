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
#   SAMPLES_PER_JOB     - Number of samples per job (default: 1, enables chunking)
#
# Sample Chunking:
#   Set SAMPLES_PER_JOB to process multiple samples per job sequentially.
#   Example: 100 samples with SAMPLES_PER_JOB=20 will create 5 jobs, each processing 20 samples.
#   Usage: sbatch --export=ALL,SAMPLES_PER_JOB=20 process_array.sh
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
# Sample Chunking Configuration
# ============================================

# Number of samples to process per job (default: 1, no chunking)
SAMPLES_PER_JOB="${SAMPLES_PER_JOB:-1}"

# Calculate which samples this job should process
# Count total samples (excluding header line)
TOTAL_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l)

# Calculate start and end sample indices for this job
START_SAMPLE=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_JOB + 1 ))
END_SAMPLE=$(( SLURM_ARRAY_TASK_ID * SAMPLES_PER_JOB ))

# Don't exceed total samples
if [[ ${END_SAMPLE} -gt ${TOTAL_SAMPLES} ]]; then
    END_SAMPLE=${TOTAL_SAMPLES}
fi

echo "========================================="
echo "Array Job Configuration"
echo "========================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Samples per job: ${SAMPLES_PER_JOB}"
echo "Total samples: ${TOTAL_SAMPLES}"
echo "Processing samples: ${START_SAMPLE} to ${END_SAMPLE}"
echo "========================================="
echo ""

# ============================================
# Source Shared Functions
# ============================================

source "${SCRIPT_DIR}/function.sh"

# ============================================
# Process Samples in Chunk
# ============================================

# Track overall status
OVERALL_SUCCESS=0
FAILED_SAMPLES=()

# Loop through samples in this chunk
for SAMPLE_IDX in $(seq ${START_SAMPLE} ${END_SAMPLE}); do
    # Get sample info (skip header line, so add 1 to index)
    SAMPLE_INFO=$(sed -n "$((SAMPLE_IDX + 1))p" "${SAMPLE_LIST}")
    
    if [[ -z "${SAMPLE_INFO}" ]]; then
        echo "WARNING: No sample information found for sample index ${SAMPLE_IDX}"
        continue
    fi
    
    # Parse sample information
    SAMPLE_NAME=$(echo "${SAMPLE_INFO}" | awk '{print $1}')
    FQ1=$(echo "${SAMPLE_INFO}" | awk '{print $2}')
    FQ2=$(echo "${SAMPLE_INFO}" | awk '{print $3}')
    REL_PATH=$(echo "${SAMPLE_INFO}" | awk '{print $4}')
    
    # Validate parsed information
    if [[ -z "${SAMPLE_NAME}" || -z "${FQ1}" || -z "${FQ2}" ]]; then
        echo "ERROR: Failed to parse sample information from: ${SAMPLE_INFO}"
        FAILED_SAMPLES+=("${SAMPLE_NAME:-UNKNOWN}")
        OVERALL_SUCCESS=1
        continue
    fi
    
    # ============================================
    # Setup Sample Output Directory
    # ============================================
    
    # Create output directory for this sample
    dataDir="${OUTPUT_BASE}/${REL_PATH}/${SAMPLE_NAME}"
    mkdir -p "${dataDir}"
    mkdir -p "${dataDir}/log"
    mkdir -p "${dataDir}/err"
    
    # Set outputDir for JET functions (function.sh expects this to be the per-sample directory)
    outputDir="${dataDir}"
    
    # ============================================
    # Process Sample
    # ============================================
    
    echo "========================================="
    echo "Processing sample ${SAMPLE_IDX}/${TOTAL_SAMPLES}: ${SAMPLE_NAME}"
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
    # Sample Summary
    # ============================================
    
    echo ""
    echo "========================================="
    echo "Processing Summary for ${SAMPLE_NAME}"
    echo "========================================="
    echo "JET Step 1 Status: $([[ ${JET_STEP1_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "JET Step 2 Status: $([[ ${JET_STEP2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "TEProf2 Status: $([[ ${TEPROF2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "========================================="
    echo ""
    
    # Track if any tool failed
    if [[ ${JET_STEP1_STATUS} -ne 0 ]] || [[ ${JET_STEP2_STATUS} -ne 0 ]] || [[ ${TEPROF2_STATUS} -ne 0 ]]; then
        echo "ERROR: One or more pipeline steps failed for ${SAMPLE_NAME}"
        FAILED_SAMPLES+=("${SAMPLE_NAME}")
        OVERALL_SUCCESS=1
    else
        echo "SUCCESS: All pipeline steps completed for ${SAMPLE_NAME}"
    fi
done

# ============================================
# Job Summary
# ============================================

echo ""
echo "========================================="
echo "Job Summary - Array Task ${SLURM_ARRAY_TASK_ID}"
echo "========================================="
echo "Processed samples: ${START_SAMPLE} to ${END_SAMPLE}"
echo "Total samples processed: $(( END_SAMPLE - START_SAMPLE + 1 ))"

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
