#!/bin/bash

# ============================================
# TE Analysis Pipeline - Sequential Sample Processing
# ============================================
# Process all samples sequentially, followed by TEProf2 aggregation
#
# Usage:
#   1. Configure paths in config.sh (copy from config_template.sh)
#   2. Submit job: sbatch process_samples.sh
#   3. Override SBATCH settings via environment variables:
#      sbatch --export=ALL,BATCH_TIME=48:00:00,BATCH_MEM=64G process_samples.sh
#
# Configuration:
#   - Sources config.sh (or config_template.sh as fallback) for pipeline settings
#   - Optionally sources batch-config.sh for SBATCH job settings
#   - Can override SBATCH settings via environment variables at submission time
#
# SBATCH Settings (can be overridden via env vars or batch-config.sh):
#   BATCH_TIME          - Job time limit (default: 48:00:00)
#   BATCH_MEM           - Memory per job (default: 32G)
#   BATCH_CPUS          - CPUs per task (default: 8)
#   BATCH_LOG_DIR       - Log directory (default: logs)
# ============================================

#SBATCH --job-name=TE_analysis
#SBATCH --output=logs/TE_analysis_%A_%a.out
#SBATCH --error=logs/TE_analysis_%A_%a.err
#SBATCH --time=48:00:00
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

# Create logs directory
mkdir -p logs

# Source shared functions
source "${SCRIPT_DIR}/function.sh"

# ============================================
# Helper Functions
# ============================================

# Function to extract sample name from fastq file
get_sample_name() {
    local fq_file=$1
    local basename=$(basename "$fq_file")
    # Remove .1.fq.gz or .2.fq.gz extension
    echo "${basename%.*.fq.gz}"
}

# Function to process a sample pair
process_sample() {
    local fq1=$1
    local fq2=$2
    local rel_path=$3  # Relative path from DATA_HOME
    
    # Extract sample name
    local sample_name=$(get_sample_name "$fq1")
    
    # Set up global variables required by shared functions
    SAMPLE_NAME="${sample_name}"
    FQ1="${fq1}"
    FQ2="${fq2}"
    
    # Create output directory maintaining original structure
    dataDir="${OUTPUT_BASE}/${rel_path}/${sample_name}"
    mkdir -p "${dataDir}"
    mkdir -p "${dataDir}/log"
    mkdir -p "${dataDir}/err"
    
    # Set outputDir for JET functions
    outputDir="${dataDir}"
    
    echo "========================================="
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
    
    # Summary
    echo ""
    echo "========================================="
    echo "Processing Summary for ${SAMPLE_NAME}"
    echo "========================================="
    echo "JET Step 1 Status: $([[ ${JET_STEP1_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "JET Step 2 Status: $([[ ${JET_STEP2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "TEProf2 Status: $([[ ${TEPROF2_STATUS} -eq 0 ]] && echo 'SUCCESS' || echo 'FAILED')"
    echo "========================================="
    
    echo "Completed processing ${SAMPLE_NAME}"
    echo ""
}

# ============================================
# Main Processing Loop
# ============================================

echo "Starting TE analysis pipeline..."
echo "Data home: ${DATA_HOME}"
echo "Output base: ${OUTPUT_BASE}"
echo ""

# Process nonReferenceTE samples
for te_type in AluY L1HS LTR5 SVA_F; do
    for coverage in 5X 10X 50X 100X 200X; do
        fq_dir="${DATA_HOME}/nonReferenceTE/${te_type}/${coverage}/fq"
        
        if [[ -d "$fq_dir" ]]; then
            echo "Processing ${te_type} ${coverage}..."
            
            # Get all R1 files
            for fq1 in "${fq_dir}"/*.1.fq.gz; do
                if [[ -f "$fq1" ]]; then
                    # Construct R2 filename
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    
                    if [[ -f "$fq2" ]]; then
                        # Relative path for output
                        rel_path="nonReferenceTE/${te_type}/${coverage}"
                        process_sample "$fq1" "$fq2" "$rel_path"
                    else
                        echo "Warning: Missing R2 file for $fq1"
                    fi
                fi
            done
        fi
    done
done

# Process referenceTE/intron samples
for coverage in 5X 10X 50X 100X 200X; do
    # Process regular intron samples
    fq_dir="${DATA_HOME}/referenceTE/intron/${coverage}/fq"
    if [[ -d "$fq_dir" ]]; then
        echo "Processing referenceTE/intron ${coverage} (regular)..."
        
        for fq1 in "${fq_dir}"/reftefu403_*_${coverage}.1.fq.gz; do
            if [[ -f "$fq1" ]]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [[ -f "$fq2" ]]; then
                    rel_path="referenceTE/intron/${coverage}/fq"
                    process_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
    
    # Process mutated intron samples
    fq_mut_dir="${DATA_HOME}/referenceTE/intron/${coverage}/fq_mut"
    if [[ -d "$fq_mut_dir" ]]; then
        echo "Processing referenceTE/intron ${coverage} (mutated)..."
        
        for fq1 in "${fq_mut_dir}"/reftefu403_mut_*_${coverage}.1.fq.gz; do
            if [[ -f "$fq1" ]]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [[ -f "$fq2" ]]; then
                    rel_path="referenceTE/intron/${coverage}/fq_mut"
                    process_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
done

# Process referenceTE/TSS samples
for coverage in 5X 10X 50X 100X 200X; do
    # Process regular TSS samples
    fq_dir="${DATA_HOME}/referenceTE/TSS/${coverage}/fq"
    if [[ -d "$fq_dir" ]]; then
        echo "Processing referenceTE/TSS ${coverage} (regular)..."
        
        for fq1 in "${fq_dir}"/reftetss270_*_${coverage}.1.fq.gz; do
            if [[ -f "$fq1" ]]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [[ -f "$fq2" ]]; then
                    rel_path="referenceTE/TSS/${coverage}/fq"
                    process_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
    
    # Process mutated TSS samples
    fq_mut_dir="${DATA_HOME}/referenceTE/TSS/${coverage}/fq_mut"
    if [[ -d "$fq_mut_dir" ]]; then
        echo "Processing referenceTE/TSS ${coverage} (mutated)..."
        
        for fq1 in "${fq_mut_dir}"/reftetss270_mut_*_${coverage}.1.fq.gz; do
            if [[ -f "$fq1" ]]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [[ -f "$fq2" ]]; then
                    rel_path="referenceTE/TSS/${coverage}/fq_mut"
                    process_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
done

echo "========================================="
echo "All individual samples processed!"
echo "========================================="

# ============================================
# TEProf2 Aggregation Step
# ============================================
# After all individual samples are processed, run TEProf2 aggregation
# This step executes run_command.sh to aggregate results across all samples
# and identify TE-derived transcripts

echo ""
echo "========================================="
echo "Starting TEProf2 Aggregation"
echo "========================================="

# Set up aggregation environment
export TEPROF2_AGGREGATION_DIR="${OUTPUT_BASE}/TEProf2_aggregated"
export RUN_COMMAND_SCRIPT="${SCRIPT_DIR}/run_command.sh"

# Note: dataDir should point to the root output directory containing all TEProf2 sample outputs
export dataDir="${OUTPUT_BASE}"

# Run aggregation
run_teprof2_aggregation
AGGREGATION_STATUS=$?

if [[ ${AGGREGATION_STATUS} -eq 0 ]]; then
    echo ""
    echo "========================================="
    echo "TEProf2 Aggregation completed successfully!"
    echo "Results are in: ${TEPROF2_AGGREGATION_DIR}"
    echo "========================================="
else
    echo ""
    echo "========================================="
    echo "TEProf2 Aggregation FAILED!"
    echo "Check logs in: ${TEPROF2_AGGREGATION_DIR}"
    echo "========================================="
fi

echo ""
echo "========================================="
echo "Pipeline execution completed!"
echo "========================================="
