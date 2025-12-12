#!/bin/bash

#SBATCH --job-name=TE_array
#SBATCH --output=logs/TE_array_%A_%a.out
#SBATCH --error=logs/TE_array_%A_%a.err
#SBATCH --array=1-400%20
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# TE Analysis Pipeline - Array Job Processing
# Configuration is loaded from the shared config.sh file

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the shared configuration file
if [ -f "${SCRIPT_DIR}/config.sh" ]; then
    source "${SCRIPT_DIR}/config.sh"
else
    echo "ERROR: Configuration file not found: ${SCRIPT_DIR}/config.sh"
    echo "Please create config.sh from config_template.sh and update the paths."
    exit 1
fi

module load singularity

# Create logs directory
mkdir -p logs

# Get sample info from array task ID
SAMPLE_INFO=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${SAMPLE_LIST})

# Parse sample information
SAMPLE_NAME=$(echo ${SAMPLE_INFO} | awk '{print $1}')
FQ1=$(echo ${SAMPLE_INFO} | awk '{print $2}')
FQ2=$(echo ${SAMPLE_INFO} | awk '{print $3}')
REL_PATH=$(echo ${SAMPLE_INFO} | awk '{print $4}')

# Create output directory for this sample
dataDir="${OUTPUT_BASE}/${REL_PATH}/${SAMPLE_NAME}"
mkdir -p "${dataDir}"
mkdir -p "${dataDir}/log"
mkdir -p "${dataDir}/err"

echo "========================================="
echo "Array Job ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing sample: ${SAMPLE_NAME}"
echo "Input files:"
echo "  R1: ${FQ1}"
echo "  R2: ${FQ2}"
echo "Output directory: ${dataDir}"
echo "========================================="

# Create metadata file for this sample
metaFile="${dataDir}/metadata.txt"
echo -e "sample\tfastq1\tfastq2" > ${metaFile}
echo -e "${SAMPLE_NAME}\t${FQ1}\t${FQ2}" >> ${metaFile}

# Source shared functions
source "${SCRIPT_DIR}/function.sh"

# Run JET Step 1
run_jet_step1
JET_STEP1_STATUS=$?

# Run JET Step 2 (only if Step 1 succeeded)
if [ ${JET_STEP1_STATUS} -eq 0 ]; then
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
echo "JET Step 1 Status: $([ ${JET_STEP1_STATUS} -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "JET Step 2 Status: $([ ${JET_STEP2_STATUS} -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "TEProf2 Status: $([ ${TEPROF2_STATUS} -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "========================================="

# Exit with error if any tool failed
if [ ${JET_STEP1_STATUS} -ne 0 ] || [ ${JET_STEP2_STATUS} -ne 0 ] || [ ${TEPROF2_STATUS} -ne 0 ]; then
    exit 1
fi

exit 0
