#!/bin/bash
# Helper script to calculate SLURM array size for sample chunking
# Usage: ./calculate_array_size.sh [samples_per_job]

set -euo pipefail

# Default sample list file
SAMPLE_LIST="${SAMPLE_LIST:-sample_list.txt}"

# Check if sample list exists
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST}"
    echo "Please generate it first with: ./generate_sample_list.sh"
    exit 1
fi

# Count total samples (excluding header line)
TOTAL_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l)

# Get samples per job from argument or use default
SAMPLES_PER_JOB="${1:-1}"

# Validate input
if [[ ! "${SAMPLES_PER_JOB}" =~ ^[0-9]+$ ]] || [[ ${SAMPLES_PER_JOB} -lt 1 ]]; then
    echo "ERROR: SAMPLES_PER_JOB must be a positive integer"
    echo "Usage: $0 [samples_per_job]"
    exit 1
fi

# Calculate number of jobs needed (ceiling division)
NUM_JOBS=$(( (TOTAL_SAMPLES + SAMPLES_PER_JOB - 1) / SAMPLES_PER_JOB ))

# Display results
echo "========================================="
echo "SLURM Array Size Calculator"
echo "========================================="
echo "Sample list: ${SAMPLE_LIST}"
echo "Total samples: ${TOTAL_SAMPLES}"
echo "Samples per job: ${SAMPLES_PER_JOB}"
echo "Jobs needed: ${NUM_JOBS}"
echo "========================================="
echo ""
echo "Configuration for process_array.sh:"
echo "-----------------------------------"
echo "Update the #SBATCH --array line to:"
echo "  #SBATCH --array=1-${NUM_JOBS}%20"
echo ""
echo "Or with custom throttle (e.g., 10 concurrent jobs):"
echo "  #SBATCH --array=1-${NUM_JOBS}%10"
echo ""
echo "Set SAMPLES_PER_JOB in batch-config.sh:"
echo "  export SAMPLES_PER_JOB=\"${SAMPLES_PER_JOB}\""
echo ""
echo "Or override at submission time:"
echo "  sbatch --export=ALL,SAMPLES_PER_JOB=${SAMPLES_PER_JOB} process_array.sh"
echo "========================================="
echo ""

# Show sample distribution
echo "Sample Distribution Across Jobs:"
echo "-----------------------------------"
for ((JOB_ID=1; JOB_ID<=NUM_JOBS && JOB_ID<=5; JOB_ID++)); do
    START_SAMPLE=$(( (JOB_ID - 1) * SAMPLES_PER_JOB + 1 ))
    END_SAMPLE=$(( JOB_ID * SAMPLES_PER_JOB ))
    
    # Don't exceed total samples
    if [[ ${END_SAMPLE} -gt ${TOTAL_SAMPLES} ]]; then
        END_SAMPLE=${TOTAL_SAMPLES}
    fi
    
    NUM_IN_JOB=$(( END_SAMPLE - START_SAMPLE + 1 ))
    echo "Job ${JOB_ID}: samples ${START_SAMPLE}-${END_SAMPLE} (${NUM_IN_JOB} samples)"
done

if [[ ${NUM_JOBS} -gt 5 ]]; then
    echo "... (showing first 5 jobs only)"
    
    # Show last job
    JOB_ID=${NUM_JOBS}
    START_SAMPLE=$(( (JOB_ID - 1) * SAMPLES_PER_JOB + 1 ))
    END_SAMPLE=$(( JOB_ID * SAMPLES_PER_JOB ))
    if [[ ${END_SAMPLE} -gt ${TOTAL_SAMPLES} ]]; then
        END_SAMPLE=${TOTAL_SAMPLES}
    fi
    NUM_IN_JOB=$(( END_SAMPLE - START_SAMPLE + 1 ))
    echo "Job ${JOB_ID}: samples ${START_SAMPLE}-${END_SAMPLE} (${NUM_IN_JOB} samples)"
fi

echo "========================================="
