#!/bin/bash

# Standalone script to run TEProf2 aggregation
# This script should be run AFTER all individual TEProf2 samples have been processed
# It executes run_command.sh to perform cross-sample analysis and identify TE-derived transcripts

#SBATCH --job-name=TEProf2_agg
#SBATCH --output=logs/TEProf2_aggregation_%j.out
#SBATCH --error=logs/TEProf2_aggregation_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# Usage:
#   ./run_teprof2_aggregation.sh [output_base_dir]
#
# If output_base_dir is not provided, uses OUTPUT_BASE from config or defaults

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source configuration if available
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
    echo "Loading configuration from config.sh..."
    source "${SCRIPT_DIR}/config.sh"
elif [[ -f "${SCRIPT_DIR}/config_template.sh" ]]; then
    echo "WARNING: Using config_template.sh (you should create config.sh)"
    source "${SCRIPT_DIR}/config_template.sh"
else
    echo "ERROR: No configuration file found!"
    exit 1
fi

# Source shared functions
source "${SCRIPT_DIR}/function.sh"

# Override OUTPUT_BASE if provided as argument
if [[ $# -ge 1 ]]; then
    OUTPUT_BASE="$1"
    echo "Using output directory from argument: ${OUTPUT_BASE}"
fi

# Validate required variables
: "${OUTPUT_BASE:?ERROR: OUTPUT_BASE not set. Provide as argument or in config.sh}"
: "${ARGUMENTS_TXT:?ERROR: ARGUMENTS_TXT not set in config.sh}"
: "${TEPROF2_CUFFMERGE_GTF:?ERROR: TEPROF2_CUFFMERGE_GTF not set in config.sh}"

# Set up aggregation environment
export TEPROF2_AGGREGATION_DIR="${OUTPUT_BASE}/TEProf2_aggregated"
export RUN_COMMAND_SCRIPT="${SCRIPT_DIR}/run_command.sh"
export dataDir="${OUTPUT_BASE}"

echo "========================================="
echo "TEProf2 Aggregation"
echo "========================================="
echo "Output base directory: ${OUTPUT_BASE}"
echo "Aggregation directory: ${TEPROF2_AGGREGATION_DIR}"
echo "Arguments file: ${ARGUMENTS_TXT}"
echo "Run command script: ${RUN_COMMAND_SCRIPT}"
echo "========================================="

# Check that we have processed samples
# Note: These find commands may take time for large directory trees
echo "Scanning for processed samples..."
sample_count=$(find "${OUTPUT_BASE}" -type f -name "*_annotated_filtered_test_all" 2>/dev/null | wc -l)
bam_count=$(find "${OUTPUT_BASE}" -type f -name "*.bam" 2>/dev/null | wc -l)

echo "Found ${sample_count} annotated sample files"
echo "Found ${bam_count} BAM files"

if [[ ${sample_count} -eq 0 ]]; then
    echo "ERROR: No processed TEProf2 samples found in ${OUTPUT_BASE}"
    echo "Please run individual sample processing first (process_samples.sh or process_array.sh)"
    exit 1
fi

if [[ ${bam_count} -eq 0 ]]; then
    echo "WARNING: No BAM files found. run_command.sh requires BAM files for quantification."
fi

echo ""
echo "Starting aggregation..."

# Create logs directory if it doesn't exist
mkdir -p "${SCRIPT_DIR}/logs"

# Run aggregation
run_teprof2_aggregation
AGGREGATION_STATUS=$?

if [ ${AGGREGATION_STATUS} -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "TEProf2 Aggregation completed successfully!"
    echo "========================================="
    echo ""
    echo "Output files in ${TEPROF2_AGGREGATION_DIR}:"
    echo "  - filter_combined_candidates.tsv: All TE-gene transcripts"
    echo "  - initial_candidate_list.tsv: Summary of unique candidates"
    echo "  - read_filtered_candidates.tsv: Candidates passing read filters"
    echo "  - candidate_transcripts.gff3: GFF3 format of detected transcripts"
    echo "  - reference_merged_candidates.gtf: Merged reference with candidates"
    echo "  - table_frac_tot_cand: Fraction of expression per candidate"
    echo "  - table_tpm_cand: TPM values per candidate"
    echo ""
    exit 0
else
    echo ""
    echo "========================================="
    echo "TEProf2 Aggregation FAILED!"
    echo "========================================="
    echo "Check logs and intermediate files in: ${TEPROF2_AGGREGATION_DIR}"
    exit 1
fi
