#!/bin/bash
set -euo pipefail

# ============================================
# TE Analysis Pipeline - Independent Batch Job Submission
# ============================================
# Submit independent SLURM jobs for batch processing (NO JOB ARRAYS)
#
# Usage:
#   ./process_array.sh --sample_list_dir <dir> [OPTIONS]
#
# Options:
#   --sample_list_dir DIR    Directory containing *.list files (required)
#   --partition NAME         SLURM partition (default: from batch-config.sh)
#   --job_name PREFIX        Job name prefix (default: TE_batch)
#   --force                  Resubmit all jobs, even if already tracked
#   --help                   Show help message
#
# Output:
#   Creates/updates <sample_list_dir>/submitted_jobs.tsv with:
#   list_file  job_id  submit_time
# ============================================

# Determine script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
SAMPLE_LIST_DIR=""
PARTITION=""
JOB_NAME="TE_batch"
FORCE=false
TOOLS=""

# Usage function
usage() {
    cat << EOF
Usage: $0 --sample_list_dir <dir> [OPTIONS]

Submit independent SLURM jobs for batch sample processing.

OPTIONS:
    --sample_list_dir DIR    Directory containing *.list files (required)
    --partition NAME         SLURM partition name (optional)
    --job_name PREFIX        Job name prefix (default: TE_batch)
    --tools TOOL [TOOL...]   Select which tools to run: JET2, TEProf2, or both
                             (default: run all tools)
    --force                  Resubmit all jobs, ignoring tracking file
    --help                   Show this help message

EXAMPLES:
    # Submit jobs for all batch files in directory
    $0 --sample_list_dir sample_lists/L1/30x

    # Submit to specific partition
    $0 --sample_list_dir sample_lists/L1/30x --partition compute

    # Run only JET2 pipeline
    $0 --sample_list_dir sample_lists/L1/30x --tools JET2

    # Run only TEProf2 pipeline
    $0 --sample_list_dir sample_lists/L1/30x --tools TEProf2

    # Run both pipelines (explicit)
    $0 --sample_list_dir sample_lists/L1/30x --tools JET2 TEProf2

    # Force resubmission of all jobs
    $0 --sample_list_dir sample_lists/L1/30x --force

OUTPUT:
    Creates/updates <sample_list_dir>/submitted_jobs.tsv with:
    list_file	job_id	submit_time

EOF
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample_list_dir)
            SAMPLE_LIST_DIR="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --job_name)
            JOB_NAME="$2"
            shift 2
            ;;
        --tools)
            shift
            # Collect all tool names until we hit another option or run out of args
            TOOLS=""
            while [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-- ]]; do
                if [[ "$1" == "JET2" ]] || [[ "$1" == "TEProf2" ]]; then
                    TOOLS="${TOOLS}${TOOLS:+ }$1"
                    shift
                else
                    echo "ERROR: Invalid tool name: $1"
                    echo "Valid tools are: JET2, TEProf2"
                    exit 1
                fi
            done
            if [[ -z "${TOOLS}" ]]; then
                echo "ERROR: --tools requires at least one tool name"
                echo "Valid tools are: JET2, TEProf2"
                exit 1
            fi
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --help|-h)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "${SAMPLE_LIST_DIR}" ]]; then
    echo "ERROR: --sample_list_dir is required"
    echo "Use --help for usage information"
    exit 1
fi

if [[ ! -d "${SAMPLE_LIST_DIR}" ]]; then
    echo "ERROR: Sample list directory not found: ${SAMPLE_LIST_DIR}"
    exit 1
fi

# Load configuration
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
    echo "Loading configuration from config.sh..."
    source "${SCRIPT_DIR}/config.sh"
elif [[ -f "${SCRIPT_DIR}/config_template.sh" ]]; then
    echo "WARNING: config.sh not found, using config_template.sh"
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

# Apply partition from batch-config.sh if not specified on command line
if [[ -z "${PARTITION}" && -n "${SLURM_PARTITION:-}" ]]; then
    PARTITION="${SLURM_PARTITION}"
fi

# Get batch settings with defaults
BATCH_TIME="${BATCH_TIME:-24:00:00}"
BATCH_MEM="${BATCH_MEM:-32G}"
BATCH_CPUS="${BATCH_CPUS:-8}"

echo "========================================="
echo "Batch Job Submission"
echo "========================================="
echo "Sample list directory: ${SAMPLE_LIST_DIR}"
echo "Job name prefix: ${JOB_NAME}"
[[ -n "${PARTITION}" ]] && echo "Partition: ${PARTITION}" || echo "Partition: (default)"
echo "Time limit: ${BATCH_TIME}"
echo "Memory: ${BATCH_MEM}"
echo "CPUs: ${BATCH_CPUS}"
[[ -n "${TOOLS}" ]] && echo "Tools to run: ${TOOLS}" || echo "Tools to run: all (JET2 and TEProf2)"
echo "Force resubmit: ${FORCE}"
echo "========================================="
echo ""

# Job tracking file
JOB_TRACKING_FILE="${SAMPLE_LIST_DIR}/submitted_jobs.tsv"

# Create tracking file with header if it doesn't exist
if [[ ! -f "${JOB_TRACKING_FILE}" ]]; then
    echo -e "list_file\tjob_id\tsubmit_time" > "${JOB_TRACKING_FILE}"
    echo "Created job tracking file: ${JOB_TRACKING_FILE}"
fi

# Find all .list files (use LC_ALL=C for reproducible sorting)
mapfile -t LIST_FILES < <(LC_ALL=C find "${SAMPLE_LIST_DIR}" -maxdepth 1 -name "*.list" -type f | LC_ALL=C sort)

if [[ ${#LIST_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .list files found in ${SAMPLE_LIST_DIR}"
    exit 1
fi

echo "Found ${#LIST_FILES[@]} batch files to process"
echo ""

# Create logs directory
mkdir -p "${SAMPLE_LIST_DIR}/logs"

# Function to parse job ID from sbatch output
# Handles various formats:
# - "Submitted batch job 12345"
# - "12345"
# - "Submitted batch job 12345 on cluster foo"
parse_job_id() {
    local output="$1"
    # Try to extract numeric job ID from various formats
    if [[ "${output}" =~ [[:space:]]([0-9]+)([[:space:]]|$) ]]; then
        echo "${BASH_REMATCH[1]}"
        return 0
    elif [[ "${output}" =~ ^([0-9]+)$ ]]; then
        echo "${output}"
        return 0
    else
        echo ""
        return 1
    fi
}

# Submit jobs
SUBMITTED_COUNT=0
SKIPPED_COUNT=0
FAILED_COUNT=0

for list_file in "${LIST_FILES[@]}"; do
    list_basename=$(basename "${list_file}")
    
    # Check if already submitted (unless --force)
    if [[ "${FORCE}" == "false" ]]; then
        if grep -q "^${list_basename}" "${JOB_TRACKING_FILE}" 2>/dev/null; then
            echo "SKIP: ${list_basename} (already tracked, use --force to resubmit)"
            SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
            continue
        fi
    fi
    
    # Count samples in this batch
    sample_count=$(grep -v "^#" "${list_file}" | wc -l)
    
    echo "Submitting job for ${list_basename} (${sample_count} samples)..."
    
    # Build sbatch command
    sbatch_cmd="sbatch"
    sbatch_cmd="${sbatch_cmd} --job-name=${JOB_NAME}_${list_basename%.list}"
    sbatch_cmd="${sbatch_cmd} --output=${SAMPLE_LIST_DIR}/logs/${list_basename%.list}_%j.out"
    sbatch_cmd="${sbatch_cmd} --error=${SAMPLE_LIST_DIR}/logs/${list_basename%.list}_%j.err"
    sbatch_cmd="${sbatch_cmd} --time=${BATCH_TIME}"
    sbatch_cmd="${sbatch_cmd} --mem=${BATCH_MEM}"
    sbatch_cmd="${sbatch_cmd} --cpus-per-task=${BATCH_CPUS}"
    
    # Add partition if specified
    if [[ -n "${PARTITION}" ]]; then
        sbatch_cmd="${sbatch_cmd} --partition=${PARTITION}"
    fi
    
    # Export TOOLS variable if specified
    if [[ -n "${TOOLS}" ]]; then
        sbatch_cmd="${sbatch_cmd} --export=ALL,TOOLS='${TOOLS}'"
    fi
    
    # Add worker script and sample list file
    sbatch_cmd="${sbatch_cmd} ${SCRIPT_DIR}/process_batch_worker.sh ${list_file}"
    
    # Submit job
    sbatch_output=$(eval "${sbatch_cmd}" 2>&1)
    submit_status=$?
    
    if [[ ${submit_status} -eq 0 ]]; then
        # Parse job ID from output
        job_id=$(parse_job_id "${sbatch_output}")
        
        if [[ -n "${job_id}" ]]; then
            # Record submission
            submit_time=$(date '+%Y-%m-%d %H:%M:%S')
            echo -e "${list_basename}\t${job_id}\t${submit_time}" >> "${JOB_TRACKING_FILE}"
            
            echo "  ✓ Submitted job ${job_id}"
            SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
        else
            echo "  ✗ Failed to parse job ID from: ${sbatch_output}"
            FAILED_COUNT=$((FAILED_COUNT + 1))
        fi
    else
        echo "  ✗ sbatch failed: ${sbatch_output}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
done

echo ""
echo "========================================="
echo "Submission Summary"
echo "========================================="
echo "Total batch files: ${#LIST_FILES[@]}"
echo "Submitted: ${SUBMITTED_COUNT}"
echo "Skipped: ${SKIPPED_COUNT}"
echo "Failed: ${FAILED_COUNT}"
echo "========================================="
echo ""

if [[ ${SUBMITTED_COUNT} -gt 0 ]]; then
    echo "Job tracking file: ${JOB_TRACKING_FILE}"
    echo ""
    echo "Check job status with:"
    echo "  squeue -u \${USER}"
    echo "  scripts/check_status.sh --sample_list_dir ${SAMPLE_LIST_DIR}"
    echo ""
fi

if [[ ${FAILED_COUNT} -gt 0 ]]; then
    echo "WARNING: ${FAILED_COUNT} job submissions failed"
    exit 1
fi

exit 0


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

# Validate that we have samples to process
if [[ ${START_SAMPLE} -gt ${TOTAL_SAMPLES} ]]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} is beyond the range of available samples"
    echo "Total samples: ${TOTAL_SAMPLES}"
    echo "Calculated start sample: ${START_SAMPLE}"
    echo "Please adjust the SLURM array range to match the number of jobs needed"
    exit 1
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
