#!/bin/bash
set -euo pipefail

# ============================================
# Job Status Report
# ============================================
# Query SLURM for job status based on submitted jobs tracking file
#
# Usage:
#   ./job_status_report.sh --job_list <tsv_file>
#   ./job_status_report.sh --sample_list_dir <dir>
# ============================================

# Default values
JOB_LIST=""
SAMPLE_LIST_DIR=""

# Usage function
usage() {
    cat << EOF
Usage: $0 --job_list <tsv_file>
   or: $0 --sample_list_dir <dir>

Query SLURM for job status and report summary.

OPTIONS:
    --job_list FILE          Path to submitted_jobs.tsv file
    --sample_list_dir DIR    Directory containing submitted_jobs.tsv
    --help                   Show this help message

EXAMPLES:
    # Check status using job list file
    $0 --job_list sample_lists/L1/30x/submitted_jobs.tsv

    # Check status using sample list directory
    $0 --sample_list_dir sample_lists/L1/30x

OUTPUT:
    Reports counts of jobs by status:
    - Total jobs
    - Running (R)
    - Pending (PD)
    - Completed (CD)
    - Failed (F, CA, TO)
    - Other states

    Exit codes:
    0 - All jobs completed successfully
    1 - Some jobs failed
    2 - Some jobs still running/pending

EOF
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --job_list)
            JOB_LIST="$2"
            shift 2
            ;;
        --sample_list_dir)
            SAMPLE_LIST_DIR="$2"
            shift 2
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

# Determine job list file
if [[ -n "${SAMPLE_LIST_DIR}" ]]; then
    JOB_LIST="${SAMPLE_LIST_DIR}/submitted_jobs.tsv"
elif [[ -z "${JOB_LIST}" ]]; then
    echo "ERROR: Either --job_list or --sample_list_dir is required"
    echo "Use --help for usage information"
    exit 1
fi

# Validate job list file exists
if [[ ! -f "${JOB_LIST}" ]]; then
    echo "ERROR: Job list file not found: ${JOB_LIST}"
    exit 1
fi

echo "========================================="
echo "Job Status Report"
echo "========================================="
echo "Job list: ${JOB_LIST}"
echo "Generated: $(date)"
echo "========================================="
echo ""

# Read job IDs from tracking file (skip header)
mapfile -t JOB_ENTRIES < <(tail -n +2 "${JOB_LIST}")

if [[ ${#JOB_ENTRIES[@]} -eq 0 ]]; then
    echo "No jobs found in tracking file"
    exit 0
fi

TOTAL_JOBS=${#JOB_ENTRIES[@]}

# Initialize counters
RUNNING=0
PENDING=0
COMPLETED=0
FAILED=0
CANCELLED=0
TIMEOUT=0
NOT_FOUND=0
OTHER=0

# Collect job IDs
declare -a JOB_IDS
for entry in "${JOB_ENTRIES[@]}"; do
    # Parse TSV: list_file, job_id, submit_time
    job_id=$(echo "${entry}" | awk -F'\t' '{print $2}')
    if [[ -n "${job_id}" ]]; then
        JOB_IDS+=("${job_id}")
    fi
done

# Query SLURM for each job
# First try squeue for running/pending jobs
declare -A JOB_STATUS

echo "Querying SLURM for job status..."

# Get status from squeue (running/pending jobs)
if command -v squeue &>/dev/null; then
    for job_id in "${JOB_IDS[@]}"; do
        status=$(squeue -j "${job_id}" -h -o "%T" 2>/dev/null || echo "")
        if [[ -n "${status}" ]]; then
            JOB_STATUS["${job_id}"]="${status}"
        fi
    done
fi

# Get status from sacct for completed/failed jobs (not in squeue)
if command -v sacct &>/dev/null; then
    for job_id in "${JOB_IDS[@]}"; do
        # Skip if already found in squeue
        if [[ -n "${JOB_STATUS[${job_id}]:-}" ]]; then
            continue
        fi
        
        # Query sacct
        status=$(sacct -j "${job_id}" -n -o "State%20" 2>/dev/null | head -1 | awk '{print $1}')
        if [[ -n "${status}" ]]; then
            JOB_STATUS["${job_id}"]="${status}"
        else
            JOB_STATUS["${job_id}"]="NOT_FOUND"
        fi
    done
fi

# Count status categories
for job_id in "${JOB_IDS[@]}"; do
    status="${JOB_STATUS[${job_id}]:-NOT_FOUND}"
    
    case "${status}" in
        RUNNING|R)
            RUNNING=$((RUNNING + 1))
            ;;
        PENDING|PD|CONFIGURING|CF|REQUEUED|RQ|RESIZING|RS|SUSPENDED|S)
            PENDING=$((PENDING + 1))
            ;;
        COMPLETED|CD)
            COMPLETED=$((COMPLETED + 1))
            ;;
        FAILED|F|NODE_FAIL|NF|BOOT_FAIL|BF|DEADLINE|DL|OUT_OF_MEMORY|OOM|PREEMPTED|PR)
            FAILED=$((FAILED + 1))
            ;;
        CANCELLED|CA|CANCELLED+|CA+)
            CANCELLED=$((CANCELLED + 1))
            ;;
        TIMEOUT|TO)
            TIMEOUT=$((TIMEOUT + 1))
            ;;
        NOT_FOUND)
            NOT_FOUND=$((NOT_FOUND + 1))
            ;;
        *)
            OTHER=$((OTHER + 1))
            ;;
    esac
done

# Display summary
echo "Status Summary:"
echo "  Total jobs: ${TOTAL_JOBS}"
echo "  Running: ${RUNNING}"
echo "  Pending: ${PENDING}"
echo "  Completed: ${COMPLETED}"
echo "  Failed: ${FAILED}"
echo "  Cancelled: ${CANCELLED}"
echo "  Timeout: ${TIMEOUT}"
echo "  Not found: ${NOT_FOUND}"
[[ ${OTHER} -gt 0 ]] && echo "  Other: ${OTHER}"
echo ""

# Show detailed status for failed/problematic jobs
if [[ ${FAILED} -gt 0 ]] || [[ ${CANCELLED} -gt 0 ]] || [[ ${TIMEOUT} -gt 0 ]]; then
    echo "Problematic Jobs:"
    echo "----------------------------------------"
    
    line_num=0
    for entry in "${JOB_ENTRIES[@]}"; do
        list_file=$(echo "${entry}" | awk -F'\t' '{print $1}')
        job_id=$(echo "${entry}" | awk -F'\t' '{print $2}')
        status="${JOB_STATUS[${job_id}]:-UNKNOWN}"
        
        case "${status}" in
            FAILED|F|NODE_FAIL|NF|BOOT_FAIL|BF|DEADLINE|DL|OUT_OF_MEMORY|OOM|PREEMPTED|PR|CANCELLED|CA|CANCELLED+|CA+|TIMEOUT|TO)
                echo "  Job ${job_id} (${list_file}): ${status}"
                ;;
        esac
    done
    echo ""
fi

# Calculate percentages
if [[ ${TOTAL_JOBS} -gt 0 ]]; then
    completed_pct=$(awk "BEGIN {printf \"%.1f\", ${COMPLETED}/${TOTAL_JOBS}*100}")
    echo "Progress: ${COMPLETED}/${TOTAL_JOBS} completed (${completed_pct}%)"
fi

echo "========================================="

# Determine exit code
if [[ ${FAILED} -gt 0 ]] || [[ ${CANCELLED} -gt 0 ]] || [[ ${TIMEOUT} -gt 0 ]]; then
    echo "Status: FAILED - Some jobs have failed"
    exit 1
elif [[ ${RUNNING} -gt 0 ]] || [[ ${PENDING} -gt 0 ]]; then
    echo "Status: IN PROGRESS - Jobs still running or pending"
    exit 2
elif [[ ${COMPLETED} -eq ${TOTAL_JOBS} ]]; then
    echo "Status: SUCCESS - All jobs completed"
    exit 0
else
    echo "Status: UNKNOWN - Cannot determine overall status"
    exit 2
fi
