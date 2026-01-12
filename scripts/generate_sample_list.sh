#!/bin/bash
set -euo pipefail

# ============================================
# Generate sample lists for batch job processing
# ============================================
# This script generates sample lists and splits them into batch files
# for independent SLURM job submission (no job arrays)
#
# Usage:
#   ./generate_sample_list.sh --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists
#   ./generate_sample_list.sh --data_home /path/to/data --outdir /path/to/lists --batch_size 10
#   ./generate_sample_list.sh --help
# ============================================

# Default values
DATA_HOME="/home/junseokp/workspaces/data/rTea-simul/sims"
OUTDIR="sample_lists"
BATCH_SIZE=""
TE_TYPE=""
COVERAGE=""

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Generate sample lists and split them into batches for SLURM job submission.

OPTIONS:
    --te_type TYPE       Filter by TE type (e.g., L1, AluY, L1HS, LTR5, SVA_F)
                         If not specified, processes all TE types
    --coverage COV       Filter by coverage (e.g., 5x, 10x, 30x, 50x, 100x, 200x)
                         If not specified, processes all coverages
    --batch_size N       Number of samples per batch file (required)
    --outdir DIR         Output directory for batch files (default: sample_lists)
    --data_home DIR      Path to data directory (default: ${DATA_HOME})
    --help               Show this help message

EXAMPLES:
    # Generate lists for L1 at 30x coverage, 20 samples per batch
    $0 --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists

    # Generate lists for all samples, 10 samples per batch
    $0 --batch_size 10 --outdir sample_lists

    # Generate lists for AluY at all coverages
    $0 --te_type AluY --batch_size 15

OUTPUT:
    Creates directory structure: <outdir>/<te_type>/<coverage>/
    Batch files named: samples_0001.list, samples_0002.list, etc.
    Each batch file contains: Sample_Name FQ1 FQ2 Output_Path

EOF
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --te_type)
            TE_TYPE="$2"
            shift 2
            ;;
        --coverage)
            COVERAGE="$2"
            shift 2
            ;;
        --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --data_home)
            DATA_HOME="$2"
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

# Validate required arguments
if [[ -z "${BATCH_SIZE}" ]]; then
    echo "ERROR: --batch_size is required"
    echo "Use --help for usage information"
    exit 1
fi

if [[ ! "${BATCH_SIZE}" =~ ^[0-9]+$ ]] || [[ ${BATCH_SIZE} -lt 1 ]]; then
    echo "ERROR: --batch_size must be a positive integer"
    exit 1
fi

# Validate data home exists
if [[ ! -d "${DATA_HOME}" ]]; then
    echo "ERROR: Data directory not found: ${DATA_HOME}"
    exit 1
fi

# Normalize coverage format (add X suffix if missing)
if [[ -n "${COVERAGE}" ]]; then
    if [[ ! "${COVERAGE}" =~ [Xx]$ ]]; then
        COVERAGE="${COVERAGE}X"
    fi
    # Normalize to uppercase X
    COVERAGE="${COVERAGE%[Xx]}X"
fi

# Normalize TE type case mapping
normalize_te_type() {
    local input="$1"
    case "${input,,}" in
        l1|l1hs)
            echo "L1HS"
            ;;
        aluy|alu)
            echo "AluY"
            ;;
        ltr5|ltr)
            echo "LTR5"
            ;;
        sva_f|sva)
            echo "SVA_F"
            ;;
        *)
            echo "${input}"
            ;;
    esac
}

if [[ -n "${TE_TYPE}" ]]; then
    TE_TYPE=$(normalize_te_type "${TE_TYPE}")
fi

echo "========================================="
echo "Sample List Generation"
echo "========================================="
echo "Data directory: ${DATA_HOME}"
echo "Output directory: ${OUTDIR}"
echo "Batch size: ${BATCH_SIZE}"
[[ -n "${TE_TYPE}" ]] && echo "TE type filter: ${TE_TYPE}" || echo "TE type filter: all"
[[ -n "${COVERAGE}" ]] && echo "Coverage filter: ${COVERAGE}" || echo "Coverage filter: all"
echo "========================================="
echo ""

# Create temporary file for all samples
TEMP_SAMPLE_LIST=$(mktemp)
trap "rm -f ${TEMP_SAMPLE_LIST}" EXIT

# Function to add sample to list
add_sample() {
    local fq1=$1
    local fq2=$2
    local rel_path=$3
    local sample_name=$(basename "$fq1" | sed 's/.1.fq.gz$//')
    
    echo "${sample_name} ${fq1} ${fq2} ${rel_path}" >> "${TEMP_SAMPLE_LIST}"
}

# Determine which TE types to process
if [[ -n "${TE_TYPE}" ]]; then
    TE_TYPES=("${TE_TYPE}")
else
    TE_TYPES=(AluY L1HS LTR5 SVA_F)
fi

# Determine which coverages to process
if [[ -n "${COVERAGE}" ]]; then
    COVERAGES=("${COVERAGE}")
else
    COVERAGES=(5X 10X 50X 100X 200X)
fi

# Process nonReferenceTE samples
echo "Scanning for samples..."
for te_type in "${TE_TYPES[@]}"; do
    for coverage in "${COVERAGES[@]}"; do
        fq_dir="${DATA_HOME}/nonReferenceTE/${te_type}/${coverage}/fq"
        
        if [[ -d "$fq_dir" ]]; then
            # Use LC_ALL=C sort for reproducible ordering
            for fq1 in $(LC_ALL=C find "${fq_dir}" -name "*.1.fq.gz" 2>/dev/null | LC_ALL=C sort); do
                if [[ -f "$fq1" ]]; then
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    if [[ -f "$fq2" ]]; then
                        rel_path="nonReferenceTE/${te_type}/${coverage}"
                        add_sample "$fq1" "$fq2" "$rel_path"
                    fi
                fi
            done
        fi
    done
done

# Process referenceTE samples (if no specific TE type filter)
if [[ -z "${TE_TYPE}" ]]; then
    for coverage in "${COVERAGES[@]}"; do
        # Regular intron samples
        fq_dir="${DATA_HOME}/referenceTE/intron/${coverage}/fq"
        if [[ -d "$fq_dir" ]]; then
            for fq1 in $(LC_ALL=C find "${fq_dir}" -name "reftefu403_*_${coverage}.1.fq.gz" 2>/dev/null | LC_ALL=C sort); do
                if [[ -f "$fq1" ]]; then
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    if [[ -f "$fq2" ]]; then
                        rel_path="referenceTE/intron/${coverage}/fq"
                        add_sample "$fq1" "$fq2" "$rel_path"
                    fi
                fi
            done
        fi
        
        # Mutated intron samples
        fq_mut_dir="${DATA_HOME}/referenceTE/intron/${coverage}/fq_mut"
        if [[ -d "$fq_mut_dir" ]]; then
            for fq1 in $(LC_ALL=C find "${fq_mut_dir}" -name "reftefu403_mut_*_${coverage}.1.fq.gz" 2>/dev/null | LC_ALL=C sort); do
                if [[ -f "$fq1" ]]; then
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    if [[ -f "$fq2" ]]; then
                        rel_path="referenceTE/intron/${coverage}/fq_mut"
                        add_sample "$fq1" "$fq2" "$rel_path"
                    fi
                fi
            done
        fi
        
        # Regular TSS samples
        fq_dir="${DATA_HOME}/referenceTE/TSS/${coverage}/fq"
        if [[ -d "$fq_dir" ]]; then
            for fq1 in $(LC_ALL=C find "${fq_dir}" -name "reftetss270_*_${coverage}.1.fq.gz" 2>/dev/null | LC_ALL=C sort); do
                if [[ -f "$fq1" ]]; then
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    if [[ -f "$fq2" ]]; then
                        rel_path="referenceTE/TSS/${coverage}/fq"
                        add_sample "$fq1" "$fq2" "$rel_path"
                    fi
                fi
            done
        fi
        
        # Mutated TSS samples
        fq_mut_dir="${DATA_HOME}/referenceTE/TSS/${coverage}/fq_mut"
        if [[ -d "$fq_mut_dir" ]]; then
            for fq1 in $(LC_ALL=C find "${fq_mut_dir}" -name "reftetss270_mut_*_${coverage}.1.fq.gz" 2>/dev/null | LC_ALL=C sort); do
                if [[ -f "$fq1" ]]; then
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    if [[ -f "$fq2" ]]; then
                        rel_path="referenceTE/TSS/${coverage}/fq_mut"
                        add_sample "$fq1" "$fq2" "$rel_path"
                    fi
                fi
            done
        fi
    done
fi

# Count total samples
total_samples=$(wc -l < "${TEMP_SAMPLE_LIST}")

if [[ ${total_samples} -eq 0 ]]; then
    echo "ERROR: No samples found matching criteria"
    exit 1
fi

echo "Found ${total_samples} samples"
echo ""

# Calculate number of batch files needed
num_batches=$(( (total_samples + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Creating ${num_batches} batch files..."

# Determine output directory structure
if [[ -n "${TE_TYPE}" && -n "${COVERAGE}" ]]; then
    batch_dir="${OUTDIR}/${TE_TYPE}/${COVERAGE}"
elif [[ -n "${TE_TYPE}" ]]; then
    batch_dir="${OUTDIR}/${TE_TYPE}/all_coverages"
elif [[ -n "${COVERAGE}" ]]; then
    batch_dir="${OUTDIR}/all_te_types/${COVERAGE}"
else
    batch_dir="${OUTDIR}/all_samples"
fi

# Create output directory
mkdir -p "${batch_dir}"

# Split samples into batch files
batch_num=1
line_count=0

# Read through temp file and split
while IFS= read -r line; do
    if [[ $((line_count % BATCH_SIZE)) -eq 0 ]]; then
        # Start new batch file
        batch_file="${batch_dir}/samples_$(printf '%04d' ${batch_num}).list"
        echo "# Sample_Name FQ1 FQ2 Output_Path" > "${batch_file}"
        echo "  Creating ${batch_file}..."
        batch_num=$((batch_num + 1))
    fi
    
    echo "${line}" >> "${batch_file}"
    line_count=$((line_count + 1))
done < "${TEMP_SAMPLE_LIST}"

echo ""
echo "========================================="
echo "Sample List Generation Complete"
echo "========================================="
echo "Total samples: ${total_samples}"
echo "Batch files created: ${num_batches}"
echo "Samples per batch: ${BATCH_SIZE}"
echo "Output directory: ${batch_dir}"
echo "========================================="
echo ""

# Show summary by checking first few batch files
echo "Sample distribution:"
for batch_file in $(LC_ALL=C ls "${batch_dir}"/samples_*.list 2>/dev/null | head -3); do
    sample_count=$(grep -v "^#" "${batch_file}" | wc -l)
    echo "  $(basename "${batch_file}"): ${sample_count} samples"
done

if [[ ${num_batches} -gt 3 ]]; then
    echo "  ..."
    last_batch=$(LC_ALL=C ls "${batch_dir}"/samples_*.list 2>/dev/null | tail -1)
    sample_count=$(grep -v "^#" "${last_batch}" | wc -l)
    echo "  $(basename "${last_batch}"): ${sample_count} samples"
fi

echo ""
echo "Next steps:"
echo "  1. Review batch files in: ${batch_dir}"
echo "  2. Submit jobs: scripts/process_array.sh --sample_list_dir ${batch_dir}"
echo ""

exit 0
