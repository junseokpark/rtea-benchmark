#!/bin/bash

#SBATCH --job-name=TE_analysis
#SBATCH --output=logs/TE_analysis_%A_%a.out
#SBATCH --error=logs/TE_analysis_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# TE Analysis Pipeline - Sequential Processing
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

# Source shared functions
source "${SCRIPT_DIR}/function.sh"

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
    SAMPLE_NAME=${sample_name}
    FQ1=${fq1}
    FQ2=${fq2}
    
    # Create output directory maintaining original structure
    dataDir="${OUTPUT_BASE}/${rel_path}/${sample_name}"
    mkdir -p "${dataDir}"
    mkdir -p "${dataDir}/log"
    mkdir -p "${dataDir}/err"
    
    echo "========================================="
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
    
    echo "Completed processing ${SAMPLE_NAME}"
    echo ""
}

# Main processing loop
echo "Starting TE analysis pipeline..."
echo "Data home: ${DATA_HOME}"
echo "Output base: ${OUTPUT_BASE}"
echo ""

# Process nonReferenceTE samples
for te_type in AluY L1HS LTR5 SVA_F; do
    for coverage in 5X 10X 50X 100X 200X; do
        fq_dir="${DATA_HOME}/nonReferenceTE/${te_type}/${coverage}/fq"
        
        if [ -d "$fq_dir" ]; then
            echo "Processing ${te_type} ${coverage}..."
            
            # Get all R1 files
            for fq1 in ${fq_dir}/*.1.fq.gz; do
                if [ -f "$fq1" ]; then
                    # Construct R2 filename
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    
                    if [ -f "$fq2" ]; then
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
    if [ -d "$fq_dir" ]; then
        echo "Processing referenceTE/intron ${coverage} (regular)..."
        
        for fq1 in ${fq_dir}/reftefu403_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [ -f "$fq2" ]; then
                    rel_path="referenceTE/intron/${coverage}/fq"
                    process_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
    
    # Process mutated intron samples
    fq_mut_dir="${DATA_HOME}/referenceTE/intron/${coverage}/fq_mut"
    if [ -d "$fq_mut_dir" ]; then
        echo "Processing referenceTE/intron ${coverage} (mutated)..."
        
        for fq1 in ${fq_mut_dir}/reftefu403_mut_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [ -f "$fq2" ]; then
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
    if [ -d "$fq_dir" ]; then
        echo "Processing referenceTE/TSS ${coverage} (regular)..."
        
        for fq1 in ${fq_dir}/reftetss270_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [ -f "$fq2" ]; then
                    rel_path="referenceTE/TSS/${coverage}/fq"
                    process_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
    
    # Process mutated TSS samples
    fq_mut_dir="${DATA_HOME}/referenceTE/TSS/${coverage}/fq_mut"
    if [ -d "$fq_mut_dir" ]; then
        echo "Processing referenceTE/TSS ${coverage} (mutated)..."
        
        for fq1 in ${fq_mut_dir}/reftetss270_mut_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                
                if [ -f "$fq2" ]; then
                    rel_path="referenceTE/TSS/${coverage}/fq_mut"
                    process_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
done

echo "========================================="
echo "All samples processed!"
echo "========================================="
