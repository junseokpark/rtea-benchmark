#!/bin/bash

#SBATCH --job-name=TE_analysis
#SBATCH --output=logs/TE_analysis_%A_%a.out
#SBATCH --error=logs/TE_analysis_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# Configuration
DATA_HOME=/home/junseokp/workspaces/data/rTea-simul
REF=${DATA_HOME}/ref
OUTPUT_BASE=${DATA_HOME}/output

JET=/home/sasidharp/jet_docker/jet.sif
TEProf2=/home/sasidharp/jet_docker/teprof2.sif

THREADS=8

module load singularity

# Create logs directory
mkdir -p logs

# Function to extract sample name from fastq file
get_sample_name() {
    local fq_file=$1
    local basename=$(basename "$fq_file")
    # Remove .1.fq.gz or .2.fq.gz extension
    echo "${basename%.*.fq.gz}"
}

# Function to process with JET
run_jet() {
    local fq1=$1
    local fq2=$2
    local output_dir=$3
    local sample_name=$4
    
    echo "Running JET on ${sample_name}..."
    mkdir -p "${output_dir}/JET"
    
    # JET Step 1: Alignment
    singularity exec ${JET} bwa mem -t ${THREADS} \
        ${REF}/reference.fa \
        ${fq1} ${fq2} | \
        singularity exec ${JET} samtools view -bS - | \
        singularity exec ${JET} samtools sort -@ ${THREADS} -o "${output_dir}/JET/${sample_name}.sorted.bam"
    
    singularity exec ${JET} samtools index "${output_dir}/JET/${sample_name}.sorted.bam"
    
    # JET Step 2: TE Detection
    singularity exec ${JET} python /path/to/JET/scripts/detect_TE.py \
        -b "${output_dir}/JET/${sample_name}.sorted.bam" \
        -r ${REF}/reference.fa \
        -t ${REF}/TE_annotation.bed \
        -o "${output_dir}/JET/${sample_name}" \
        -p ${THREADS}
    
    echo "JET completed for ${sample_name}"
}

# Function to process with TEProf2
run_teprof2() {
    local fq1=$1
    local fq2=$2
    local output_dir=$3
    local sample_name=$4
    
    echo "Running TEProf2 on ${sample_name}..."
    mkdir -p "${output_dir}/TEProf2"
    
    # TEProf2 processing
    singularity exec ${TEProf2} teprof2 \
        --fq1 ${fq1} \
        --fq2 ${fq2} \
        --ref ${REF}/reference.fa \
        --te-annot ${REF}/TE_annotation.gtf \
        --output "${output_dir}/TEProf2/${sample_name}" \
        --threads ${THREADS}
    
    echo "TEProf2 completed for ${sample_name}"
}

# Function to process a sample pair
process_sample() {
    local fq1=$1
    local fq2=$2
    local rel_path=$3  # Relative path from DATA_HOME
    
    # Extract sample name
    local sample_name=$(get_sample_name "$fq1")
    
    # Create output directory maintaining original structure
    local output_dir="${OUTPUT_BASE}/${rel_path}"
    mkdir -p "${output_dir}"
    
    echo "========================================="
    echo "Processing sample: ${sample_name}"
    echo "Input files:"
    echo "  R1: ${fq1}"
    echo "  R2: ${fq2}"
    echo "Output directory: ${output_dir}"
    echo "========================================="
    
    # Run JET
    run_jet "${fq1}" "${fq2}" "${output_dir}" "${sample_name}"
    
    # Run TEProf2
    run_teprof2 "${fq1}" "${fq2}" "${output_dir}" "${sample_name}"
    
    echo "Completed processing ${sample_name}"
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
