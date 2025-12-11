#!/bin/bash

#SBATCH --job-name=TE_array
#SBATCH --output=logs/TE_array_%A_%a.out
#SBATCH --error=logs/TE_array_%A_%a.err
#SBATCH --array=1-400%20
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# Configuration
DATA_HOME=/home/junseokp/workspaces/data/rTea-simul
REF=${DATA_HOME}/ref
OUTPUT_BASE=${DATA_HOME}/output

JET=/home/sasidharp/jet_docker/jet.sif
TEProf2=/home/sasidharp/jet_docker/teprof2.sif

THREADS=8
SAMPLE_LIST="sample_list.txt"

module load singularity

# Create logs directory
mkdir -p logs

# Get sample info from array task ID
# Skip header line, get line number matching array task ID
SAMPLE_INFO=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${SAMPLE_LIST})

# Parse sample information
SAMPLE_NAME=$(echo ${SAMPLE_INFO} | awk '{print $1}')
FQ1=$(echo ${SAMPLE_INFO} | awk '{print $2}')
FQ2=$(echo ${SAMPLE_INFO} | awk '{print $3}')
REL_PATH=$(echo ${SAMPLE_INFO} | awk '{print $4}')

# Create output directory
OUTPUT_DIR="${OUTPUT_BASE}/${REL_PATH}"
mkdir -p "${OUTPUT_DIR}/JET" "${OUTPUT_DIR}/TEProf2"

echo "========================================="
echo "Array Job ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing sample: ${SAMPLE_NAME}"
echo "Input files:"
echo "  R1: ${FQ1}"
echo "  R2: ${FQ2}"
echo "Output directory: ${OUTPUT_DIR}"
echo "========================================="

# Function to run JET
run_jet() {
    echo "[$(date)] Starting JET analysis..."
    
    # JET Step 1: Alignment with BWA-MEM
    echo "[$(date)] Step 1: Alignment..."
    singularity exec ${JET} bwa mem -t ${THREADS} \
        ${REF}/reference.fa \
        ${FQ1} ${FQ2} | \
        singularity exec ${JET} samtools view -bS - | \
        singularity exec ${JET} samtools sort -@ ${THREADS} \
        -o "${OUTPUT_DIR}/JET/${SAMPLE_NAME}.sorted.bam"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: BWA alignment failed for ${SAMPLE_NAME}"
        return 1
    fi
    
    # Index BAM file
    echo "[$(date)] Indexing BAM file..."
    singularity exec ${JET} samtools index "${OUTPUT_DIR}/JET/${SAMPLE_NAME}.sorted.bam"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: BAM indexing failed for ${SAMPLE_NAME}"
        return 1
    fi
    
    # JET Step 2: TE Detection
    echo "[$(date)] Step 2: TE Detection..."
    singularity exec ${JET} python /JET/identify_polymorphic_insertions.py \
        -i "${OUTPUT_DIR}/JET/${SAMPLE_NAME}.sorted.bam" \
        -r ${REF}/reference.fa \
        -g ${REF}/gene_annotation.gtf \
        -t ${REF}/TE_annotation.bed \
        -o "${OUTPUT_DIR}/JET/${SAMPLE_NAME}" \
        -p ${THREADS}
    
    if [ $? -ne 0 ]; then
        echo "ERROR: JET TE detection failed for ${SAMPLE_NAME}"
        return 1
    fi
    
    echo "[$(date)] JET analysis completed successfully"
    return 0
}

# Function to run TEProf2
run_teprof2() {
    echo "[$(date)] Starting TEProf2 analysis..."
    
    # TEProf2 processing
    singularity exec ${TEProf2} teprof2 \
        --fq1 ${FQ1} \
        --fq2 ${FQ2} \
        --ref ${REF}/reference.fa \
        --te-annot ${REF}/TE_annotation.gtf \
        --gene-annot ${REF}/gene_annotation.gtf \
        --output-dir "${OUTPUT_DIR}/TEProf2/${SAMPLE_NAME}" \
        --threads ${THREADS} \
        --min-mapq 20 \
        --min-base-quality 20
    
    if [ $? -ne 0 ]; then
        echo "ERROR: TEProf2 analysis failed for ${SAMPLE_NAME}"
        return 1
    fi
    
    echo "[$(date)] TEProf2 analysis completed successfully"
    return 0
}

# Run JET
run_jet
JET_STATUS=$?

# Run TEProf2
run_teprof2
TEPROF2_STATUS=$?

# Summary
echo ""
echo "========================================="
echo "Processing Summary for ${SAMPLE_NAME}"
echo "========================================="
echo "JET Status: $([ ${JET_STATUS} -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "TEProf2 Status: $([ ${TEPROF2_STATUS} -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "========================================="

# Exit with error if either tool failed
if [ ${JET_STATUS} -ne 0 ] || [ ${TEPROF2_STATUS} -ne 0 ]; then
    exit 1
fi

exit 0
