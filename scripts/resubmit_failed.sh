#!/bin/bash

# Resubmit failed or incomplete samples
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

FAILED_LIST="failed_samples_$(date +%Y%m%d_%H%M%S).txt"

echo "# Failed/Incomplete Samples - Generated $(date)" > ${FAILED_LIST}
echo "# Format: Sample_Name FQ1 FQ2 Output_Path Tool_Status" >> ${FAILED_LIST}

# Function to check if JET completed successfully
check_jet() {
    local output_dir=$1
    local sample_name=$2
    
    if [ -f "${output_dir}/JET/${sample_name}.sorted.bam" ] && \
       [ -f "${output_dir}/JET/${sample_name}.sorted.bam.bai" ]; then
        return 0  # Success
    else
        return 1  # Failed or incomplete
    fi
}

# Function to check if TEProf2 completed successfully
check_teprof2() {
    local output_dir=$1
    local sample_name=$2
    
    if [ -d "${output_dir}/TEProf2/${sample_name}" ] && \
       [ "$(ls -A ${output_dir}/TEProf2/${sample_name} 2>/dev/null)" ]; then
        return 0  # Success
    else
        return 1  # Failed or incomplete
    fi
}

echo "Identifying failed/incomplete samples..."

failed_count=0
while IFS=' ' read -r sample_name fq1 fq2 rel_path; do
    # Skip header
    if [[ "$sample_name" == "#"* ]]; then
        continue
    fi
    
    output_dir="${OUTPUT_BASE}/${rel_path}"
    
    # Check both tools
    check_jet "$output_dir" "$sample_name"
    jet_ok=$?
    
    check_teprof2 "$output_dir" "$sample_name"
    teprof2_ok=$?
    
    # If either failed, add to failed list
    if [ $jet_ok -ne 0 ] || [ $teprof2_ok -ne 0 ]; then
        status="JET:"
        [ $jet_ok -eq 0 ] && status="${status}OK" || status="${status}FAILED"
        status="${status},TEProf2:"
        [ $teprof2_ok -eq 0 ] && status="${status}OK" || status="${status}FAILED"
        
        echo "${sample_name} ${fq1} ${fq2} ${rel_path} ${status}" >> ${FAILED_LIST}
        failed_count=$((failed_count + 1))
    fi
    
done < ${SAMPLE_LIST}

echo "Found ${failed_count} failed/incomplete samples"
echo "Failed list saved to: ${FAILED_LIST}"

if [ $failed_count -eq 0 ]; then
    echo "All samples completed successfully!"
    rm ${FAILED_LIST}
    exit 0
fi

# Create resubmission script
RESUBMIT_SCRIPT="resubmit_failed_$(date +%Y%m%d_%H%M%S).sh"

cat > ${RESUBMIT_SCRIPT} << 'EOF'
#!/bin/bash

#SBATCH --job-name=TE_resubmit
#SBATCH --output=logs/TE_resubmit_%A_%a.out
#SBATCH --error=logs/TE_resubmit_%A_%a.err
#SBATCH --array=1-ARRAY_SIZE%10
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# Resubmission script for failed samples
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

FAILED_LIST="FAILED_LIST_FILE"

module load singularity

mkdir -p logs

# Get sample info from array task ID
SAMPLE_INFO=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${FAILED_LIST})

# Parse sample information
SAMPLE_NAME=$(echo ${SAMPLE_INFO} | awk '{print $1}')
FQ1=$(echo ${SAMPLE_INFO} | awk '{print $2}')
FQ2=$(echo ${SAMPLE_INFO} | awk '{print $3}')
REL_PATH=$(echo ${SAMPLE_INFO} | awk '{print $4}')
TOOL_STATUS=$(echo ${SAMPLE_INFO} | awk '{print $5}')

OUTPUT_DIR="${OUTPUT_BASE}/${REL_PATH}"
mkdir -p "${OUTPUT_DIR}/JET" "${OUTPUT_DIR}/TEProf2"

echo "========================================="
echo "Resubmitting: ${SAMPLE_NAME}"
echo "Status: ${TOOL_STATUS}"
echo "========================================="

# Check if JET needs to be rerun
if [[ $TOOL_STATUS == *"JET:FAILED"* ]]; then
    echo "[$(date)] Running JET..."
    
    # Clean up previous failed attempt
    rm -rf "${OUTPUT_DIR}/JET/${SAMPLE_NAME}"*
    
    # Run JET
    singularity exec ${JET} bwa mem -t ${threads} \
        ${REF_DIR}/reference.fa \
        ${FQ1} ${FQ2} | \
        singularity exec ${JET} samtools view -bS - | \
        singularity exec ${JET} samtools sort -@ ${threads} \
        -o "${OUTPUT_DIR}/JET/${SAMPLE_NAME}.sorted.bam"
    
    singularity exec ${JET} samtools index "${OUTPUT_DIR}/JET/${SAMPLE_NAME}.sorted.bam"
    
    singularity exec ${JET} python /JET/identify_polymorphic_insertions.py \
        -i "${OUTPUT_DIR}/JET/${SAMPLE_NAME}.sorted.bam" \
        -r ${REF_DIR}/reference.fa \
        -g ${REF_DIR}/gene_annotation.gtf \
        -t ${REF_DIR}/TE_annotation.bed \
        -o "${OUTPUT_DIR}/JET/${SAMPLE_NAME}" \
        -p ${threads}
fi

# Check if TEProf2 needs to be rerun
if [[ $TOOL_STATUS == *"TEProf2:FAILED"* ]]; then
    echo "[$(date)] Running TEProf2..."
    
    # Clean up previous failed attempt
    rm -rf "${OUTPUT_DIR}/TEProf2/${SAMPLE_NAME}"
    
    # Run TEProf2
    singularity exec ${TEProf2} teprof2 \
        --fq1 ${FQ1} \
        --fq2 ${FQ2} \
        --ref ${REF_DIR}/reference.fa \
        --te-annot ${REF_DIR}/TE_annotation.gtf \
        --gene-annot ${REF_DIR}/gene_annotation.gtf \
        --output-dir "${OUTPUT_DIR}/TEProf2/${SAMPLE_NAME}" \
        --threads ${threads} \
        --min-mapq 20 \
        --min-base-quality 20
fi

echo "[$(date)] Resubmission completed for ${SAMPLE_NAME}"
EOF

# Update placeholders in resubmission script
sed -i "s/ARRAY_SIZE/$failed_count/" ${RESUBMIT_SCRIPT}
sed -i "s|FAILED_LIST_FILE|${FAILED_LIST}|" ${RESUBMIT_SCRIPT}

chmod +x ${RESUBMIT_SCRIPT}

echo ""
echo "Resubmission script created: ${RESUBMIT_SCRIPT}"
echo ""
echo "To resubmit failed samples, run:"
echo "  sbatch ${RESUBMIT_SCRIPT}"
echo ""
echo "Or to run a specific sample manually:"
echo "  bash ${RESUBMIT_SCRIPT}"
echo "  (after commenting out SBATCH directives)"
