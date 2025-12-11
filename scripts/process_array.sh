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
OUTPUT_BASE=${DATA_HOME}/output

# JET Configuration
JETProjectDir=/path/to/JET  # UPDATE THIS PATH
samtoolsBinDir=/path/to/samtools/bin  # UPDATE THIS PATH
starBinDir=/path/to/STAR/bin  # UPDATE THIS PATH
readLength=150  # UPDATE if different
organism="human"  # UPDATE if different
genome="hg38"  # UPDATE if different
database="database_name"  # UPDATE THIS
refDir=${DATA_HOME}/ref
fastaFile=${refDir}/reference.fa
gtfGeneFile=${refDir}/gene_annotation.gtf
starIndexesDir=${refDir}/STAR_indexes
repeatsFile=${refDir}/repeats.txt
gffFile=${refDir}/TE_annotation.gff
RlibDir=/path/to/R/library  # UPDATE THIS PATH
threads=8

# TEProf2 Configuration
TEProf2=/home/sasidharp/jet_docker/teprof2.sif

SAMPLE_LIST="sample_list.txt"

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

# Function to run JET Step 1
run_jet_step1() {
    echo "[$(date)] Starting JET Step 1 - STAR Alignment..."
    
    outputsDir="${dataDir}/output"
    logDir="${dataDir}/log"
    logFile="${logDir}/step1_multisample_running_$(date +'%Y%m%d').log"
    
    mkdir -p ${outputsDir}
    
    echo -e "\e[1m${dataDir}\t${metaFile}\e[0m" > "${logFile}"
    
    # Execute JET Step 1
    executeCMD="${JETProjectDir}/Step1_pipelineJETs_STAR.sh \
        --samtools ${samtoolsBinDir} \
        --star ${starBinDir} \
        --read-length ${readLength} \
        --organism ${organism} \
        --genome ${genome} \
        --database ${database} \
        --ref-dir ${refDir} \
        --fasta ${fastaFile} \
        --gtf ${gtfGeneFile} \
        --threads ${threads} \
        --meta ${metaFile} \
        --data-dir ${dataDir} \
        --output ${outputsDir}"
    
    echo $executeCMD >> "${logFile}"
    eval $executeCMD
    
    if [ $? -ne 0 ]; then
        echo "ERROR: JET Step 1 failed for ${SAMPLE_NAME}" | tee -a "${logFile}"
        return 1
    fi
    
    echo "[$(date)] JET Step 1 completed successfully" | tee -a "${logFile}"
    return 0
}

# Function to run JET Step 2
run_jet_step2() {
    echo "[$(date)] Starting JET Step 2 - R Analysis..."
    
    outputsDir="${dataDir}/output"
    logDir="${dataDir}/log"
    ErrorDir="${dataDir}/err"
    logFile="${logDir}/step2_multisample_running_$(date +'%Y%m%d').log"
    
    echo -e "\e[1m${dataDir}\t${metaFile}\e[0m" > "${logFile}"
    
    # Execute JET Step 2
    executeCMD="${JETProjectDir}/Step2_pipelineJETs_R.sh \
        --jetprojectdir ${JETProjectDir} \
        --data-dir ${dataDir} \
        --outputs-dir ${outputsDir} \
        --log-dir ${logDir} \
        --star-dir ${starIndexesDir} \
        --metadata ${metaFile} \
        --error-dir ${ErrorDir} \
        --read-length ${readLength} \
        --organism ${organism} \
        --genome ${genome} \
        --database ${database} \
        --rlib-dir ${RlibDir} \
        --repeats-file ${repeatsFile} \
        --gff-file ${gffFile}"
    
    echo $executeCMD >> "${logFile}"
    eval $executeCMD
    
    if [ $? -ne 0 ]; then
        echo "ERROR: JET Step 2 failed for ${SAMPLE_NAME}" | tee -a "${logFile}"
        return 1
    fi
    
    echo "[$(date)] JET Step 2 completed successfully" | tee -a "${logFile}"
    return 0
}

# Function to run TEProf2
run_teprof2() {
    echo "[$(date)] Starting TEProf2 analysis..."
    
    teprof2_output="${dataDir}/TEProf2"
    mkdir -p "${teprof2_output}"
    
    # TODO: Update with actual TEProf2 command structure
    # This is a placeholder - please provide the actual TEProf2 command
    singularity exec ${TEProf2} teprof2 \
        --fq1 ${FQ1} \
        --fq2 ${FQ2} \
        --ref ${refDir}/reference.fa \
        --te-annot ${refDir}/TE_annotation.gtf \
        --gene-annot ${refDir}/gene_annotation.gtf \
        --output-dir "${teprof2_output}" \
        --threads ${threads} \
        --min-mapq 20 \
        --min-base-quality 20
    
    if [ $? -ne 0 ]; then
        echo "ERROR: TEProf2 analysis failed for ${SAMPLE_NAME}"
        return 1
    fi
    
    echo "[$(date)] TEProf2 analysis completed successfully"
    return 0
}

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
