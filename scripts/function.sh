#!/bin/bash

# Common functions for TE analysis pipeline
# This file contains shared functions used by process_samples.sh and process_array.sh

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
