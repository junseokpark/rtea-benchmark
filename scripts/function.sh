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
    
    # Execute JET Step 1 using singularity
    executeCMD="singularity exec ${JET2} /JET/Step1_pipelineJETs_STAR.sh \
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
    
    # Execute JET Step 2 using singularity
    executeCMD="singularity exec ${JET2} /JET/Step2_pipelineJETs_R.sh \
        --jetprojectdir /JET \
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
  set -euo pipefail

  echo "[$(date)] Starting TEProf2 for ${SAMPLE_NAME}"

  # ---------- Required variables ----------
  # Inputs
  : "${FQ1:?Need FQ1}"
  : "${FQ2:?Need FQ2}"
  : "${SAMPLE_NAME:?Need SAMPLE_NAME}"

  # Threads
  threads="${threads:-16}"

  # References
  : "${refDir:?Need refDir}"
  STAR_INDEX="${STAR_INDEX:-${refDir}/STAR_hg38_index}"   # prebuilt STAR index dir
  GENCODE_GTF="${GENCODE_GTF:-${refDir}/gencode.gtf}"     # used by StringTie + cuffmerge
  ARGUMENTS_TXT="${ARGUMENTS_TXT:-${refDir}/arguments.txt}"  # TEProf2 arguments.txt

  # Tools (must be in PATH)
  command -v STAR &>/dev/null
  command -v samtools &>/dev/null
  command -v stringtie &>/dev/null
  command -v gffread &>/dev/null
  command -v cuffmerge &>/dev/null
  command -v rmskhg38_annotate_gtf_update_test_tpm.py &>/dev/null
  command -v annotationtpmprocess.py &>/dev/null
  command -v filterReadCandidates.R &>/dev/null

  # ---------- Outputs ----------
  outRoot="${dataDir:-.}/TEProf2/${SAMPLE_NAME}"
  mkdir -p "${outRoot}"
  cd "${outRoot}"

  # ---------- Step 0: Align FASTQ -> sorted BAM ----------
  echo "[$(date)] Step 0: Alignment (STAR) -> BAM"

  STAR \
    --runThreadN "${threads}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${FQ1}" "${FQ2}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${SAMPLE_NAME}." \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM XS

  BAM="${outRoot}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam"

  samtools index -@ "${threads}" "${BAM}"

  # ---------- Step 1: Assemble transcripts -> sample GTF ----------
  echo "[$(date)] Step 1: Transcript assembly (StringTie) -> GTF"

  SAMPLE_GTF="${outRoot}/${SAMPLE_NAME}.stringtie.gtf"

  stringtie "${BAM}" \
    -p "${threads}" \
    -G "${GENCODE_GTF}" \
    -o "${SAMPLE_GTF}"

  # ---------- Step 2: TEProf2 annotation (normal) ----------
  echo "[$(date)] Step 2: TEProf2 RepeatMasker annotation"

  rmskhg38_annotate_gtf_update_test_tpm.py \
    "${SAMPLE_GTF}" \
    "${ARGUMENTS_TXT}"

  annotated_gtf="${SAMPLE_GTF}_annotated_filtered_test_all"
  [[ -f "${annotated_gtf}" ]] || { echo "ERROR: missing ${annotated_gtf}"; return 1; }

  # ---------- Step 3: TPM processing ----------
  echo "[$(date)] Step 3: TPM processing"
  annotationtpmprocess.py "${annotated_gtf}"

  tpm_processed="${annotated_gtf}_annotation_tpm_processed"
  [[ -f "${tpm_processed}" ]] || { echo "ERROR: missing ${tpm_processed}"; return 1; }

  # ---------- Step 4: Filter read candidates ----------
  echo "[$(date)] Step 4: Filter read candidates"
  filterReadCandidates.R "${tpm_processed}" "${BAM}"

  filtered_output="${tpm_processed}_filtered"
  [[ -f "${filtered_output}" ]] || { echo "ERROR: missing ${filtered_output}"; return 1; }

  echo "[$(date)] TEProf2 completed successfully for ${SAMPLE_NAME}"
  return 0
}
