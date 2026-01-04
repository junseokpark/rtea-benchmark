#!/bin/bash

# Common functions for TE analysis pipeline
# This file contains shared functions used by process_samples.sh and process_array.sh

# Function to extract name prefix (everything before first underscore, or full name if no underscore)
extract_name_prefix() {
    local sample_name="$1"
    echo "${sample_name%%_*}"
}

# Function to run JET Step 1
run_jet_step1() {
    echo "[$(date)] Starting JET Step 1 - STAR Alignment..."
    
    outputsDir="${outputDir}/output"
    logDir="${outputDir}/log"
    logFile="${logDir}/step1_multisample_running_$(date +'%Y%m%d').log"
    
    mkdir -p "${outputsDir}"
    
    echo -e "\e[1m${dataDir}\e[0m" > "${logFile}"
    
    # Execute JET Step 1 using singularity
    executeCMD="singularity exec --bind \"${JET2_localPath}:${JET2_localPath}\" \
        --bind \"${dataDir}:${dataDir}\" \
        \"${JET2}\" \"${JET2_localPath}/Step1_pipelineJETs_STAR.sh\" \
        --samtools \"${samtoolsBinDir}\" \
        --star \"${starBinDir}\" \
        --read-length \"${readLength}\" \
        --organism \"${organism}\" \
        --genome \"${genome}\" \
        --database \"${database}\" \
        --ref-dir \"${refDir}\" \
        --fasta \"${fastaFile}\" \
        --gtf \"${gtfGeneFile}\" \
        --threads \"${threads}\" \
        --rna-sample \"${rnaSample}\" \
        --name \"${name}\" \
        --name-prefix \"${namePrefix}\" \
        --data-dir \"${dataDir}\" \
        --output \"${outputsDir}\""
    
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
    
    echo -e "\e[1m${dataDir}\e[0m" > "${logFile}"
    
    # Execute JET Step 2 using singularity
    executeCMD="singularity exec \"${JET2}\" \"${JET2_localPath}/Step2_pipelineJETs_R.sh\" \
        --jetprojectdir \"${JET2_localPath}\" \
        --data-dir \"${dataDir}\" \
        --outputs-dir \"${outputsDir}\" \
        --log-dir \"${logDir}\" \
        --star-dir \"${starIndexesDir}\" \
        --rna-sample \"${rnaSample}\" \
        --name \"${name}\" \
        --name-prefix \"${namePrefix}\" \
        --error-dir \"${ErrorDir}\" \
        --read-length \"${readLength}\" \
        --organism \"${organism}\" \
        --genome \"${genome}\" \
        --database \"${database}\" \
        --rlib-dir \"${RlibDir}\" \
        --repeats-file \"${repeatsFile}\" \
        --gff-file \"${gffFile}\""
    
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
  command -v STAR &>/dev/null || { echo "ERROR: STAR not found in PATH"; return 1; }
  command -v samtools &>/dev/null || { echo "ERROR: samtools not found in PATH"; return 1; }
  command -v stringtie &>/dev/null || { echo "ERROR: stringtie not found in PATH"; return 1; }
  command -v rmskhg38_annotate_gtf_update_test_tpm.py &>/dev/null || { echo "ERROR: rmskhg38_annotate_gtf_update_test_tpm.py not found in PATH"; return 1; }
  command -v annotationtpmprocess.py &>/dev/null || { echo "ERROR: annotationtpmprocess.py not found in PATH"; return 1; }
  command -v filterReadCandidates.R &>/dev/null || { echo "ERROR: filterReadCandidates.R not found in PATH"; return 1; }

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

# Function to run TEProf2 aggregation across all samples
# This should be called AFTER all individual samples have been processed
# It executes run_command.sh which performs cross-sample analysis
run_teprof2_aggregation() {
  set -euo pipefail

  echo "[$(date)] Starting TEProf2 aggregation across all samples"

  # ---------- Required variables ----------
  : "${TEPROF2_AGGREGATION_DIR:?Need TEPROF2_AGGREGATION_DIR - directory containing all TEProf2 sample outputs}"
  : "${ARGUMENTS_TXT:?Need ARGUMENTS_TXT - path to arguments.txt file}"
  : "${RUN_COMMAND_SCRIPT:?Need RUN_COMMAND_SCRIPT - path to run_command.sh script}"
  : "${TEPROF2_CUFFMERGE_GTF:?Need TEPROF2_CUFFMERGE_GTF - path to GENCODE GTF for cuffmerge}"

  # Check for required tools (must be in PATH or accessible via container)
  echo "Checking for required tools..."
  required_tools=(
    "aggregateProcessedAnnotation.R"
    "filterReadCandidates.R"
    "mergeAnnotationProcess.R"
    "finalStatisticsOutput.R"
    "rmskhg38_annotate_gtf_update_test_tpm_cuff.py"
    "commandsmax_speed.py"
    "stringtieExpressionFrac.py"
    "samtools"
    "stringtie"
    "gffread"
    "cuffmerge"
  )
  
  missing_tools=()
  for tool in "${required_tools[@]}"; do
    if ! command -v "$tool" &>/dev/null; then
      missing_tools+=("$tool")
    fi
  done
  
  if [ ${#missing_tools[@]} -gt 0 ]; then
    echo "WARNING: The following tools are not in PATH:"
    for tool in "${missing_tools[@]}"; do
      echo "  - $tool"
    done
    echo "These tools must be available in PATH or run this within the TEProf2 container."
    echo "Continuing anyway - run_command.sh will fail if tools are missing."
  fi

  # Check if run_command.sh exists and is executable
  if [[ ! -f "${RUN_COMMAND_SCRIPT}" ]]; then
    echo "ERROR: run_command.sh not found at ${RUN_COMMAND_SCRIPT}"
    return 1
  fi

  if [[ ! -x "${RUN_COMMAND_SCRIPT}" ]]; then
    echo "Making run_command.sh executable..."
    chmod +x "${RUN_COMMAND_SCRIPT}"
  fi

  # Create aggregation directory
  mkdir -p "${TEPROF2_AGGREGATION_DIR}"
  cd "${TEPROF2_AGGREGATION_DIR}"

  echo "[$(date)] Working directory: ${TEPROF2_AGGREGATION_DIR}"

  # Copy arguments.txt to current directory (required by run_command.sh)
  if [[ ! -f "./arguments.txt" ]]; then
    echo "Copying arguments.txt to aggregation directory..."
    cp "${ARGUMENTS_TXT}" ./arguments.txt
  fi

  # Symlink all BAM files from individual sample directories to current directory
  # run_command.sh expects BAM files to be accessible via: find ./ -name "*bam"
  echo "[$(date)] Symlinking BAM files from sample directories..."
  
  # Find all TEProf2 sample output directories containing BAM files
  # Pattern: look for .bam files in subdirectories
  bam_count=0
  find "${dataDir}" -type f -name "*.bam" 2>/dev/null | while read bam_file; do
    bam_basename=$(basename "${bam_file}")
    if [[ ! -e "./${bam_basename}" ]]; then
      ln -s "${bam_file}" "./${bam_basename}"
      bam_count=$((bam_count + 1))
    fi
  done

  echo "[$(date)] Linked BAM files for aggregation"

  # Copy or symlink all *_annotated_filtered_test_all files (output from per-sample processing)
  echo "[$(date)] Collecting annotated GTF files from samples..."
  find "${dataDir}" -type f -name "*_annotated_filtered_test_all" -exec ln -sf {} . \;

  # Execute run_command.sh with proper environment
  echo "[$(date)] Executing run_command.sh for cross-sample aggregation..."

  # Build optional arguments for aggregateProcessedAnnotation.R
  AGG_ARGS="-a ./arguments.txt"
  
  [[ -n "${TEPROF2_AGG_EXT_TREATMENT:-}" ]] && AGG_ARGS+=" -e ${TEPROF2_AGG_EXT_TREATMENT}"
  [[ -n "${TEPROF2_AGG_EXON1_LENGTH_MAX:-}" ]] && AGG_ARGS+=" -l ${TEPROF2_AGG_EXON1_LENGTH_MAX}"
  [[ -n "${TEPROF2_AGG_EXON_SKIP_MAX:-}" ]] && AGG_ARGS+=" -s ${TEPROF2_AGG_EXON_SKIP_MAX}"
  [[ -n "${TEPROF2_AGG_SAMPLE_TOTAL_MIN:-}" ]] && AGG_ARGS+=" -n ${TEPROF2_AGG_SAMPLE_TOTAL_MIN}"
  [[ -n "${TEPROF2_AGG_TREATMENT_TOTAL_MIN:-}" ]] && AGG_ARGS+=" -t ${TEPROF2_AGG_TREATMENT_TOTAL_MIN}"
  [[ -n "${TEPROF2_AGG_TREATMENT_EXCLUSIVE:-}" ]] && AGG_ARGS+=" -x ${TEPROF2_AGG_TREATMENT_EXCLUSIVE}"
  [[ -n "${TEPROF2_AGG_KEEP_NONE:-}" ]] && AGG_ARGS+=" -k ${TEPROF2_AGG_KEEP_NONE}"
  [[ -n "${TEPROF2_AGG_FILTER_FOR_TES:-}" ]] && AGG_ARGS+=" -f ${TEPROF2_AGG_FILTER_FOR_TES}"

  echo "Aggregation arguments: ${AGG_ARGS}"

  # Build optional arguments for filterReadCandidates.R
  FILTER_ARGS=""
  [[ -n "${TEPROF2_FILTER_MIN_READS_IN_TE:-}" ]] && FILTER_ARGS+=" -r ${TEPROF2_FILTER_MIN_READS_IN_TE}"
  [[ -n "${TEPROF2_FILTER_MIN_START_READ:-}" ]] && FILTER_ARGS+=" -s ${TEPROF2_FILTER_MIN_START_READ}"
  [[ -n "${TEPROF2_FILTER_EXONIZATION_MAX_PERCENT:-}" ]] && FILTER_ARGS+=" -e ${TEPROF2_FILTER_EXONIZATION_MAX_PERCENT}"
  [[ -n "${TEPROF2_FILTER_DISTANCE_TE:-}" ]] && FILTER_ARGS+=" -d ${TEPROF2_FILTER_DISTANCE_TE}"

  echo "Filter arguments: ${FILTER_ARGS}"

  # Note: run_command.sh needs to be modified to accept and use these parameters
  # For now, it uses defaults. The script would need to be edited to pass AGG_ARGS and FILTER_ARGS

  # Update the hardcoded gencode path in run_command.sh if needed
  # The script has: ../genome_46/gencode.v46.basic.annotation.sorted.gtf
  # We need to ensure this path exists or modify the script
  gencode_parent_dir=$(dirname "${TEPROF2_CUFFMERGE_GTF}")
  gencode_parent_parent_dir=$(dirname "${gencode_parent_dir}")
  
  # Create a symlink to make the hardcoded path work
  if [[ ! -e "../genome_46" ]]; then
    mkdir -p ../genome_46
    ln -sf "${TEPROF2_CUFFMERGE_GTF}" ../genome_46/gencode.v46.basic.annotation.sorted.gtf
  fi

  # Execute run_command.sh
  bash "${RUN_COMMAND_SCRIPT}"

  if [ $? -ne 0 ]; then
    echo "ERROR: run_command.sh failed during aggregation" | tee -a aggregation.log
    return 1
  fi

  echo "[$(date)] TEProf2 aggregation completed successfully"
  
  # Output files created by run_command.sh:
  echo "Output files:"
  echo "  - filter_combined_candidates.tsv: All TE-gene transcripts"
  echo "  - initial_candidate_list.tsv: Summary of unique candidates"
  echo "  - read_filtered_candidates.tsv: Candidates passing read filters"
  echo "  - candidate_transcripts.gff3: GFF3 format of detected transcripts"
  echo "  - reference_merged_candidates.gtf: Merged reference with candidates"
  echo "  - table_frac_tot_cand: Fraction of expression per candidate"
  echo "  - table_tpm_cand: TPM values per candidate"
  echo "  - Final statistics and annotations"

  return 0
}
