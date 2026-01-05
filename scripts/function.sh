#!/bin/bash

# Common functions for TE analysis pipeline
# This file contains shared functions used by process_samples.sh and process_array.sh

# Function to extract name prefix (everything before first underscore, or full name if no underscore)
extract_name_prefix() {
    local sample_name="$1"
    echo "${sample_name%%_*}"
}

# Function to filter broken FASTQ records
# Checks that sequence and quality lines have the same length
# Returns the path to the filtered file (or original if no filtering needed)
filter_broken_fastq() {
    local fq="$1"

    if [[ ! -f "$fq" ]]; then
        echo "ERROR: FASTQ not found: $fq" >&2
        return 1
    fi

    # Determine base name and check if gzipped - handle multiple extensions
    local base_out=""
    local is_gzipped=false
    
    if [[ "$fq" == *.fq.gz ]]; then
        is_gzipped=true
        base_out="${fq%.fq.gz}"
    elif [[ "$fq" == *.fastq.gz ]]; then
        is_gzipped=true
        base_out="${fq%.fastq.gz}"
    elif [[ "$fq" == *.fq ]]; then
        base_out="${fq%.fq}"
    elif [[ "$fq" == *.fastq ]]; then
        base_out="${fq%.fastq}"
    else
        echo "ERROR: Unrecognized FASTQ file extension: $fq" >&2
        return 1
    fi

    local out="${base_out}.filtered.fq"
    local removed="${base_out}.removed.seqs"
    
    # If original was gzipped, output gzipped filtered file too
    local out_final="$out"
    if [[ "$is_gzipped" == true ]]; then
        out_final="${out}.gz"
    fi

    # Remove old output files if they exist
    [[ -f "$out" ]] && rm "$out"
    [[ -f "$out_final" ]] && rm "$out_final"
    [[ -f "$removed" ]] && rm "$removed"

    # Process the FASTQ file
    if [[ "$is_gzipped" == true ]]; then
        zcat "$fq" | awk '
            NR % 4 == 1 { h = $0 }
            NR % 4 == 2 { s = $0; sl = length($0) }
            NR % 4 == 3 { p = $0 }
            NR % 4 == 0 {
                q = $0; ql = length($0)
                if (sl == ql) {
                    print h >> out
                    print s >> out
                    print p >> out
                    print q >> out
                } else {
                    print h >> bad
                    print s >> bad
                    print p >> bad
                    print q >> bad
                }
            }
        ' out="$out" bad="$removed"
    else
        awk '
            NR % 4 == 1 { h = $0 }
            NR % 4 == 2 { s = $0; sl = length($0) }
            NR % 4 == 3 { p = $0 }
            NR % 4 == 0 {
                q = $0; ql = length($0)
                if (sl == ql) {
                    print h >> out
                    print s >> out
                    print p >> out
                    print q >> out
                } else {
                    print h >> bad
                    print s >> bad
                    print p >> bad
                    print q >> bad
                }
            }
        ' out="$out" bad="$removed" "$fq"
    fi

    # Gzip the filtered output if original was gzipped
    if [[ "$is_gzipped" == true && -f "$out" ]]; then
        gzip "$out"
    fi

    # Check if any sequences were removed
    if [[ -f "$removed" && -s "$removed" ]]; then
        echo "DONE:" >&2
        echo "  valid FASTQ   → $out_final" >&2
        echo "  removed reads → $removed" >&2
        echo "$out_final"  # Return filtered file path
    else
        # No broken sequences found, remove empty files and use original
        [[ -f "$removed" ]] && rm "$removed"
        [[ -f "$out" ]] && rm "$out"
        [[ -f "$out_final" ]] && rm "$out_final"
        echo "No broken FASTQ records found in: $fq" >&2
        echo "$fq"  # Return original file path
    fi
    
    return 0
}

# Function to run JET Step 1
run_jet_step1() {
    echo "[$(date)] Starting JET Step 1 - STAR Alignment..."
    
    outputsDir="${outputDir}/output"
    logDir="${outputDir}/log"
    logFile="${logDir}/step1_multisample_running_$(date +'%Y%m%d').log"
    
    mkdir -p "${outputsDir}"
    
    echo -e "\e[1m${dataDir}\e[0m" > "${logFile}"
    
    # Filter FASTQ files if FQ1 and FQ2 are set
    if [[ -n "${FQ1}" && -f "${FQ1}" ]]; then
        echo "[$(date)] Filtering FQ1: ${FQ1}" | tee -a "${logFile}"
        if ! filtered_fq1=$(filter_broken_fastq "${FQ1}"); then
            echo "ERROR: Failed to filter FQ1: ${FQ1}" | tee -a "${logFile}"
            return 1
        fi
        export FQ1="${filtered_fq1}"
        # rnaSample is used by JET and should match FQ1
        export rnaSample="${filtered_fq1}"
    fi
    
    if [[ -n "${FQ2}" && -f "${FQ2}" ]]; then
        echo "[$(date)] Filtering FQ2: ${FQ2}" | tee -a "${logFile}"
        if ! filtered_fq2=$(filter_broken_fastq "${FQ2}"); then
            echo "ERROR: Failed to filter FQ2: ${FQ2}" | tee -a "${logFile}"
            return 1
        fi
        export FQ2="${filtered_fq2}"
    fi
    
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

  # Filter FASTQ files before processing
  echo "[$(date)] Filtering FQ1: ${FQ1}"
  if ! filtered_fq1=$(filter_broken_fastq "${FQ1}"); then
      echo "ERROR: Failed to filter FQ1: ${FQ1}"
      return 1
  fi
  FQ1="${filtered_fq1}"
  
  echo "[$(date)] Filtering FQ2: ${FQ2}"
  if ! filtered_fq2=$(filter_broken_fastq "${FQ2}"); then
      echo "ERROR: Failed to filter FQ2: ${FQ2}"
      return 1
  fi
  FQ2="${filtered_fq2}"

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
  while read bam_file; do
    bam_basename=$(basename "${bam_file}")
    if [[ ! -e "./${bam_basename}" ]]; then
      ln -s "${bam_file}" "./${bam_basename}"
      bam_count=$((bam_count + 1))
    fi
  done < <(find "${dataDir}" -type f -name "*.bam" 2>/dev/null)

  echo "[$(date)] Linked ${bam_count} BAM files for aggregation"

  # Copy or symlink all *_annotated_filtered_test_all files (output from per-sample processing)
  echo "[$(date)] Collecting annotated GTF files from samples..."
  find "${dataDir}" -type f -name "*_annotated_filtered_test_all" -exec ln -sf {} . \;

  # Execute run_command.sh with proper environment
  echo "[$(date)] Executing run_command.sh for cross-sample aggregation..."

  # NOTE: The original run_command.sh script has hardcoded paths and doesn't accept parameters.
  # The optional configuration variables in config_template.sh are documented for reference,
  # but to use non-default values, users must either:
  # 1. Modify run_command.sh directly, OR
  # 2. Call the R/Python scripts individually with custom parameters
  #
  # Example of calling aggregateProcessedAnnotation.R with custom parameters:
  # ./aggregateProcessedAnnotation.R -a ./arguments.txt -e "treatment" -l 2588 -s 2 -n 1
  
  echo "Using default parameters from run_command.sh"
  echo "To customize, edit run_command.sh or call R/Python scripts directly"

  # Handle the hardcoded gencode path in run_command.sh
  # Line 14 of run_command.sh has: ../genome_46/gencode.v46.basic.annotation.sorted.gtf
  # Create this directory structure to make it work
  echo "[$(date)] Setting up reference path for cuffmerge..."
  
  parent_dir=$(dirname "${TEPROF2_AGGREGATION_DIR}")
  genome_ref_dir="${parent_dir}/genome_46"
  
  if [[ ! -e "${genome_ref_dir}/gencode.v46.basic.annotation.sorted.gtf" ]]; then
    echo "Creating genome_46 directory structure..."
    mkdir -p "${genome_ref_dir}"
    
    if [[ -f "${TEPROF2_CUFFMERGE_GTF}" ]]; then
      ln -sf "${TEPROF2_CUFFMERGE_GTF}" "${genome_ref_dir}/gencode.v46.basic.annotation.sorted.gtf"
      echo "Symlinked ${TEPROF2_CUFFMERGE_GTF} to ${genome_ref_dir}/gencode.v46.basic.annotation.sorted.gtf"
    else
      echo "WARNING: TEPROF2_CUFFMERGE_GTF not found: ${TEPROF2_CUFFMERGE_GTF}"
      echo "You may need to create ${genome_ref_dir}/gencode.v46.basic.annotation.sorted.gtf manually"
    fi
  fi

  # Execute run_command.sh
  echo "[$(date)] Running run_command.sh..."
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
