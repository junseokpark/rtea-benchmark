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

# Function to run JET Step 1 - STAR Alignment
# 
# Updated to pass FASTQ files directly as arguments per revised JET pipeline
# 
# Required environment variables:
#   FQ1             - Path to first FASTQ file (read 1) 
#   FQ2             - Path to second FASTQ file (read 2)
#   outputDir       - Base output directory
#   dataDir         - Data directory
#   JET2            - Path to JET singularity image
#   JET2_localPath  - Path to JET pipeline scripts
#   samtoolsBinDir  - Path to samtools binary (within container)
#   starBinDir      - Path to STAR binary (within container)
#   readLength      - Read length (e.g., 100, 150)
#   organism        - Organism name (e.g., Human, Mouse)
#   genome          - Genome version (e.g., hg38, mm10)
#   database        - Database name (e.g., ensembl)
#   refDir          - Reference directory path
#   fastaFile       - Reference FASTA filename (in refDir)
#   gtfGeneFile     - Gene annotation GTF filename (in refDir)
#   threads         - Number of CPU threads
#   SAMPLE_NAME     - Sample name (for error reporting)
#
run_jet_step1() {
    echo "[$(date)] Starting JET Step 1 - STAR Alignment..."
    
    outputsDir="${outputDir}/output"
    logDir="${outputDir}/log"
    logFile="${logDir}/step1_multisample_running_$(date +'%Y%m%d').log"
    
    mkdir -p "${outputsDir}"
    mkdir -p "${logDir}"
    
    echo -e "\e[1m${dataDir}\e[0m" > "${logFile}"
    
    # Filter FASTQ files if FQ1 and FQ2 are set
    if [[ -n "${FQ1}" && -f "${FQ1}" ]]; then
        echo "[$(date)] Filtering FQ1: ${FQ1}" | tee -a "${logFile}"
        if ! filtered_fq1=$(filter_broken_fastq "${FQ1}"); then
            echo "ERROR: Failed to filter FQ1: ${FQ1}" | tee -a "${logFile}"
            return 1
        fi
        export FQ1="${filtered_fq1}"
    fi
    
    if [[ -n "${FQ2}" && -f "${FQ2}" ]]; then
        echo "[$(date)] Filtering FQ2: ${FQ2}" | tee -a "${logFile}"
        if ! filtered_fq2=$(filter_broken_fastq "${FQ2}"); then
            echo "ERROR: Failed to filter FQ2: ${FQ2}" | tee -a "${logFile}"
            return 1
        fi
        export FQ2="${filtered_fq2}"
    fi
    
    # Determine bind paths for singularity
    # Need to bind directories containing FQ1, FQ2, output, and reference files
    local fq1_dir=$(dirname "${FQ1}")
    local fq2_dir=$(dirname "${FQ2}")
    local bind_paths="${JET2_localPath}:${JET2_localPath},${outputsDir}:${outputsDir},${refDir}:${refDir},${fq1_dir}:${fq1_dir}"
    if [[ "${fq2_dir}" != "${fq1_dir}" ]]; then
        bind_paths="${bind_paths},${fq2_dir}:${fq2_dir}"
    fi
    
    # Execute JET Step 1 using singularity with revised arguments
    # Note: FASTQ files are now passed directly via --fq1 and --fq2
    executeCMD="singularity exec --bind \"${bind_paths}\" \
        \"${JET2}\" \"${JET2_localPath}/Step1_pipelineJETs_STAR.sh\" \
        --samtools \"${samtoolsBinDir}\" \
        --star \"${starBinDir}\" \
        --read-length \"${readLength}\" \
        --organism \"${organism}\" \
        --genome \"${genome}\" \
        --database \"${database}\" \
        --data-dir \"${dataDir}\" \
        --ref-dir \"${refDir}\" \
        --fasta \"${fastaFile}\" \
        --gtf \"${gtfGeneFile}\" \
        --fq1 \"${FQ1}\" \
        --fq2 \"${FQ2}\" \
        --threads \"${threads}\" \
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

# Function to run JET Step 2 - R Analysis
#
# Updated with revised arguments per new JET pipeline specification
#
# Required environment variables:
#   outputDir        - Base output directory
#   FQ1              - Path to FASTQ file from Step 1 (used as --step1-fq1)
#   JET2             - Path to JET singularity image
#   JET2_localPath   - Path to JET pipeline scripts
#   starIndexesDir   - Path to STAR indexes directory
#   readLength       - Read length (default: 100)
#   organism         - Organism name (e.g., Human, Mouse)
#   genome           - Genome version (e.g., hg38, mm10)
#   database         - Database name (e.g., ensembl)
#   RlibDir          - Path to R library directory (within container)
#   repeatsFile      - Path to repeats file
#   gffFile          - Path to GFF annotation file
#   minJunction      - Minimum junction size (optional, default: 2e7)
#   SAMPLE_NAME      - Sample name (for error reporting)
#
run_jet_step2() {
    echo "[$(date)] Starting JET Step 2 - R Analysis..."
    
    outputsDir="${outputDir}/output"
    logDir="${outputDir}/log"
    logFile="${logDir}/step2_multisample_running_$(date +'%Y%m%d').log"
    
    mkdir -p "${logDir}"
    
    echo -e "\e[1m${outputsDir}\e[0m" > "${logFile}"
    
    # Determine step1-fq1 path - use filtered FQ1 if available, otherwise use original
    local step1_fq1="${FQ1}"
    
    # Set default minJunction if not provided
    local minJunction="${minJunction:-2e7}"
    
    # Determine bind paths for singularity
    local bind_paths="${JET2_localPath}:${JET2_localPath},${outputsDir}:${outputsDir},${starIndexesDir}:${starIndexesDir},${repeatsFile}:${repeatsFile},${gffFile}:${gffFile}"
    if [[ -n "${step1_fq1}" ]]; then
        local fq1_dir=$(dirname "${step1_fq1}")
        bind_paths="${bind_paths},${fq1_dir}:${fq1_dir}"
    fi
    
    # Execute JET Step 2 using singularity with revised arguments
    executeCMD="singularity exec --bind \"${bind_paths}\" \
        \"${JET2}\" \"${JET2_localPath}/Step2_pipelineJETs_R.sh\" \
        --jetprojectdir \"${JET2_localPath}\" \
        --outputs-dir \"${outputsDir}\" \
        --star-dir \"${starIndexesDir}\" \
        --read-length \"${readLength}\" \
        --organism \"${organism}\" \
        --genome \"${genome}\" \
        --database \"${database}\" \
        --rlib-dir \"${RlibDir}\" \
        --repeats-file \"${repeatsFile}\" \
        --gff-file \"${gffFile}\" \
        --min-junction \"${minJunction}\" \
        --step1-fq1 \"${step1_fq1}\""
    
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
