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

    local base_out="" is_gzipped=false
    if [[ "$fq" == *.fq.gz ]]; then
        is_gzipped=true; base_out="${fq%.fq.gz}"
    elif [[ "$fq" == *.fastq.gz ]]; then
        is_gzipped=true; base_out="${fq%.fastq.gz}"
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
    local out_final="$out"
    [[ "$is_gzipped" == true ]] && out_final="${out}.gz"

    rm -f "$out" "$out_final" "$removed"

    # Use gawk (preferred). Works with awk too on most systems.
    local reader="cat"
    [[ "$is_gzipped" == true ]] && reader="zcat"

    $reader "$fq" | awk -v out="$out" -v bad="$removed" '
        function flush_bad(reason) {
            # Write the 4 lines we *think* are a record to bad, plus a comment line.
            # (Comment is safe for inspection; if you want strict FASTQ in removed, delete comment line.)
            print "#REASON=" reason >> bad
            if (h!="") print h >> bad
            if (s!="") print s >> bad
            if (p!="") print p >> bad
            if (q!="") print q >> bad
        }

        BEGIN { state=0; h=s=p=q="" }

        {
            line=$0
            if (state==0) {
                h=line; s=p=q=""
                state=1
            } else if (state==1) {
                s=line
                state=2
            } else if (state==2) {
                p=line
                state=3
            } else if (state==3) {
                q=line
                # Validate structure
                if (substr(h,1,1) != "@") {
                    flush_bad("header_not_at")
                } else if (substr(p,1,1) != "+") {
                    flush_bad("plus_line_missing_or_shifted")
                } else if (length(s) != length(q)) {
                    flush_bad("seq_qual_length_mismatch seqLen=" length(s) " qualLen=" length(q))
                } else {
                    print h >> out
                    print s >> out
                    print p >> out
                    print q >> out
                }
                state=0
            }
        }

        END {
            # If file ends mid-record, dump partial record
            if (state != 0) {
                flush_bad("truncated_record_EOF")
            }
        }
    '

    # gzip output if needed
    if [[ "$is_gzipped" == true && -f "$out" ]]; then
        gzip -f "$out"
    fi

    # Return filtered or original
    if [[ -s "$removed" ]]; then
        echo "DONE:" >&2
        echo "  valid FASTQ   → $out_final" >&2
        echo "  removed reads → $removed" >&2
        echo "$out_final"
    else
        rm -f "$removed" "$out" "$out_final"
        echo "No broken FASTQ records found in: $fq" >&2
        echo "$fq"
    fi
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
    
    # Validate required variables
    if [[ -z "${FQ1}" || -z "${FQ2}" ]]; then
        echo "ERROR: FQ1 and FQ2 must be set" >&2
        return 1
    fi
    
    if [[ ! -f "${FQ1}" ]]; then
        echo "ERROR: FQ1 file not found: ${FQ1}" >&2
        return 1
    fi
    
    if [[ ! -f "${FQ2}" ]]; then
        echo "ERROR: FQ2 file not found: ${FQ2}" >&2
        return 1
    fi
    
    outputsDir="${outputDir}/output"
    logDir="${outputDir}/log"
    logFile="${logDir}/step1_multisample_running_$(date +'%Y%m%d').log"
    
    mkdir -p "${outputsDir}"
    mkdir -p "${logDir}"
    
    echo -e "\e[1m${dataDir}\e[0m" > "${logFile}"
    
    # # Filter FASTQ files
    # echo "[$(date)] Filtering FQ1: ${FQ1}" | tee -a "${logFile}"
    # if ! filtered_fq1=$(filter_broken_fastq "${FQ1}"); then
    #     echo "ERROR: Failed to filter FQ1: ${FQ1}" | tee -a "${logFile}"
    #     return 1
    # fi
    # export FQ1="${filtered_fq1}"
    
    # echo "[$(date)] Filtering FQ2: ${FQ2}" | tee -a "${logFile}"
    # if ! filtered_fq2=$(filter_broken_fastq "${FQ2}"); then
    #     echo "ERROR: Failed to filter FQ2: ${FQ2}" | tee -a "${logFile}"
    #     return 1
    # fi
    # export FQ2="${filtered_fq2}"
    
    # Determine bind paths for singularity
    # Need to bind directories containing FQ1, FQ2, output, and reference files
    local fq1_dir=$(dirname "${FQ1}")
    local fq2_dir=$(dirname "${FQ2}")
    
    # Validate directories exist
    if [[ ! -d "${fq1_dir}" ]]; then
        echo "ERROR: FQ1 directory not found: ${fq1_dir}" | tee -a "${logFile}"
        return 1
    fi
    
    if [[ ! -d "${fq2_dir}" ]]; then
        echo "ERROR: FQ2 directory not found: ${fq2_dir}" | tee -a "${logFile}"
        return 1
    fi
    
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
    
    # Validate required variables
    if [[ -z "${repeatsFile}" ]]; then
        echo "ERROR: repeatsFile must be set" >&2
        return 1
    fi
    
    if [[ -z "${gffFile}" ]]; then
        echo "ERROR: gffFile must be set" >&2
        return 1
    fi
    
    if [[ ! -f "${repeatsFile}" ]]; then
        echo "ERROR: repeatsFile not found: ${repeatsFile}" >&2
        return 1
    fi
    
    if [[ ! -f "${gffFile}" ]]; then
        echo "ERROR: gffFile not found: ${gffFile}" >&2
        return 1
    fi
    
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
    # Bind parent directories of files, not the files themselves
    local repeats_dir=$(dirname "${repeatsFile}")
    local gff_dir=$(dirname "${gffFile}")
    
    # Validate directories exist
    if [[ ! -d "${repeats_dir}" ]]; then
        echo "ERROR: repeatsFile directory not found: ${repeats_dir}" | tee -a "${logFile}"
        return 1
    fi
    
    if [[ ! -d "${gff_dir}" ]]; then
        echo "ERROR: gffFile directory not found: ${gff_dir}" | tee -a "${logFile}"
        return 1
    fi
    
    local bind_paths="${JET2_localPath}:${JET2_localPath},${outputsDir}:${outputsDir},${starIndexesDir}:${starIndexesDir},${repeats_dir}:${repeats_dir}"
    
    # Add gff directory if different from repeats directory
    if [[ "${gff_dir}" != "${repeats_dir}" ]]; then
        bind_paths="${bind_paths},${gff_dir}:${gff_dir}"
    fi
    
    # Add FASTQ directory if step1_fq1 is set and valid
    if [[ -n "${step1_fq1}" ]]; then
        if [[ -f "${step1_fq1}" ]]; then
            local fq1_dir=$(dirname "${step1_fq1}")
            # Only add if different from already bound paths and directory exists
            if [[ -d "${fq1_dir}" && "${fq1_dir}" != "${repeats_dir}" && "${fq1_dir}" != "${gff_dir}" ]]; then
                bind_paths="${bind_paths},${fq1_dir}:${fq1_dir}"
            fi
        else
            echo "WARNING: step1_fq1 file not found: ${step1_fq1}" | tee -a "${logFile}"
        fi
    fi
    
    # Optional env flag for Singularity
    R_LIBS_ENV_OPT=""
    if [[ -n "${JET2_R_LIBS_USER:-}" && -d "${JET2_R_LIBS_USER}" ]]; then
        R_LIBS_ENV_OPT="--env R_LIBS_USER=${JET2_R_LIBS_USER}"
    fi

    executeCMD="singularity exec ${R_LIBS_ENV_OPT} --bind \"${bind_paths}\" \
        \"${JET2}\" \"${JET2_localPath}/Step2_pipelineJETs_R.sh\" \
        --jetprojectdir \"${JET2_localPath}\" \
        --outputs-dir \"${outputsDir}\" \
        --star-dir \"${starIndexesDir}\" \
        --read-length \"${readLength}\" \
        --organism \"${organism}\" \
        --genome \"${genome}\" \
        --database \"${database}\" \
        --data-dir \"${dataDir}\" \
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

# Function to run TEProf2 - Complete TE expression profiling pipeline
#
# Performs TE (transposable element) expression analysis including:
# 1. FASTQ filtering and quality control
# 2. STAR alignment to reference genome
# 3. StringTie transcript assembly
# 4. TEProf2 RepeatMasker annotation
# 5. TPM (transcripts per million) processing
# 6. Read candidate filtering
#
# Required environment variables:
#   FQ1                     - Path to first FASTQ file (read 1), gzipped (.fq.gz or .fastq.gz)
#   FQ2                     - Path to second FASTQ file (read 2), gzipped (.fq.gz or .fastq.gz)
#   SAMPLE_NAME             - Sample name/identifier (used for output naming and error reporting)
#   refDir                  - Base reference directory path (used as default location for reference files)
#   TEProf2                 - Path to TEProf2 singularity container image (e.g., /path/to/teprof2.sif)
#   TEProf2_Local_Path      - Path to TEProf2 scripts directory containing:
#                             * rmskhg38_annotate_gtf_update_test_tpm.py
#                             * annotationtpmprocess.py
#                             * filterReadCandidates.R
#
# Optional environment variables (with defaults):
#   threads                 - Number of CPU threads for parallel processing (default: 16)
#   STAR_INDEX              - Path to prebuilt STAR genome index directory (default: ${refDir}/STAR_hg38_index)
#   GENCODE_GTF             - Path to GENCODE annotation GTF file for StringTie (default: ${refDir}/gencode.gtf)
#   ARGUMENTS_TXT           - Path to TEProf2 arguments.txt configuration file (default: ${refDir}/arguments.txt)
#                             This file must contain paths to RepeatMasker annotations and other TEProf2 references
#   OUTPUT_BASE             - Data output base directory (default: current directory '.')
#                             Output will be written to ${OUTPUT_BASE}/TEProf2/${SAMPLE_NAME}
#   TEProf2_local_STAR_Path - Path to STAR binary within container (optional override, default: uses 'STAR' in PATH)
#   STAR_READ_CMD           - Command to read FASTQ files (default: zcat for gzipped files)
#   STAR_EXTRA_ARGS         - Additional arguments to pass to STAR (default: empty)
#   STRINGTIE_EXTRA_ARGS    - Additional arguments to pass to StringTie (default: empty)
#   RMSK_ANNOTATE_SCRIPT    - Name of RepeatMasker annotation script (default: rmskhg38_annotate_gtf_update_test_tpm.py)
#   TPM_PROCESS_SCRIPT      - Name of TPM processing script (default: annotationtpmprocess.py)
#   FILTER_CANDIDATES_SCRIPT - Name of read candidate filtering script (default: filterReadCandidates.R)
#
# Required files (paths specified via environment variables):
#   ARGUMENTS_TXT           - TEProf2 configuration file with tab-delimited entries:
#                             * rmsk: path to RepeatMasker BED6 file (tabix-indexed)
#                             * rmskannotationfile: TE subfamily/class/family descriptions
#                             * gencodeplusdic: GENCODE plus strand dictionary
#                             * gencodeminusdic: GENCODE minus strand dictionary
#
# Output files (created in ${OUTPUT_BASE}/TEProf2/${SAMPLE_NAME}/):
#   ${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam      - Sorted BAM alignment file
#   ${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam.bai  - BAM index file
#   ${SAMPLE_NAME}.stringtie.gtf                      - Assembled transcripts
#   ${SAMPLE_NAME}.stringtie.gtf_annotated_filtered_test_all        - Annotated with RepeatMasker TEs
#   ${SAMPLE_NAME}.stringtie.gtf_annotated_filtered_test_all_annotation_tpm_processed  - TPM processed
#   ${SAMPLE_NAME}.stringtie.gtf_annotated_filtered_test_all_annotation_tpm_processed_filtered  - Final filtered candidates
#
run_teprof2() {
  set -euo pipefail

  echo "[$(date)] Starting TEProf2 for ${SAMPLE_NAME}"

  # ---------- Required variables ----------
  : "${FQ1:?Need FQ1}"
  : "${FQ2:?Need FQ2}"
  : "${SAMPLE_NAME:?Need SAMPLE_NAME}"
  : "${refDir:?Need refDir}"

  # Threads
  threads="${threads:-16}"

  # References
  STAR_INDEX="${STAR_INDEX:-${refDir}/STAR_hg38_index}"
  GENCODE_GTF="${GENCODE_GTF:-${refDir}/gencode.gtf}"
  ARGUMENTS_TXT="${ARGUMENTS_TXT:-${refDir}/arguments.txt}"

  # TEProf2 Singularity container and paths
  : "${TEProf2:?Need TEProf2 singularity container path}"
  : "${TEProf2_Local_Path:?Need TEProf2_Local_Path for TEProf2 scripts}"

  # Script path parameterization
  STAR_READ_CMD="${STAR_READ_CMD:-zcat}"
  STAR_EXTRA_ARGS="${STAR_EXTRA_ARGS:-}"
  STRINGTIE_EXTRA_ARGS="${STRINGTIE_EXTRA_ARGS:-}"
  RMSK_ANNOTATE_SCRIPT="${RMSK_ANNOTATE_SCRIPT:-rmskhg38_annotate_gtf_update_test_tpm.py}"
  TPM_PROCESS_SCRIPT="${TPM_PROCESS_SCRIPT:-annotationtpmprocess.py}"
  FILTER_CANDIDATES_SCRIPT="${FILTER_CANDIDATES_SCRIPT:-filterReadCandidates.R}"

  # ---------- Outputs ----------
  outRoot="${OUTPUT_BASE:-.}/TEProf2/${SAMPLE_NAME}"
  mkdir -p "${outRoot}" || { echo "ERROR: Failed to create output directory: ${outRoot}"; return 1; }
  cd "${outRoot}" || { echo "ERROR: Failed to change to output directory: ${outRoot}"; return 1; }

  # Determine bind paths for singularity
  local fq1_dir=$(dirname "${FQ1}")
  local fq2_dir=$(dirname "${FQ2}")
  local bind_paths="${TEProf2_Local_Path},${outRoot},${refDir},${fq1_dir}"
  
  [[ "${fq2_dir}" != "${fq1_dir}" ]] && bind_paths="${bind_paths},${fq2_dir}"
  
  # Add additional directories if they are different from refDir and not already in bind_paths
  for path in "${STAR_INDEX}" "${GENCODE_GTF}" "${ARGUMENTS_TXT}"; do
    local dir=$(dirname "$path")
    if [[ "$dir" != "${refDir}" ]] && [[ ! "$bind_paths" =~ (^|,)"$dir"(,|$) ]]; then
      bind_paths="${bind_paths},$dir"
    fi
  done

  # Determine STAR binary path
  local star_cmd="STAR"
  [[ -n "${TEProf2_local_STAR_Path:-}" ]] && star_cmd="${TEProf2_local_STAR_Path}/STAR"

  # ---------- Step 0: Align FASTQ -> sorted BAM ----------
  echo "[$(date)] [Step 0] STAR alignment"

  BAM="${outRoot}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam"

  singularity exec --bind "${bind_paths}" "${TEProf2}" bash -c '
    source activate teprof2 && \
    '"${star_cmd}"' \
      --runThreadN '"${threads}"' \
      --genomeDir "'"${STAR_INDEX}"'" \
      --readFilesIn "'"${FQ1}"'" "'"${FQ2}"'" \
      --readFilesCommand '"${STAR_READ_CMD}"' \
      --outFileNamePrefix "'"${SAMPLE_NAME}"'." \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM XS \
      '"${STAR_EXTRA_ARGS}"'
  '

  [[ -f "${BAM}" ]] || { echo "ERROR: STAR alignment failed - BAM not created: ${BAM}"; return 1; }

  echo "[$(date)] [Step 0.5] BAM indexing"
  singularity exec --bind "${bind_paths}" "${TEProf2}" bash -c '
    source activate teprof2 && \
    samtools index -@ '"${threads}"' "'"${BAM}"'"
  '

  # ---------- Step 1: Assemble transcripts -> sample GTF ----------
  echo "[$(date)] [Step 1] Transcript assembly (StringTie) -> GTF"

  SAMPLE_GTF="${outRoot}/${SAMPLE_NAME}.stringtie.gtf"

  singularity exec --bind "${bind_paths}" "${TEProf2}" bash -c '
    source activate teprof2 && \
    stringtie "'"${BAM}"'" \
      -p '"${threads}"' \
      -G "'"${GENCODE_GTF}"'" \
      -o "'"${SAMPLE_GTF}"'" \
      '"${STRINGTIE_EXTRA_ARGS}"'
  '

  [[ -f "${SAMPLE_GTF}" ]] || { echo "ERROR: StringTie assembly failed - GTF not created: ${SAMPLE_GTF}"; return 1; }

  # ---------- Step 2: TEProf2 annotation ----------
  echo "[$(date)] [Step 2] TEProf2 RepeatMasker annotation"

  singularity exec --bind "${bind_paths}" "${TEProf2}" bash -c '
    source activate teprof2 && \
    '"${TEProf2_Local_Path}"'/'"${RMSK_ANNOTATE_SCRIPT}"' \
      "'"${SAMPLE_GTF}"'" \
      "'"${ARGUMENTS_TXT}"'"
  '

  annotated_gtf="${SAMPLE_GTF}_annotated_filtered_test_all"
  [[ -f "${annotated_gtf}" ]] || { echo "ERROR: missing ${annotated_gtf}"; return 1; }

  # ---------- Step 3: TPM processing ----------
  echo "[$(date)] [Step 3] TPM processing"
  
  singularity exec --bind "${bind_paths}" "${TEProf2}" bash -c '
    source activate teprof2 && \
    '"${TEProf2_Local_Path}"'/'"${TPM_PROCESS_SCRIPT}"' "'"${annotated_gtf}"'"
  '

  tpm_processed="${annotated_gtf}_annotation_tpm_processed"
  [[ -f "${tpm_processed}" ]] || { echo "ERROR: missing ${tpm_processed}"; return 1; }

  # ---------- Step 4: Filter read candidates ----------
  echo "[$(date)] [Step 4] Filter read candidates"
  
  singularity exec --bind "${bind_paths}" "${TEProf2}" bash -c '
    source activate teprof2 && \
    '"${TEProf2_Local_Path}"'/'"${FILTER_CANDIDATES_SCRIPT}"' "'"${tpm_processed}"'" "'"${BAM}"'"
  '

  filtered_output="${tpm_processed}_filtered"
  [[ -f "${filtered_output}" ]] || { echo "ERROR: missing ${filtered_output}"; return 1; }

  echo "[$(date)] TEProf2 completed successfully for ${SAMPLE_NAME}"
  return 0
}

# Function to run TEProf2 aggregation across all samples
#
# Aggregates results from individual TEProf2 sample analyses to identify
# TE-derived transcripts across the entire dataset. This function MUST be
# called AFTER all individual samples have been processed with run_teprof2().
#
# The aggregation pipeline performs the following steps:
# 1. Aggregates processed annotations from all samples
# 2. Filters candidates based on read support across samples
# 3. Merges candidate transcripts with reference annotations using cuffmerge
# 4. Quantifies expression levels (TPM and fraction) across all samples
# 5. Generates final statistics and candidate lists
#
# Required environment variables:
#   TEPROF2_AGGREGATION_DIR  - Output directory for aggregation results
#                              All aggregated output files will be written here
#   TEPROF2_ARGUMENTS_FILE   - Path to TEProf2 arguments.txt configuration file
#                              Same file used in individual sample processing
#                              Must contain tab-delimited entries:
#                              * rmsk: RepeatMasker BED6 file (tabix-indexed)
#                              * rmskannotationfile: TE descriptions (subfamily, class, family)
#                              * gencodeplusdic: GENCODE plus strand dictionary
#                              * gencodeminusdic: GENCODE minus strand dictionary
#   RUN_COMMAND_SCRIPT       - Path to run_command.sh script (included in repository)
#                              This script orchestrates the aggregation pipeline
#                              Location: scripts/run_command.sh
#   TEPROF2_CUFFMERGE_GTF    - Path to GENCODE GTF file for cuffmerge reference merging
#                              Used by cuffmerge to merge candidate transcripts with reference
#                              (e.g., gencode.v46.basic.annotation.sorted.gtf)
#   dataDir                  - Base directory containing all processed TEProf2 sample outputs
#                              Function will search recursively for:
#                              * *.bam files (alignment files from each sample)
#                              * *_annotated_filtered_test_all files (processed annotations)
#
# Optional environment variables:
#   TEProf2_Local_Path       - Path to TEProf2 scripts directory (optional)
#                              If set, tools will be located at ${TEProf2_Local_Path}/<tool>
#                              If not set, tools must be in system PATH
#                              Required tools/scripts:
#                              * aggregateProcessedAnnotation.R
#                              * filterReadCandidates.R
#                              * mergeAnnotationProcess.R
#                              * finalStatisticsOutput.R
#                              * rmskhg38_annotate_gtf_update_test_tpm_cuff.py
#                              * commandsmax_speed.py
#                              * stringtieExpressionFrac.py
#
# Required tools (must be in PATH or accessible via TEProf2 container):
#   System tools:
#     * samtools              - BAM file manipulation
#     * stringtie             - Transcript assembly and quantification
#     * gffread               - GFF/GTF format conversion
#     * cuffmerge             - Merging transcript assemblies with reference
#
# Configuration note:
#   The aggregation parameters are HARDCODED in run_command.sh and cannot be
#   configured via command-line arguments. To customize filtering parameters
#   (e.g., min reads, exon skip max, treatment filtering), you must either:
#   1. Edit run_command.sh directly to pass custom parameters to R/Python scripts, OR
#   2. Run the R/Python scripts individually with custom parameters
#
#   For reference, common configurable parameters (see config_template.sh):
#     * TEPROF2_AGG_EXON1_LENGTH_MAX (default: 2588)
#     * TEPROF2_AGG_EXON_SKIP_MAX (default: 2)
#     * TEPROF2_AGG_SAMPLE_TOTAL_MIN (default: 1)
#     * TEPROF2_FILTER_MIN_READS_IN_TE (default: 10)
#     * TEPROF2_FILTER_MIN_START_READ (default: 1)
#     * TEPROF2_FILTER_DISTANCE_TE (default: 2500)
#
# Output files (created in TEPROF2_AGGREGATION_DIR):
#   filter_combined_candidates.tsv       - All TE-gene transcript candidates
#   initial_candidate_list.tsv           - Summary of unique candidates
#   read_filtered_candidates.tsv         - Candidates passing read support filters
#   candidate_transcripts.gff3           - GFF3 format of detected transcripts
#   candidate_transcripts.gtf            - GTF format of detected transcripts
#   reference_merged_candidates.gtf      - Candidates merged with reference annotations
#   reference_merged_candidates.gff3     - GFF3 version of merged candidates
#   table_frac_tot_cand                  - Fraction of expression per candidate across samples
#   table_tpm_cand                       - TPM values per candidate across samples
#   table_i_all                          - Intron coverage table
#   Final_output_unique_novel_start_annotated  - Final statistics and annotations
#
run_teprof2_aggregation() {
  set -euo pipefail

  echo "[$(date)] Starting TEProf2 aggregation across all samples"

  # ---------- Required variables ----------
  : "${TEPROF2_AGGREGATION_DIR:?Need TEPROF2_AGGREGATION_DIR - directory containing all TEProf2 sample outputs}"
  : "${TEPROF2_ARGUMENTS_FILE:?Need TEPROF2_ARGUMENTS_FILE - path to arguments.txt file}"
  : "${RUN_COMMAND_SCRIPT:?Need RUN_COMMAND_SCRIPT - path to run_command.sh script}"
  : "${TEPROF2_CUFFMERGE_GTF:?Need TEPROF2_CUFFMERGE_GTF - path to GENCODE GTF for cuffmerge}"

  # Check for required tools (must be in PATH or accessible via container)
  echo "Checking for required tools..."
  
  # TEProf2-specific scripts that should use TEProf2_Local_Path if set
  teprof2_scripts=(
    "aggregateProcessedAnnotation.R"
    "filterReadCandidates.R"
    "mergeAnnotationProcess.R"
    "finalStatisticsOutput.R"
    "rmskhg38_annotate_gtf_update_test_tpm_cuff.py"
    "commandsmax_speed.py"
    "stringtieExpressionFrac.py"
  )
  
  # System tools that should be in PATH
  system_tools=(
    "samtools"
    "stringtie"
    "gffread"
    "cuffmerge"
  )
  
  missing_tools=()
  
  # Check TEProf2 scripts with conditional path prefix
  for tool in "${teprof2_scripts[@]}"; do
    if [[ -n "${TEProf2_Local_Path}" ]]; then
      # When TEProf2_Local_Path is set, check if tool exists and is executable there
      if [[ ! -x "${TEProf2_Local_Path}/${tool}" ]]; then
        missing_tools+=("$tool")
      fi
    else
      # When TEProf2_Local_Path is not set, check if tool is in PATH
      if ! command -v "$tool" &>/dev/null; then
        missing_tools+=("$tool")
      fi
    fi
  done
  
  # Check system tools in PATH
  for tool in "${system_tools[@]}"; do
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
    cp "${TEPROF2_ARGUMENTS_FILE}" ./arguments.txt
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
