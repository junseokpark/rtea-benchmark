#!/bin/bash
# JET and TEProf2 Pipeline Configuration
# Copy this file to config.sh and fill in your actual paths

# ============================================
# REQUIRED: UPDATE ALL PATHS BELOW
# ============================================

# Data Directories
export DATA_HOME="/home/junseokp/workspaces/data/rTea-simul/sims"
export OUTPUT_BASE="/home/junseokp/workspaces/projects/rtea/results"
export REF_DIR="/home/junseokp/workspaces/data/rTea-simul/ref"



# ============================================
# JET Configuration
# ============================================
# JET Singularity Image Path
export JET2="/home/sasidharp/jet_docker/jet.sif"  # REQUIRED: Path to JET singularity image (e.g., /home/user/jet_docker/jet.sif)
export JET2_localPath="/home/junseokp/workspaces/projects/rtea/tools/JET_identification_pipeline"

# Tool Paths (within JET singularity image)
export samtoolsBinDir="/usr/local/bin"  # Path within singularity image
export starBinDir="/usr/local/bin"  # Path within singularity image

# Directories
export outputDir=${OUTPUT_BASE}
export dataDir=${DATA_HOME}
export refDir="${REF_DIR}/hg38"

# Sequencing Parameters
export readLength=150  # UPDATE if your read length is different
export organism="human"  # Options: human, mouse, etc.
export genome="hg38"  # Genome version
export database="database_name"  # UPDATE: Database name for JET
export threads=8

# JET Sample Parameters (set these when running JET functions)
# These should be set per sample in your processing scripts:
# export rnaSample="path/to/sample.fastq.gz"  # Path to RNA sample file
# export name="sample_name"                    # Sample name
# export namePrefix="prefix"                   # Sample name prefix
# export outputDir="${TEST_LOG_DIR}/test_output"
# export SAMPLE_NAME="sim200_AluY_blood"
# export rnaSample="nonReferenceTE/AluY/5X/fq"

# Reference Files for JET
export fastaFile="Homo_sapiens_assembly38.fasta"  # REQUIRED
export gtfGeneFile="gencode.v46.annotation.gtf"  # REQUIRED
export starIndexesDir="star/idx"  # REQUIRED: STAR index directory
export repeatsFile="${REF_DIR}/hg38.RepeatMasker-4.0.6-Dfam-2.0.reformat.txt"  # REQUIRED: Repeat elements file - Check Point
export gffFile="${REF_DIR}/hg38.RepeatMasker-4.0.6-Dfam-2.0.annotation.gff"  # REQUIRED: TE annotation in GFF format

# R Configuration (within JET singularity image)
export RlibDir="/usr/local/lib/R/library"  # Path within singularity image

# ============================================
# TEProf2 Configuration
# ============================================

export TEProf2="/home/sasidharp/jet_docker/teprof2.sif"  # Singularity container path

# Reference Files for TEProf2 (can be same as JET if compatible)
export TEPROF2_REF="${REF_DIR}/Homo_sapiens_assembly38.fasta"
export TEPROF2_TE_ANNOT="${REF_DIR}/TE_annotation.gtf"
export TEPROF2_GENE_ANNOT="${REF_DIR}/gencode.v46.annotation.gtf"

# TEProf2 Pipeline-specific References (REQUIRED for per-sample processing)
export STAR_INDEX="${REF_DIR}/STAR_hg38_index"  # REQUIRED: Prebuilt STAR index directory for TEProf2
export GENCODE_GTF="${REF_DIR}/gencode.gtf"     # REQUIRED: GENCODE annotation GTF for StringTie
export ARGUMENTS_TXT="${REF_DIR}/arguments.txt" # REQUIRED: TEProf2 arguments.txt file containing reference paths

# ============================================
# TEProf2 Aggregation (run_command.sh) Configuration
# ============================================
# After all individual samples are processed, run_command.sh aggregates results
# across all samples to identify TE-derived transcripts.
#
# IMPORTANT: The aggregation step requires additional TEProf2 tools to be available:
# - aggregateProcessedAnnotation.R, filterReadCandidates.R, mergeAnnotationProcess.R, finalStatisticsOutput.R
# - rmskhg38_annotate_gtf_update_test_tpm_cuff.py, commandsmax_speed.py, stringtieExpressionFrac.py
# - samtools, stringtie, gffread, cuffmerge
#
# These tools should either be:
# 1. In your PATH, OR
# 2. Run the aggregation within the TEProf2 singularity container, e.g.:
#    singularity exec --bind /path/to/data:/path/to/data $TEProf2 bash run_teprof2_aggregation.sh

# REQUIRED for run_command.sh:
# The ARGUMENTS_TXT file (above) must contain these tab-delimited entries:
# rmsk	<path_to_repeatmasker.bed.gz>       - Tabix-formatted RepeatMasker BED6 file
# rmskannotationfile	<path_to_description.lst>  - Tab-delimited: subfamily, class, family
# gencodeplusdic	<path_to_gencode_plus.dic>    - Dictionary for plus strand elements
# gencodeminusdic	<path_to_gencode_minus.dic>  - Dictionary for minus strand elements

# Optional entries for ARGUMENTS_TXT:
# focusgenes	<path_to_gene_list.txt>   - Filter to specific genes (one per line)
# plusintron	<path_to_introns_plus.gz>  - Tabix file of plus strand introns
# minusintron	<path_to_introns_minus.gz> - Tabix file of minus strand introns

# Location where run_command.sh should be executed
# This should be the directory containing all processed TEProf2 sample subdirectories
export TEPROF2_AGGREGATION_DIR="${OUTPUT_BASE}/TEProf2_aggregated"

# Path to run_command.sh script (included in this repository)
export RUN_COMMAND_SCRIPT="${SCRIPT_DIR:-/path/to/scripts}/run_command.sh"

# Reference GTF for cuffmerge step in run_command.sh
# This path is hardcoded in run_command.sh as ../genome_46/gencode.v46.basic.annotation.sorted.gtf
# Update this to match your actual reference location
export TEPROF2_CUFFMERGE_GTF="${REF_DIR}/gencode.v46.basic.annotation.sorted.gtf"

# ============================================
# TEProf2 Optional Arguments for aggregateProcessedAnnotation.R
# ============================================
# These parameters control filtering during the aggregation step

# Treatment label for identifying treatment samples (default: '')
# If not specified, all samples are considered as treatment
export TEPROF2_AGG_EXT_TREATMENT=""

# Maximum length of exon 1 in bp (default: 2588)
# Based on 99th percentile of GENCODE v25 transcripts
export TEPROF2_AGG_EXON1_LENGTH_MAX=2588

# Maximum exon skipping events allowed (default: 2)
# Helps filter out potential genomic contamination/intron retention artifacts
export TEPROF2_AGG_EXON_SKIP_MAX=2

# Minimum number of samples that must contain a candidate (default: 1)
# Increase for larger studies to focus on recurrent events
export TEPROF2_AGG_SAMPLE_TOTAL_MIN=1

# Minimum number of treatment samples required (default: 0)
# Useful for treatment-specific analysis
export TEPROF2_AGG_TREATMENT_TOTAL_MIN=0

# Treatment-exclusive filtering (default: no)
# Set to 'yes' to only keep candidates absent in untreated samples
# Options: yes, no
export TEPROF2_AGG_TREATMENT_EXCLUSIVE="no"

# Keep transcripts that don't splice into genes (default: no)
# TE-derived transcripts without gene splicing are filtered by default
# Options: yes, no
export TEPROF2_AGG_KEEP_NONE="no"

# Filter for TEs only, removing other repeats (default: yes)
# RepeatMasker includes non-TE repeats; this filters to TEs only
# Options: yes, no
# NOTE: Setting to 'no' may break downstream steps
export TEPROF2_AGG_FILTER_FOR_TES="yes"

# ============================================
# TEProf2 Optional Arguments for filterReadCandidates.R
# ============================================
# These parameters filter candidates based on read support

# Minimum paired-end reads starting within TE per file (default: 10)
export TEPROF2_FILTER_MIN_READS_IN_TE=10

# Minimum paired-end reads spanning TE-to-gene distance across all files (default: 1)
export TEPROF2_FILTER_MIN_START_READ=1

# Maximum percentage of files with reads suggesting exonization (default: 0.15)
# Reads starting in gene and going to TE suggest exonization rather than promoter activity
export TEPROF2_FILTER_EXONIZATION_MAX_PERCENT=0.15

# Minimum distance between TE and transcript start in bp (default: 2500)
# Helps remove noise from annotated promoters
export TEPROF2_FILTER_DISTANCE_TE=2500

# ============================================
# TEProf2 StringTie Quantification Parameters
# ============================================
# Parameters used during transcript-level expression quantification

# Minimum read coverage (default: 1)
export TEPROF2_STRINGTIE_MIN_COVERAGE=1

# Minimum assembled transcript length in bp (default: 100)
export TEPROF2_STRINGTIE_MIN_TRANSCRIPT_LENGTH=100

# Number of threads for StringTie (default: 2)
export TEPROF2_STRINGTIE_THREADS=2

# MAPQ filter for uniquely mapped reads (default: 255 for STAR, 60 for HISAT2)
# TEProf2 processes repetitive elements, so unique mapping is important
export TEPROF2_MAPQ_FILTER=255

# TEProf2 Basic Parameters (for per-sample processing)
export TEPROF2_MIN_MAPQ=20
export TEPROF2_MIN_BASE_QUAL=20

# ============================================
# Compute Resources
# ============================================

export THREADS=8
export MEMORY="32G"
export WALLTIME="24:00:00"
export MAX_CONCURRENT_JOBS=20  # For array jobs: --array=1-N%20

# ============================================
# Module Loading (if needed)
# ============================================

# Uncomment and modify if you need to load modules
# module load singularity
# module load samtools
# module load STAR
# module load R

# ============================================
# Validation
# ============================================

# Function to validate configuration
validate_config() {
    local errors=0
    
    echo "Validating configuration..."
    
    # Check JET singularity image
    if [ ! -f "$JET2" ]; then
        echo "ERROR: JET2 singularity image not found: $JET2"
        errors=$((errors + 1))
    fi
    
    # Check reference files
    if [ ! -f "$fastaFile" ]; then
        echo "ERROR: Reference FASTA not found: $fastaFile"
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$gtfGeneFile" ]; then
        echo "ERROR: Gene GTF not found: $gtfGeneFile"
        errors=$((errors + 1))
    fi
    
    if [ ! -d "$starIndexesDir" ]; then
        echo "ERROR: STAR index directory not found: $starIndexesDir"
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$repeatsFile" ]; then
        echo "ERROR: Repeats file not found: $repeatsFile"
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$gffFile" ]; then
        echo "ERROR: GFF file not found: $gffFile"
        errors=$((errors + 1))
    fi
    
    # Check TEProf2
    if [ ! -f "$TEProf2" ]; then
        echo "ERROR: TEProf2 container not found: $TEProf2"
        errors=$((errors + 1))
    fi
    
    # Check TEProf2-specific reference files
    if [ ! -d "$STAR_INDEX" ]; then
        echo "ERROR: STAR index directory not found: $STAR_INDEX"
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$GENCODE_GTF" ]; then
        echo "ERROR: GENCODE GTF not found: $GENCODE_GTF"
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$ARGUMENTS_TXT" ]; then
        echo "ERROR: TEProf2 arguments.txt not found: $ARGUMENTS_TXT"
        errors=$((errors + 1))
    fi
    
    # Check TEProf2 aggregation requirements
    if [ ! -f "$TEPROF2_CUFFMERGE_GTF" ]; then
        echo "WARNING: TEPROF2_CUFFMERGE_GTF not found: $TEPROF2_CUFFMERGE_GTF"
        echo "         This is required for TEProf2 aggregation step"
    fi
    
    # Check data directory
    if [ ! -d "$DATA_HOME" ]; then
        echo "ERROR: DATA_HOME not found: $DATA_HOME"
        errors=$((errors + 1))
    fi
    
    if [ $errors -eq 0 ]; then
        echo "✓ Configuration validation passed!"
        return 0
    else
        echo "✗ Configuration validation failed with $errors error(s)"
        return 1
    fi
}

# ============================================
# Usage Example
# ============================================
# 
# 1. Copy this file: cp config_template.sh config.sh
# 2. Edit config.sh with your actual paths
# 3. Source in your scripts: source config.sh
# 4. Validate: validate_config
#
# ============================================

echo "Configuration template loaded."
echo "Remember to:"
echo "  1. Copy to config.sh"
echo "  2. Update all paths marked REQUIRED"
echo "  3. Run validate_config to check"
