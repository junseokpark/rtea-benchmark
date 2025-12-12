#!/bin/bash
# JET and TEProf2 Pipeline Configuration Template
# 
# This is a template file. To use it:
#   1. Copy this file to config.sh: cp config_template.sh config.sh
#   2. Edit config.sh and update all paths marked with "UPDATE" or "REQUIRED"
#   3. All pipeline scripts will automatically source config.sh
#
# NOTE: Do not edit this template directly. Always work with config.sh
# so your custom configuration is not overwritten by repository updates.

# ============================================
# REQUIRED: UPDATE ALL PATHS BELOW
# ============================================

# Data Directories
export DATA_HOME="/home/junseokp/workspaces/data/rTea-simul"
export OUTPUT_BASE="${DATA_HOME}/output"
export REF_DIR="${DATA_HOME}/ref"

# Sample List File
export SAMPLE_LIST="sample_list.txt"

# ============================================
# JET Configuration
# ============================================

# JET Installation Path
export JETProjectDir="/path/to/JET"  # REQUIRED: Path to JET installation

# Tool Paths
export samtoolsBinDir="/path/to/samtools/bin"  # REQUIRED
export starBinDir="/path/to/STAR/bin"  # REQUIRED

# Sequencing Parameters
export readLength=150  # UPDATE if your read length is different
export organism="human"  # Options: human, mouse, etc.
export genome="hg38"  # Genome version
export database="database_name"  # UPDATE: Database name for JET

# Reference Files for JET
export refDir="${REF_DIR}"  # Reference directory (for compatibility with existing scripts)
export fastaFile="${REF_DIR}/reference.fa"  # REQUIRED
export gtfGeneFile="${REF_DIR}/gene_annotation.gtf"  # REQUIRED
export starIndexesDir="${REF_DIR}/STAR_indexes"  # REQUIRED: STAR index directory
export repeatsFile="${REF_DIR}/repeats.txt"  # REQUIRED: Repeat elements file
export gffFile="${REF_DIR}/TE_annotation.gff"  # REQUIRED: TE annotation in GFF format

# R Configuration
export RlibDir="/path/to/R/library"  # REQUIRED: R library path for JET Step 2

# ============================================
# Singularity Container Paths
# ============================================

export JET="/home/sasidharp/jet_docker/jet.sif"  # JET Singularity container path
export TEProf2="/home/sasidharp/jet_docker/teprof2.sif"  # TEProf2 Singularity container path

# ============================================
# TEProf2 Parameters
# ============================================

# Reference Files for TEProf2 (can be same as JET if compatible)
export TEPROF2_REF="${REF_DIR}/reference.fa"
export TEPROF2_TE_ANNOT="${REF_DIR}/TE_annotation.gtf"
export TEPROF2_GENE_ANNOT="${REF_DIR}/gene_annotation.gtf"

# TEProf2 Parameters
export TEPROF2_MIN_MAPQ=20
export TEPROF2_MIN_BASE_QUAL=20

# ============================================
# Compute Resources
# ============================================

export threads=8  # Thread count (lowercase for compatibility with existing scripts)
export THREADS=8  # Thread count (uppercase alternative)
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
    
    # Check JET paths
    if [ ! -d "$JETProjectDir" ]; then
        echo "ERROR: JETProjectDir not found: $JETProjectDir"
        errors=$((errors + 1))
    fi
    
    if [ ! -d "$samtoolsBinDir" ]; then
        echo "ERROR: samtoolsBinDir not found: $samtoolsBinDir"
        errors=$((errors + 1))
    fi
    
    if [ ! -d "$starBinDir" ]; then
        echo "ERROR: starBinDir not found: $starBinDir"
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
    
    if [ ! -d "$RlibDir" ]; then
        echo "ERROR: R library directory not found: $RlibDir"
        errors=$((errors + 1))
    fi
    
    # Check TEProf2
    if [ ! -f "$TEProf2" ]; then
        echo "ERROR: TEProf2 container not found: $TEProf2"
        errors=$((errors + 1))
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
