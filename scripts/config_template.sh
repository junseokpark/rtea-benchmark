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
export JET2_localPath="/home/junseokp/workspaces/tools/JET_identification_pipeline"

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

# TEProf2 Pipeline-specific References
export STAR_INDEX="${REF_DIR}/STAR_hg38_index"  # REQUIRED: Prebuilt STAR index directory for TEProf2
export GENCODE_GTF="${REF_DIR}/gencode.gtf"     # REQUIRED: GENCODE annotation GTF for StringTie
export ARGUMENTS_TXT="${REF_DIR}/arguments.txt" # REQUIRED: TEProf2 arguments.txt file

# TEProf2 Parameters
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
