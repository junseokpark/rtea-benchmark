# TE Analysis Pipeline - Complete Package

## Overview

This package contains a complete pipeline for processing transposable element (TE) sequencing data using JET and TEProf2 tools. The pipeline handles ~200 samples across different TE types, tissue types, and coverage levels.

## Package Contents

All scripts are located in the `scripts/` directory.

### Main Scripts

1. **pipeline.sh** - Master control script
   - Simplifies the entire workflow with easy commands
   - Provides color-coded output and status checking
   - Usage: `./pipeline.sh [command]`
   - Commands: help, setup, check-config, list, submit, status, resubmit, validate

2. **process_array.sh** - SLURM array job script (RECOMMENDED)
   - Processes all samples in parallel using STAR aligner
   - Uses array jobs for efficient HPC resource utilization
   - Automatically scales to your sample count
   - Runs JET (Step 1 & 2) and TEProf2 on each sample

3. **process_samples.sh** - Sequential processing script
   - Alternative to array jobs
   - Processes samples one by one
   - Useful for smaller datasets or non-HPC environments

### Configuration

4. **config_template.sh** - Configuration template
   - Template for environment-specific settings
   - Copy to config.sh and customize paths
   - Contains JET and TEProf2 configuration
   - Includes validation function

### Helper Scripts

5. **generate_sample_list.sh** - Sample discovery
   - Scans your data directory
   - Creates sample_list.txt with all paired-end FASTQ files
   - Organizes by TE type, tissue, and coverage

6. **check_status.sh** - Progress monitoring
   - Generates detailed status reports
   - Tracks JET Step 1 (STAR alignment) separately from Step 2 (R analysis)
   - Shows completion percentage for each tool
   - Breaks down results by category and coverage
   - Identifies failed samples

7. **resubmit_failed.sh** - Failure recovery
   - Identifies incomplete or failed samples
   - Creates a resubmission script automatically
   - Allows targeted reprocessing without redoing successful samples

8. **runTools.sh** - Tool execution wrapper
   - Helper script for running tools in containers

### Documentation

7. **README.md** - Comprehensive documentation
   - Detailed usage instructions
   - Configuration guidelines with config_template.sh
   - Troubleshooting tips
   - Expected output formats

9. **SUMMARY.md** (this file) - Complete package overview
   - Lists all scripts and their purposes
   - Quick start guide
   - Resource requirements

## Quick Start (4 Steps)

### Step 1: Configure
```bash
cd scripts
cp config_template.sh config.sh
# Edit config.sh with your paths
nano config.sh
```

### Step 2: Setup
```bash
chmod +x pipeline.sh
./pipeline.sh setup
./pipeline.sh check-config
```

### Step 3: Generate Sample List and Submit
```bash
./pipeline.sh list
./pipeline.sh submit
```

### Step 4: Monitor and Verify
```bash
# Check status periodically
./pipeline.sh status

# After completion
./pipeline.sh validate
```

## Directory Structure

### Input Data
```
/home/junseokp/workspaces/data/rTea-simul/
├── nonReferenceTE/
│   ├── AluY/
│   ├── L1HS/
│   ├── LTR5/
│   └── SVA_F/
│       └── {5X,10X,50X,100X,200X}/fq/*.fq.gz
└── referenceTE/
    ├── intron/{coverage}/fq/*.fq.gz
    └── TSS/{coverage}/fq/*.fq.gz
```

### Output Results
```
/home/junseokp/workspaces/data/rTea-simul/output/
├── nonReferenceTE/{TE}/{coverage}/
│   ├── JET/{sample}/
│   │   ├── {sample}.sorted.bam
│   │   ├── polymorphic_insertions.bed
│   │   └── TE_expression.tsv
│   └── TEProf2/{sample}/
│       ├── TE_counts.tsv
│       └── insertion_sites.bed
└── referenceTE/{type}/{coverage}/{fq|fq_mut}/
    └── [same structure]
```

## Sample Categories

### nonReferenceTE (100 samples)
- **TE types:** AluY, L1HS, LTR5, SVA_F
- **Tissues:** blood, brain, colon, esophagus, lung
- **Coverage:** 5X, 10X, 50X, 100X, 200X
- **Pattern:** sim200_{TE}_{tissue}.{1,2}.fq.gz

### referenceTE/intron (50 samples)
- **Replicates:** 1-5
- **Coverage:** 5X, 10X, 50X, 100X, 200X
- **Conditions:** regular + mutated
- **Pattern:** reftefu403_{mut_}_{replicate}_{coverage}.{1,2}.fq.gz

### referenceTE/TSS (50 samples)
- **Replicates:** 1-5
- **Coverage:** 5X, 10X, 50X, 100X, 200X
- **Conditions:** regular + mutated
- **Pattern:** reftetss270_{mut_}_{replicate}_{coverage}.{1,2}.fq.gz

## Resource Requirements

### Per Sample
- **Memory:** 32 GB
- **CPUs:** 8 cores
- **Time:** 2-7 hours
- **Disk:** 5-20 GB

### Total Project
- **Samples:** ~200
- **Disk Space:** 1-4 TB
- **CPU Hours:** 3,200-11,200
- **Wall Time (array):** 10-20 hours
- **Wall Time (sequential):** 17-58 days

## Workflow Diagram

```
┌─────────────────┐
│  Generate List  │  ./pipeline.sh list
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Submit Jobs    │  ./pipeline.sh submit
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Processing     │  JET + TEProf2 (parallel)
│  (10-20 hours)  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Check Status   │  ./pipeline.sh status
└────────┬────────┘
         │
         ├─── All Success ────────▶ Done!
         │
         └─── Some Failed ────────▶ ./pipeline.sh resubmit
```

## Key Features

✓ **Automated Sample Discovery** - Finds all FASTQ pairs automatically
✓ **Parallel Processing** - Uses SLURM array jobs for efficiency
✓ **Progress Tracking** - Real-time status monitoring
✓ **Error Recovery** - Automatic identification and resubmission of failures
✓ **Organized Output** - Results mirror input structure
✓ **Validation** - Built-in checks for output integrity
✓ **Comprehensive Logging** - Detailed logs for each sample

## Configuration

Before running, verify these paths in the scripts:

```bash
DATA_HOME=/home/junseokp/workspaces/data/rTea-simul
REF=${DATA_HOME}/ref
OUTPUT_BASE=${DATA_HOME}/output

JET=/home/sasidharp/jet_docker/jet.sif
TEProf2=/home/sasidharp/jet_docker/teprof2.sif
```

Required reference files:
- `${REF}/reference.fa`
- `${REF}/TE_annotation.bed` or `.gtf`
- `${REF}/gene_annotation.gtf`

## Common Commands

```bash
# Initial setup
./pipeline.sh setup
./pipeline.sh check-config

# Generate sample list
./pipeline.sh list

# Submit processing jobs
./pipeline.sh submit         # Array job (recommended)
./pipeline.sh submit-seq     # Sequential (slower)

# Monitor progress
squeue -u $USER              # Check job queue
./pipeline.sh status         # Detailed status report
./pipeline.sh summary        # Quick summary

# Handle failures
./pipeline.sh resubmit       # Resubmit failed samples

# Validation
./pipeline.sh validate       # Verify outputs

# Cleanup (use with caution!)
./pipeline.sh clean-logs     # Remove log files
./pipeline.sh clean          # Remove all outputs
```

## Troubleshooting

### Jobs Not Starting
- Check cluster load: `sinfo`, `squeue`
- Verify SLURM account and partition settings

### Out of Memory
- Increase `--mem` in script headers
- Reduce concurrent jobs (lower `%N` in --array)

### Missing Dependencies
- Load singularity: `module load singularity`
- Verify container paths exist
- Check reference file paths

### Tool Errors
- Check tool-specific logs in `logs/` directory
- Refer to tool documentation:
  - JET: https://github.com/junseokpark/JET_identification_pipeline
  - TEProf2: https://github.com/junseokpark/TEProf2Paper

## Support

For pipeline issues or questions:
1. Check QUICKSTART.txt for common workflows
2. Review README.md for detailed documentation
3. Examine log files in `logs/` directory
4. Verify configuration with `./pipeline.sh check-config`

## Version

- Pipeline Version: 1.0
- Created: November 2025
- Compatible with: SLURM-based HPC systems
- Tools: JET, TEProf2
- Container: Singularity/Apptainer

---

**Ready to begin?** Run: `./pipeline.sh help`
