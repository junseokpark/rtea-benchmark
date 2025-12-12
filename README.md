# TE Analysis Pipeline - Processing Guide

This pipeline processes transposable element (TE) sequencing data using JET and TEProf2 tools.

## Directory Structure

```
DATA_HOME/
├── nonReferenceTE/
│   ├── AluY/
│   ├── L1HS/
│   ├── LTR5/
│   └── SVA_F/
│       └── {5X,10X,50X,100X,200X}/
│           └── fq/
│               └── sim200_{TE}_{tissue}.{1,2}.fq.gz
└── referenceTE/
    ├── intron/
    │   └── {5X,10X,50X,100X,200X}/
    │       ├── fq/
    │       │   └── reftefu403_{replicate}_{coverage}.{1,2}.fq.gz
    │       └── fq_mut/
    │           └── reftefu403_mut_{replicate}_{coverage}.{1,2}.fq.gz
    └── TSS/
        └── {5X,10X,50X,100X,200X}/
            ├── fq/
            │   └── reftetss270_{replicate}_{coverage}.{1,2}.fq.gz
            └── fq_mut/
                └── reftetss270_mut_{replicate}_{coverage}.{1,2}.fq.gz
```

## Output Structure

Results will be organized as:
```
output/
├── nonReferenceTE/{TE_type}/{coverage}/
│   ├── JET/{sample_name}/
│   └── TEProf2/{sample_name}/
└── referenceTE/{intron|TSS}/{coverage}/{fq|fq_mut}/
    ├── JET/{sample_name}/
    └── TEProf2/{sample_name}/
```

## Prerequisites

1. **Required modules:**
   - singularity

2. **Required singularity containers:**
   - JET: `/home/sasidharp/jet_docker/jet.sif`
   - TEProf2: `/home/sasidharp/jet_docker/teprof2.sif`

3. **Reference files:**
   - `${REF}/reference.fa` - Reference genome
   - `${REF}/TE_annotation.gff` - TE annotations in GFF format
   - `${REF}/gene_annotation.gtf` - Gene annotations
   - `${REF}/repeats.txt` - Repeat elements file
   - `${REF}/STAR_indexes/` - STAR genome indexes directory

4. **JET dependencies:**
   - JET installation directory
   - Samtools binary directory
   - STAR binary directory
   - R library path

## Configuration

All pipeline scripts now load their configuration from a shared `config.sh` file, which centralizes all environment variables and paths in one place.

Before running the pipeline, you should configure the paths specific to your environment:

**Step 1: Create configuration file**
```bash
cd scripts
cp config_template.sh config.sh
```

**Step 2: Edit config.sh and update these required paths:**
- `DATA_HOME` - Base directory for input data
- `JETProjectDir` - Path to JET installation
- `samtoolsBinDir` - Path to samtools binary directory
- `starBinDir` - Path to STAR binary directory
- `RlibDir` - Path to R library for JET Step 2
- `TEProf2` - Path to TEProf2 singularity container
- Update reference file paths if different from defaults
- Adjust sequencing parameters (read length, organism, genome)

**Step 3: Validate configuration**
```bash
source config.sh
validate_config
```

**Note:** The `config.sh` file is user-specific and is excluded from version control (via `.gitignore`). This allows you to maintain your custom configuration without affecting the repository. All pipeline scripts (`check_status.sh`, `generate_sample_list.sh`, `process_array.sh`, `process_samples.sh`, `resubmit_failed.sh`, `pipeline.sh`, and `runTools.sh`) automatically source this configuration file.

## Usage

### Option 1: Array Job (Recommended for HPC)

This approach runs samples in parallel using SLURM array jobs.

**Step 1: Generate sample list**
```bash
chmod +x generate_sample_list.sh
./generate_sample_list.sh
```

This creates `sample_list.txt` with all samples to process.

**Step 2: Adjust array size**
Check the number of samples:
```bash
wc -l sample_list.txt
```

Edit `process_array.sh` and update the array size:
```bash
#SBATCH --array=1-N%20
# where N = total number of samples (minus 1 for header)
# %20 limits to 20 concurrent jobs
```

**Step 3: Submit array job**
```bash
chmod +x process_array.sh
sbatch process_array.sh
```

**Monitor jobs:**
```bash
squeue -u $USER
sacct -j JOBID --format=JobID,JobName,State,ExitCode
```

### Option 2: Sequential Processing

Process all samples one by one (slower but simpler):

```bash
chmod +x process_samples.sh
sbatch process_samples.sh
# or run interactively:
# ./process_samples.sh
```

## Customization

### Adjust Resource Requirements

Edit SBATCH parameters in the scripts:
```bash
#SBATCH --mem=32G          # Memory per job
#SBATCH --cpus-per-task=8  # CPU cores
#SBATCH --time=24:00:00    # Time limit
```

### Modify Tool Parameters

**For JET:**
Edit the `run_jet_step1()` and `run_jet_step2()` functions in the scripts to adjust:
- STAR alignment parameters
- Read length and organism settings
- Detection thresholds
- Output formats

**For TEProf2:**
Edit the `run_teprof2()` function to adjust:
- MAPQ threshold (--min-mapq)
- Base quality threshold (--min-base-quality)
- Other tool-specific parameters

### Update Reference Paths

Configuration paths are now managed in `scripts/config.sh`. See the Configuration section above for details.

## Monitoring and Troubleshooting

### Check Progress

Use the improved status checker that monitors both JET steps separately:
```bash
cd scripts
chmod +x check_status.sh
./check_status.sh
```

The status report will show:
- JET Step 1 (STAR alignment) completion status
- JET Step 2 (R analysis) completion status
- TEProf2 completion status
- Breakdown by TE type and coverage
- List of failed samples

You can also check manually:
```bash
# Count completed alignments (JET Step 1)
find output -name "Aligned.sortedByCoord.out.bam" | wc -l

# Count completed JET Step 2 results
find output -name "*_TE_insertions.bed" | wc -l

# Count TEProf2 results
find output -type d -name "TEProf2" | wc -l

# Check for errors in logs
grep -i "error" logs/*.err
grep -i "failed" logs/*.err
```

### Common Issues

1. **Missing reference files:**
   - Ensure all reference files exist in `${REF}/` directory
   - Check file permissions

2. **Singularity container not found:**
   - Verify container paths
   - Check singularity module is loaded

3. **Out of memory:**
   - Increase `--mem` in SBATCH directives
   - Reduce number of concurrent jobs (lower the %N in --array)

4. **Input files not found:**
   - Verify DATA_HOME path
   - Check file naming conventions match the script patterns

## Sample Count Summary

Based on the directory structure:

- **nonReferenceTE:**
  - 4 TE types × 5 coverage levels × 5 tissues = 100 samples
  
- **referenceTE/intron:**
  - 5 coverage levels × 5 replicates × 2 conditions (regular + mutated) = 50 samples
  
- **referenceTE/TSS:**
  - 5 coverage levels × 5 replicates × 2 conditions (regular + mutated) = 50 samples

**Total: ~200 samples** (each with paired-end reads)

## Expected Runtime

- Per sample (estimate):
  - JET: 1-4 hours depending on coverage
  - TEProf2: 1-3 hours depending on coverage
  
- Total with array job (20 concurrent): ~10-20 hours
- Sequential processing: ~400-1400 hours

## Output Files

### JET outputs:
```
output/{path}/JET/{sample}/
├── {sample}.sorted.bam
├── {sample}.sorted.bam.bai
├── polymorphic_insertions.bed
├── TE_expression.tsv
└── summary_statistics.txt
```

### TEProf2 outputs:
```
output/{path}/TEProf2/{sample}/
├── TE_counts.tsv
├── TE_expression.tsv
├── insertion_sites.bed
└── quality_metrics.txt
```

## Support

For tool-specific issues:
- JET: https://github.com/junseokpark/JET_identification_pipeline
- TEProf2: https://github.com/junseokpark/TEProf2Paper

## Version Info

- Pipeline version: 1.0
- Last updated: 2025
