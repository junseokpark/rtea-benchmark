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

4. **TEProf2-specific reference files (required for aggregation):**
   - `${REF}/arguments.txt` - Tab-delimited file with paths to TEProf2 references
   - RepeatMasker bed6.gz file (tabix-indexed)
   - RepeatMasker annotation file (subfamily, class, family mapping)
   - GENCODE plus/minus strand dictionaries
   - Optional: Gene filter list, intron annotations

   **arguments.txt format:**
   ```
   rmsk	/path/to/repeatmasker.bed.gz
   rmskannotationfile	/path/to/repeatmasker_description_uniq.lst
   gencodeplusdic	/path/to/genecode_plus.dic
   gencodeminusdic	/path/to/genecode_minus.dic
   focusgenes	/path/to/gene_list.txt (optional)
   plusintron	/path/to/introns_plus.gz (optional)
   minusintron	/path/to/introns_minus.gz (optional)
   ```

5. **JET dependencies:**
   - JET singularity image (e.g., `/home/sasidharp/jet_docker/jet.sif`)
   - Tool paths are configured within the singularity image

## Configuration

Before running the pipeline, you should configure the paths specific to your environment:

**Step 1: Create configuration file**
```bash
cd scripts
cp config_template.sh config.sh
```

**Step 2: Edit config.sh and update these required paths:**
- `JET2` - Path to JET singularity image (e.g., /path/to/jet.sif)
- `samtoolsBinDir` - Path to samtools binary directory within singularity image (default: /opt/samtools/bin)
- `starBinDir` - Path to STAR binary directory within singularity image (default: /opt/STAR/bin)
- `RlibDir` - Path to R library within singularity image (default: /opt/R/library)
- Update reference file paths if different from defaults
- Adjust sequencing parameters (read length, organism, genome)

**Step 3: Validate configuration**
```bash
source config.sh
validate_config
```

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
# Subtract 1 for header to get total sample count
```

**Option A: Process one sample per job (default)**

Edit `process_array.sh` and update the array size:
```bash
#SBATCH --array=1-N%20
# where N = total number of samples (minus 1 for header)
# %20 limits to 20 concurrent jobs
```

**Option B: Process multiple samples per job (chunking)**

Use sample chunking to process multiple samples per job sequentially. This is useful for:
- Reducing the number of jobs submitted to the scheduler
- Grouping samples to match partition requirements (e.g., defq15*, defq610)
- Better workload control and resource utilization

**Use the calculator helper script (recommended):**
```bash
# Make the script executable (first time only)
chmod +x calculate_array_size.sh

# Calculate the required array size automatically
./calculate_array_size.sh 20  # For 20 samples per job

# This will show:
# - Total samples and jobs needed
# - SBATCH --array configuration
# - SAMPLES_PER_JOB setting
# - Sample distribution across jobs
```

**Or calculate manually:**
```bash
# Formula: Number of jobs = ceiling(total_samples / samples_per_job)
# Example: 100 samples with 20 samples per job = 5 jobs
```

Configure via batch-config.sh:
```bash
cp batch-config.sh.template batch-config.sh
# Edit batch-config.sh:
export SAMPLES_PER_JOB="20"  # Process 20 samples per job
export BATCH_ARRAY="1-5"      # 100 samples / 20 per job = 5 jobs
```

Or override at submission time:
```bash
# For 100 samples with 20 samples per job:
sbatch --export=ALL,SAMPLES_PER_JOB=20 process_array.sh
# Note: You still need to adjust #SBATCH --array in the script to 1-5
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

**Sample Chunking Examples:**
```bash
# Example 1: 100 samples, 20 per job = 5 jobs
export SAMPLES_PER_JOB="20"
# Set #SBATCH --array=1-5%5 in process_array.sh

# Example 2: 200 samples, 10 per job = 20 jobs  
export SAMPLES_PER_JOB="10"
# Set #SBATCH --array=1-20%10 in process_array.sh

# Example 3: 150 samples, 25 per job = 6 jobs
export SAMPLES_PER_JOB="25"
# Set #SBATCH --array=1-6%6 in process_array.sh
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
TEProf2 runs in two phases:

1. **Per-sample processing** (handled by `run_teprof2()` in function.sh):
   - STAR alignment
   - StringTie transcript assembly
   - RepeatMasker annotation
   - TPM processing
   - Read candidate filtering

2. **Cross-sample aggregation** (handled by `run_teprof2_aggregation()` via run_command.sh):
   - Aggregates results across all samples
   - Identifies TE-derived transcripts
   - Performs read-based filtering
   - Merges candidates with reference annotations
   - Calculates expression levels

**TEProf2 Configuration Options:**

All TEProf2 configuration options are documented in `config_template.sh` with inline comments. Key options include:

**Required References:**
- `ARGUMENTS_TXT`: Path to arguments.txt file (contains paths to RepeatMasker, GENCODE dictionaries, etc.)
- `STAR_INDEX`: STAR genome index directory
- `GENCODE_GTF`: GENCODE annotation for StringTie and cuffmerge
- `TEPROF2_CUFFMERGE_GTF`: Reference GTF for merging candidates

**Aggregation Filtering Options:**
- `TEPROF2_AGG_EXON1_LENGTH_MAX`: Maximum exon 1 length (default: 2588)
- `TEPROF2_AGG_EXON_SKIP_MAX`: Maximum exon skipping events (default: 2)
- `TEPROF2_AGG_SAMPLE_TOTAL_MIN`: Minimum samples per candidate (default: 1)
- `TEPROF2_AGG_TREATMENT_TOTAL_MIN`: Minimum treatment samples (default: 0)
- `TEPROF2_AGG_TREATMENT_EXCLUSIVE`: Keep only treatment-exclusive candidates (default: no)
- `TEPROF2_AGG_KEEP_NONE`: Keep transcripts without gene splicing (default: no)
- `TEPROF2_AGG_FILTER_FOR_TES`: Filter to TEs only, removing other repeats (default: yes)

**Read Filtering Options:**
- `TEPROF2_FILTER_MIN_READS_IN_TE`: Minimum reads starting in TE (default: 10)
- `TEPROF2_FILTER_MIN_START_READ`: Minimum spanning reads (default: 1)
- `TEPROF2_FILTER_EXONIZATION_MAX_PERCENT`: Maximum exonization percentage (default: 0.15)
- `TEPROF2_FILTER_DISTANCE_TE`: Minimum TE-to-gene distance (default: 2500)

**StringTie Quantification:**
- `TEPROF2_STRINGTIE_MIN_COVERAGE`: Minimum coverage (default: 1)
- `TEPROF2_STRINGTIE_MIN_TRANSCRIPT_LENGTH`: Minimum length (default: 100)
- `TEPROF2_MAPQ_FILTER`: MAPQ for unique mapping (default: 255 for STAR, 60 for HISAT2)

**Running TEProf2 Aggregation:**

The aggregation step runs automatically after all samples in `process_samples.sh`.

To run aggregation separately:
```bash
# After all samples are processed:
cd scripts
./run_teprof2_aggregation.sh [output_base_dir]

# Or submit as SLURM job:
sbatch run_teprof2_aggregation.sh

# If tools are not in PATH, run within TEProf2 container:
singularity exec --bind ${OUTPUT_BASE}:${OUTPUT_BASE} \
  --bind ${REF_DIR}:${REF_DIR} \
  ${TEProf2} bash run_teprof2_aggregation.sh
```

**Important Notes:**
- The aggregation step requires TEProf2 R scripts and Python tools to be available in PATH
- If using a container, ensure all paths (output directories, reference files) are bind-mounted
- The `arguments.txt` file must contain correct paths to all reference files
- BAM files from individual sample processing are required for expression quantification
```

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

# Count TEProf2 per-sample results
find output -type d -name "TEProf2" | wc -l

# Check if TEProf2 aggregation completed
ls -l output/TEProf2_aggregated/filter_combined_candidates.tsv 2>/dev/null

# Check for errors in logs
grep -i "error" logs/*.err
grep -i "failed" logs/*.err
```

### TEProf2 Pipeline Stages

TEProf2 runs in two distinct phases:

1. **Per-sample processing** (runs for each sample):
   - Completed as part of `process_samples.sh` or `process_array.sh`
   - Creates individual sample outputs in `output/{path}/TEProf2/{sample}/`
   - Monitor with: `find output -name "*_annotated_filtered_test_all" | wc -l`

2. **Cross-sample aggregation** (runs once after all samples):
   - Automatically executed after all samples in `process_samples.sh`
   - For array jobs, run manually after all jobs complete:
     ```bash
     cd scripts
     ./run_teprof2_aggregation.sh
     ```
   - Creates aggregated results in `output/TEProf2_aggregated/`
   - Key output files:
     - `filter_combined_candidates.tsv`: All TE-gene transcripts
     - `read_filtered_candidates.tsv`: Filtered candidates
     - `reference_merged_candidates.gtf`: Merged reference with candidates
     - `table_tpm_cand`: Expression levels

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
  - TEProf2 per-sample: 1-3 hours depending on coverage
  
- TEProf2 aggregation (once, after all samples): 2-6 hours depending on sample count
  
- Total with array job (20 concurrent): ~10-20 hours + aggregation time
- Sequential processing: ~400-1400 hours + aggregation time

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

### TEProf2 per-sample outputs:
```
output/{path}/TEProf2/{sample}/
├── {sample}.Aligned.sortedByCoord.out.bam
├── {sample}.Aligned.sortedByCoord.out.bam.bai
├── {sample}.stringtie.gtf
├── {sample}.stringtie.gtf_annotated_filtered_test_all
├── {sample}.stringtie.gtf_annotated_filtered_test_all_annotation_tpm_processed
└── {sample}.stringtie.gtf_annotated_filtered_test_all_annotation_tpm_processed_filtered
```

### TEProf2 aggregated outputs (cross-sample analysis):
```
output/TEProf2_aggregated/
├── arguments.txt                        # Reference paths
├── filter_combined_candidates.tsv       # All TE-gene transcripts
├── initial_candidate_list.tsv           # Summary of unique candidates
├── read_filtered_candidates.tsv         # Candidates passing read filters
├── candidate_transcripts.gff3           # Detected transcripts (GFF3)
├── candidate_transcripts.gtf            # Detected transcripts (GTF)
├── reference_merged_candidates.gtf      # Merged with reference annotation
├── reference_merged_candidates.gff3     # Merged reference (GFF3)
├── reference_merged_candidates.gff3_annotated_filtered_test_all  # Final annotations
├── table_frac_tot_cand                  # Expression fraction per candidate
├── table_tpm_cand                       # TPM values per candidate
├── table_i_all                          # Intron coverage across samples
└── filterreadstats/                     # Read statistics per candidate
```

## Support

For tool-specific issues:
- JET: https://github.com/junseokpark/JET_identification_pipeline
- TEProf2: https://github.com/junseokpark/TEProf2Paper

## Version Info

- Pipeline version: 1.0
- Last updated: 2025
