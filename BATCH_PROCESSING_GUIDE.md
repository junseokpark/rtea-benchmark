# SLURM Batch Processing - Usage Guide

This guide explains how to use the refactored pipeline scripts for batch processing without SLURM job arrays.

## Overview

The refactored pipeline supports:
- Splitting samples into batch files
- Submitting independent SLURM jobs (no job arrays)
- Job tracking via TSV files
- Status monitoring via SLURM queries
- Configurable partitions and resource allocation

## Quick Start

### 1. Setup

Make scripts executable:
```bash
cd scripts
./pipeline.sh setup
```

### 2. Generate Sample Lists

Generate batch lists for specific TE type and coverage:
```bash
./generate_sample_list.sh \
  --te_type L1 \
  --coverage 30x \
  --batch_size 20 \
  --outdir sample_lists
```

This creates:
- `sample_lists/L1/30x/samples_0001.list`
- `sample_lists/L1/30x/samples_0002.list`
- ... (one file per batch)

### 3. Submit Jobs

Submit all batch jobs to SLURM:
```bash
./process_array.sh \
  --sample_list_dir sample_lists/L1/30x \
  --partition compute
```

This creates `sample_lists/L1/30x/submitted_jobs.tsv` tracking all submitted jobs.

### 4. Check Status

Query SLURM for job status:
```bash
./check_status.sh --sample_list_dir sample_lists/L1/30x
```

Or use the dedicated status report:
```bash
./job_status_report.sh --sample_list_dir sample_lists/L1/30x
```

### 5. Full Workflow (One Command)

Run the entire workflow:
```bash
./pipeline.sh run \
  --te_type L1 \
  --coverage 30x \
  --batch_size 20 \
  --outdir sample_lists \
  --partition compute
```

## Command Reference

### generate_sample_list.sh

Generates sample lists and splits them into batch files.

**Required Options:**
- `--batch_size N` - Number of samples per batch file

**Optional Options:**
- `--te_type TYPE` - Filter by TE type (L1, AluY, L1HS, LTR5, SVA_F)
- `--coverage COV` - Filter by coverage (5x, 10x, 30x, 50x, 100x, 200x)
- `--outdir DIR` - Output directory (default: sample_lists)
- `--data_home DIR` - Data directory path

**Examples:**
```bash
# All samples, 10 per batch
./generate_sample_list.sh --batch_size 10 --outdir sample_lists

# Specific TE type at all coverages
./generate_sample_list.sh --te_type AluY --batch_size 15

# Specific TE type and coverage
./generate_sample_list.sh --te_type L1 --coverage 30x --batch_size 20
```

### process_array.sh

Submits independent SLURM jobs for each batch file.

**Required Options:**
- `--sample_list_dir DIR` - Directory containing *.list files

**Optional Options:**
- `--partition NAME` - SLURM partition
- `--job_name PREFIX` - Job name prefix (default: TE_batch)
- `--force` - Resubmit all jobs (ignore tracking file)

**Examples:**
```bash
# Basic submission
./process_array.sh --sample_list_dir sample_lists/L1/30x

# With specific partition
./process_array.sh --sample_list_dir sample_lists/L1/30x --partition compute

# Force resubmission
./process_array.sh --sample_list_dir sample_lists/L1/30x --force
```

### job_status_report.sh

Queries SLURM for job status and generates report.

**Options:**
- `--job_list FILE` - Path to submitted_jobs.tsv
- `--sample_list_dir DIR` - Directory containing submitted_jobs.tsv

**Exit Codes:**
- 0 - All jobs completed successfully
- 1 - Some jobs failed
- 2 - Jobs still running/pending

**Examples:**
```bash
# Using sample list directory
./job_status_report.sh --sample_list_dir sample_lists/L1/30x

# Using job list file directly
./job_status_report.sh --job_list sample_lists/L1/30x/submitted_jobs.tsv
```

### check_status.sh

Unified status checking (supports both job tracking and legacy mode).

**Options:**
- `--sample_list_dir DIR` - Check SLURM job status
- `--job_list FILE` - Check SLURM job status (TSV file)
- `--sample_list FILE` - Check output files (legacy mode)

**Examples:**
```bash
# Check SLURM job status
./check_status.sh --sample_list_dir sample_lists/L1/30x

# Legacy mode: check output files
./check_status.sh --sample_list sample_list.txt
```

### pipeline.sh

Master orchestration script.

**Commands:**
- `setup` - Make scripts executable
- `list` - Generate sample lists
- `submit` - Submit batch jobs
- `status` - Check job status
- `run` - Full workflow (generate + submit)
- `check-config` - Verify configuration
- `clean-logs` - Remove log files
- `summary` - Processing summary (legacy)

**Examples:**
```bash
# Full workflow
./pipeline.sh run --te_type L1 --coverage 30x --batch_size 20 --partition compute

# Just generate lists
./pipeline.sh list --te_type L1 --coverage 30x --batch_size 20

# Just submit jobs
./pipeline.sh submit --sample_list_dir sample_lists/L1/30x --partition compute

# Check status
./pipeline.sh status --sample_list_dir sample_lists/L1/30x
```

## Configuration

### batch-config.sh

Copy and customize the template:
```bash
cp batch-config.sh.template batch-config.sh
```

Edit `batch-config.sh` to set defaults:
```bash
export BATCH_TIME="24:00:00"
export BATCH_MEM="32G"
export BATCH_CPUS="8"
export SLURM_PARTITION="compute"
```

These settings are used as defaults but can be overridden via command-line arguments.

## Job Tracking

All submitted jobs are tracked in `<sample_list_dir>/submitted_jobs.tsv`:

```
list_file	job_id	submit_time
samples_0001.list	123456	2026-01-09 10:00:00
samples_0002.list	123457	2026-01-09 10:00:01
...
```

This file is:
- Created automatically by `process_array.sh`
- Used by `job_status_report.sh` and `check_status.sh` for status queries
- Updated with each submission
- Used to prevent duplicate submissions (unless `--force` is used)

## Output Structure

```
sample_lists/
└── L1/
    └── 30x/
        ├── samples_0001.list          # Batch 1 samples
        ├── samples_0002.list          # Batch 2 samples
        ├── ...
        ├── submitted_jobs.tsv         # Job tracking
        └── logs/
            ├── samples_0001_12345.out # Job output
            ├── samples_0001_12345.err # Job errors
            └── ...
```

## Workflow Examples

### Example 1: Process All L1 Samples at 30x Coverage

```bash
# Generate lists (20 samples per batch)
./generate_sample_list.sh \
  --te_type L1 \
  --coverage 30x \
  --batch_size 20 \
  --outdir sample_lists

# Submit to compute partition
./process_array.sh \
  --sample_list_dir sample_lists/L1/30x \
  --partition compute

# Monitor status
watch -n 60 './job_status_report.sh --sample_list_dir sample_lists/L1/30x'
```

### Example 2: Process All Samples (All TE Types, All Coverages)

```bash
# Generate lists
./generate_sample_list.sh \
  --batch_size 10 \
  --outdir sample_lists

# Submit jobs
./process_array.sh \
  --sample_list_dir sample_lists/all_samples \
  --partition compute
```

### Example 3: One-Command Workflow

```bash
./pipeline.sh run \
  --te_type AluY \
  --coverage 50x \
  --batch_size 25 \
  --outdir sample_lists \
  --partition normal
```

## Troubleshooting

### No samples found

Check that:
- `DATA_HOME` is correct in config.sh
- Data directory structure matches expected format
- TE type and coverage values match actual directories

### Jobs not submitting

Check that:
- sbatch is available (`which sbatch`)
- partition name is correct
- batch-config.sh has valid settings

### Jobs failing

Check logs in `<sample_list_dir>/logs/`:
```bash
ls -lh sample_lists/L1/30x/logs/
cat sample_lists/L1/30x/logs/samples_0001_12345.err
```

### Resubmit failed jobs

Use the `--force` flag:
```bash
./process_array.sh --sample_list_dir sample_lists/L1/30x --force
```

Or manually remove entries from `submitted_jobs.tsv` and resubmit.

## Key Differences from Job Arrays

| Feature | Old (Job Arrays) | New (Independent Jobs) |
|---------|-----------------|----------------------|
| Submission | Single sbatch with array | One sbatch per batch |
| Job IDs | Single ID with array indices | Multiple independent IDs |
| Tracking | Array task tracking | TSV file tracking |
| Resubmission | Resubmit entire array | Resubmit individual batches |
| Flexibility | Fixed array size | Dynamic batch sizes |
| Dependencies | Array-level only | Job-level possible |

## Migration Notes

If migrating from the old job array system:

1. The old `sample_list.txt` format is still compatible
2. Use `generate_sample_list.sh` to split into batches
3. Job tracking is now explicit via `submitted_jobs.tsv`
4. Status checking uses SLURM queries instead of file checking
5. No need to calculate array sizes manually

## Performance Considerations

**Batch Size Selection:**
- Smaller batches (5-10 samples): More granular control, faster start
- Medium batches (15-25 samples): Good balance for most cases
- Larger batches (30+ samples): Fewer jobs, reduced scheduler overhead

**Recommended:**
- For quick testing: 5-10 samples per batch
- For production: 15-20 samples per batch
- For very large datasets: 25-30 samples per batch

## Support

For issues or questions:
1. Check this documentation
2. Review script help: `./script_name.sh --help`
3. Check SLURM logs in the logs directory
4. Verify configuration with `./pipeline.sh check-config`
