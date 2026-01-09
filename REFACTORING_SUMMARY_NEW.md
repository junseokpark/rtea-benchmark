# Pipeline Refactoring Summary

## Overview

This document summarizes the major refactoring of the TE Analysis Pipeline scripts to support SLURM batch job submission **without using job arrays**.

## Key Changes

### 1. New Architecture: Independent Jobs Instead of Job Arrays

**Before:**
- Single `sbatch` command with `--array=1-N`
- Fixed array size determined at submission
- Single job ID with array task indices
- Sample processing driven by `SLURM_ARRAY_TASK_ID`

**After:**
- Multiple independent `sbatch` commands
- One job per batch file
- Multiple independent job IDs
- Sample processing driven by batch file content
- Job tracking via TSV file

### 2. Modified Scripts

#### `generate_sample_list.sh` - Completely Refactored

**Old Behavior:**
- Hardcoded paths and arrays
- Generated single `sample_list.txt` file
- No filtering options
- No batching support

**New Behavior:**
- Accepts CLI parameters (--te_type, --coverage, --batch_size, --outdir, --data_home)
- Filters samples based on criteria
- Splits samples into multiple batch files
- Creates deterministic, sortable filenames (samples_0001.list, samples_0002.list, etc.)
- Organizes output by TE type and coverage

**Usage:**
```bash
./generate_sample_list.sh --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists
```

#### `process_array.sh` - Completely Rewritten

**Old Behavior:**
- SLURM batch script with job array directives
- Processed samples based on `SLURM_ARRAY_TASK_ID`
- Hardcoded SBATCH settings

**New Behavior:**
- Shell script that submits independent jobs
- Reads directory of `.list` files
- Submits one job per batch file using `process_batch_worker.sh`
- Tracks submitted jobs in `submitted_jobs.tsv`
- Supports --partition, --force, --job_name options
- Robust job ID parsing from various sbatch formats

**Usage:**
```bash
./process_array.sh --sample_list_dir sample_lists/L1/30x --partition compute
```

#### `process_batch_worker.sh` - NEW

**Purpose:**
- Worker script that processes all samples in a single batch file
- Submitted as an independent SLURM job
- Replaces the array job processing logic

**Features:**
- Takes sample list file as argument
- Loads configuration
- Processes each sample sequentially
- Reports failures and successes
- Uses existing JET and TEProf2 functions

#### `job_status_report.sh` - NEW (replaces calculate_array_size.sh)

**Old Script (calculate_array_size.sh):**
- Calculated array size needed for samples
- Provided guidance for setting --array parameter

**New Script (job_status_report.sh):**
- Reads `submitted_jobs.tsv`
- Queries SLURM for each job's status (squeue and sacct)
- Reports counts by status (running, pending, completed, failed, etc.)
- Exit codes indicate overall status
- Supports both --job_list and --sample_list_dir options

**Usage:**
```bash
./job_status_report.sh --sample_list_dir sample_lists/L1/30x
```

#### `check_status.sh` - Enhanced

**Old Behavior:**
- Checked output files for sample completion
- Generated processing status reports
- No SLURM integration

**New Behavior:**
- Supports both new job tracking mode and legacy file checking mode
- New mode: Delegates to `job_status_report.sh` for SLURM queries
- Legacy mode: Maintains backward compatibility for file-based checking
- Accepts --job_list, --sample_list_dir, or --sample_list

**Usage:**
```bash
# New mode: Check SLURM jobs
./check_status.sh --sample_list_dir sample_lists/L1/30x

# Legacy mode: Check output files
./check_status.sh --sample_list sample_list.txt
```

#### `pipeline.sh` - Updated

**New Features:**
- New `run` command for full workflow orchestration
- New `list` command with parameter pass-through
- New `submit` command with parameter pass-through
- Enhanced `status` command supporting new modes
- Removed obsolete commands (submit-seq, validate, clean, resubmit)

**Usage:**
```bash
# Full workflow
./pipeline.sh run --te_type L1 --coverage 30x --batch_size 20 --partition compute

# Individual steps
./pipeline.sh list --te_type L1 --coverage 30x --batch_size 20
./pipeline.sh submit --sample_list_dir sample_lists/L1/30x
./pipeline.sh status --sample_list_dir sample_lists/L1/30x
```

#### `batch-config.sh.template` - Updated

**New Variables:**
- `SLURM_PARTITION` - Partition name configuration
- Removed job array-specific settings (BATCH_ARRAY, BATCH_ARRAY_THROTTLE, SAMPLES_PER_JOB)
- Updated documentation for new architecture

### 3. New Features

#### Job Tracking System

- TSV file format: `<sample_list_dir>/submitted_jobs.tsv`
- Columns: `list_file`, `job_id`, `submit_time`
- Prevents duplicate submissions (unless --force used)
- Enables status queries without hardcoded job IDs
- Supports historical tracking

Example:
```tsv
list_file	job_id	submit_time
samples_0001.list	123456	2026-01-09 10:00:00
samples_0002.list	123457	2026-01-09 10:00:01
```

#### Flexible Filtering

- Filter by TE type: L1, AluY, L1HS, LTR5, SVA_F
- Filter by coverage: 5x, 10x, 30x, 50x, 100x, 200x
- Process all samples or specific subsets
- Configurable batch sizes

#### Robust Error Handling

All scripts now include:
- `set -euo pipefail` for safe execution
- Input validation with helpful error messages
- Argument parsing with `--help` support
- Reproducible file discovery with `LC_ALL=C sort`

#### Partition Support

- Configurable via `--partition` flag
- Default from `batch-config.sh` (SLURM_PARTITION)
- Falls back to SLURM default if not specified

### 4. Workflow Comparison

#### Old Workflow
```bash
# 1. Generate single sample list
./generate_sample_list.sh

# 2. Calculate array size
./calculate_array_size.sh 20

# 3. Manually edit process_array.sh to set --array=1-N

# 4. Submit array job
sbatch process_array.sh

# 5. Check status (file-based)
./check_status.sh
```

#### New Workflow
```bash
# One-command approach
./pipeline.sh run --te_type L1 --coverage 30x --batch_size 20 --partition compute

# OR step-by-step approach

# 1. Generate batch lists
./generate_sample_list.sh --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists

# 2. Submit jobs (no manual editing needed)
./process_array.sh --sample_list_dir sample_lists/L1/30x --partition compute

# 3. Check status (SLURM-based)
./check_status.sh --sample_list_dir sample_lists/L1/30x
```

### 5. Benefits of New Architecture

1. **No Job Arrays Required**
   - Works on clusters with job array restrictions
   - Simpler job management

2. **Flexible Resubmission**
   - Resubmit individual failed batches
   - No need to resubmit entire array

3. **Dynamic Batch Sizes**
   - Change batch size without recalculating arrays
   - Adapt to workload requirements

4. **Better Tracking**
   - Explicit job tracking via TSV files
   - Historical record of submissions
   - Easier debugging

5. **Improved UX**
   - Clear command-line interface
   - Helpful error messages
   - Comprehensive help documentation

6. **Configuration Flexibility**
   - Filter samples by TE type and coverage
   - Specify partition at submission time
   - Override defaults easily

7. **Reproducibility**
   - Deterministic file ordering (LC_ALL=C sort)
   - Consistent batch file naming
   - Tracked job submissions

## File Structure Changes

### Old Structure
```
scripts/
├── generate_sample_list.sh     # Hardcoded, single output
├── process_array.sh            # Job array script
├── calculate_array_size.sh     # Array size calculator
├── check_status.sh             # File-based checking
└── pipeline.sh                 # Basic orchestration
```

### New Structure
```
scripts/
├── generate_sample_list.sh     # CLI-based, batched output
├── process_array.sh            # Job submission script
├── process_batch_worker.sh     # Worker script (NEW)
├── job_status_report.sh        # SLURM status queries (NEW)
├── check_status.sh             # Dual-mode: SLURM + legacy
├── pipeline.sh                 # Enhanced orchestration
└── batch-config.sh.template    # Updated configuration
```

### Output Structure

Old:
```
sample_list.txt       # Single file
logs/                 # Job array logs
```

New:
```
sample_lists/
└── L1/
    └── 30x/
        ├── samples_0001.list
        ├── samples_0002.list
        ├── ...
        ├── submitted_jobs.tsv      # NEW: Job tracking
        └── logs/
            ├── samples_0001_12345.out
            ├── samples_0001_12345.err
            └── ...
```

## Migration Guide

### For Existing Users

1. **Update Scripts**
   - Pull latest changes
   - Run `./pipeline.sh setup` to make scripts executable

2. **Configuration**
   - Copy `batch-config.sh.template` to `batch-config.sh`
   - Set `SLURM_PARTITION` if needed
   - Remove old job array settings

3. **Generate New Batch Lists**
   ```bash
   ./generate_sample_list.sh --batch_size 20 --outdir sample_lists
   ```

4. **Submit Using New Method**
   ```bash
   ./process_array.sh --sample_list_dir sample_lists/all_samples
   ```

5. **Check Status**
   ```bash
   ./check_status.sh --sample_list_dir sample_lists/all_samples
   ```

### Backward Compatibility

- Legacy `sample_list.txt` format still works
- `check_status.sh --sample_list` provides legacy file checking
- Old output directory structure is compatible

## Testing

All scripts have been validated:

✅ Help messages work correctly
✅ Error handling for missing arguments
✅ Sample list generation with test data
✅ Batch file splitting (tested with 5 samples, batch size 2)
✅ Job tracking file format
✅ Status reporting (with mock job data)

## Documentation

New comprehensive documentation:
- `BATCH_PROCESSING_GUIDE.md` - Complete usage guide
- Updated help text in all scripts
- In-script comments and examples

## Summary

This refactoring successfully transforms the TE Analysis Pipeline from a job array-based system to an independent job submission system, providing:

- Greater flexibility in job management
- Better tracking and debugging capabilities
- Improved user experience with CLI interfaces
- Enhanced configurability
- Maintained backward compatibility where possible

All requirements from the original issue have been implemented and tested.
