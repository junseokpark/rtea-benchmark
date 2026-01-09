# How to Use the Refactored Pipeline

## Quick Start

### 1. Setup
```bash
cd scripts
./pipeline.sh setup
```

### 2. Configure (Optional)
Copy and edit the batch configuration:
```bash
cp batch-config.sh.template batch-config.sh
# Edit batch-config.sh to set your defaults
```

### 3. Run Full Workflow (Recommended)
```bash
# Process L1 samples at 30x coverage, 20 samples per batch
./pipeline.sh run --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists --partition compute
```

### 4. Or Run Step-by-Step

**Generate Sample Lists:**
```bash
./generate_sample_list.sh --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists
```

**Submit Jobs:**
```bash
./process_array.sh --sample_list_dir sample_lists/L1/30x --partition compute
```

**Check Status:**
```bash
./check_status.sh --sample_list_dir sample_lists/L1/30x
```

## Key Features

✅ **No Job Arrays** - Uses independent SLURM jobs instead of job arrays
✅ **Flexible Filtering** - Filter by TE type, coverage, or process all samples
✅ **Configurable Batching** - Set batch size to control samples per job
✅ **Partition Support** - Specify SLURM partition at submission time
✅ **Job Tracking** - Automatic tracking in `submitted_jobs.tsv` files
✅ **Status Monitoring** - Query SLURM for real-time job status
✅ **Robust Error Handling** - Safe execution with validation and helpful errors
✅ **Reproducible** - Deterministic file ordering and naming

## Example Workflows

### Example 1: All Samples with Small Batches
```bash
./pipeline.sh run --batch_size 10 --outdir sample_lists --partition compute
```

### Example 2: Specific TE Type, All Coverages
```bash
./pipeline.sh run --te_type AluY --batch_size 15 --outdir sample_lists
```

### Example 3: Custom Data Directory
```bash
./generate_sample_list.sh \
  --te_type L1 \
  --coverage 50x \
  --batch_size 25 \
  --outdir /path/to/output \
  --data_home /path/to/data

./process_array.sh \
  --sample_list_dir /path/to/output/L1/50x \
  --partition normal
```

## Job Tracking

All submitted jobs are tracked in `<sample_list_dir>/submitted_jobs.tsv`:

```tsv
list_file	job_id	submit_time
samples_0001.list	123456	2026-01-09 10:00:00
samples_0002.list	123457	2026-01-09 10:00:01
samples_0003.list	123458	2026-01-09 10:00:02
```

This file:
- Prevents duplicate submissions
- Enables status checking
- Provides submission history
- Can be manually edited to remove/add jobs

## Monitoring Jobs

### Check Status Periodically
```bash
# Basic status check
./check_status.sh --sample_list_dir sample_lists/L1/30x

# Detailed status report
./job_status_report.sh --sample_list_dir sample_lists/L1/30x

# Watch status (updates every 60 seconds)
watch -n 60 './job_status_report.sh --sample_list_dir sample_lists/L1/30x'
```

### Check with SLURM Commands
```bash
# View your jobs
squeue -u $USER

# View job details
sacct -j <job_id> --format=JobID,JobName,State,ExitCode,Elapsed

# View job logs
tail -f sample_lists/L1/30x/logs/samples_0001_*.out
```

## Resubmitting Jobs

### Resubmit All Jobs
```bash
./process_array.sh --sample_list_dir sample_lists/L1/30x --force
```

### Resubmit Specific Batches
Remove the corresponding entries from `submitted_jobs.tsv`, then:
```bash
./process_array.sh --sample_list_dir sample_lists/L1/30x
```

## Output Structure

```
sample_lists/
└── L1/                          # TE type
    └── 30x/                     # Coverage
        ├── samples_0001.list    # Batch 1: 20 samples
        ├── samples_0002.list    # Batch 2: 20 samples
        ├── samples_0003.list    # Batch 3: remaining samples
        ├── submitted_jobs.tsv   # Job tracking
        └── logs/                # Job logs
            ├── samples_0001_123456.out
            ├── samples_0001_123456.err
            ├── samples_0002_123457.out
            ├── samples_0002_123457.err
            └── ...
```

## Batch Size Recommendations

- **Small batches (5-10 samples)**: Fast start, more jobs, easier debugging
- **Medium batches (15-25 samples)**: Good balance (recommended)
- **Large batches (30+ samples)**: Fewer jobs, less scheduler overhead

## Troubleshooting

### No samples found
- Verify `DATA_HOME` in `config.sh`
- Check data directory structure matches expected format
- Verify TE type and coverage values

### Jobs not submitting
- Check sbatch is available: `which sbatch`
- Verify partition name exists: `sinfo`
- Review batch-config.sh settings

### Jobs failing
- Check logs: `ls -lh sample_lists/L1/30x/logs/`
- View errors: `cat sample_lists/L1/30x/logs/samples_0001_*.err`
- Check SLURM status: `sacct -j <job_id> --format=JobID,State,ExitCode`

## Documentation

- **BATCH_PROCESSING_GUIDE.md** - Comprehensive usage guide
- **REFACTORING_SUMMARY_NEW.md** - Technical documentation
- **All scripts** - Built-in help: `./script_name.sh --help`

## Migration from Old System

1. Update scripts: `git pull`
2. Run setup: `./pipeline.sh setup`
3. Update config: Copy and edit `batch-config.sh.template`
4. Generate new lists: `./generate_sample_list.sh --batch_size 20 ...`
5. Submit with new method: `./process_array.sh --sample_list_dir ...`

Old `sample_list.txt` files still work with legacy mode:
```bash
./check_status.sh --sample_list sample_list.txt
```

## Summary of Changes

- ❌ **Removed**: Job arrays (`--array=1-N`)
- ❌ **Removed**: Hardcoded TE type/coverage arrays
- ❌ **Removed**: Manual array size calculations
- ✅ **Added**: Independent job submission
- ✅ **Added**: Flexible CLI interfaces
- ✅ **Added**: Job tracking system
- ✅ **Added**: SLURM status monitoring
- ✅ **Added**: Comprehensive documentation

## Getting Help

1. Check documentation: `BATCH_PROCESSING_GUIDE.md`
2. Use help flag: `./script_name.sh --help`
3. Verify configuration: `./pipeline.sh check-config`
4. Review logs in the logs directory

---

**Note**: All changes maintain backward compatibility where possible. Legacy workflows using `sample_list.txt` still function with the `--sample_list` option.
