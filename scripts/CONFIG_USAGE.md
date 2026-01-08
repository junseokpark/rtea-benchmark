# Configuration and SBATCH Settings Usage Guide

This guide explains how to configure and customize the TE analysis pipeline scripts (`process_samples.sh` and `process_array.sh`).

## Overview

The refactored scripts now use a centralized configuration system with two types of configuration files:

1. **Pipeline Configuration** (`config.sh` / `config_template.sh`) - Required settings for JET and TEProf2
2. **Batch Job Configuration** (`batch-config.sh`) - Optional SLURM SBATCH settings

## Pipeline Configuration

### Setup

1. Copy the template to create your configuration:
   ```bash
   cd scripts
   cp config_template.sh config.sh
   ```

2. Edit `config.sh` and update all paths marked as `REQUIRED`:
   - `DATA_HOME` - Directory containing your input FASTQ files
   - `OUTPUT_BASE` - Directory for pipeline outputs
   - `REF_DIR` - Directory containing reference files
   - `JET2` - Path to JET singularity image
   - `TEProf2` - Path to TEProf2 singularity image
   - Reference files (FASTA, GTF, GFF, etc.)

3. The scripts will automatically:
   - Load `config.sh` if it exists
   - Fall back to `config_template.sh` if `config.sh` is not found (with a warning)
   - Exit with an error if neither file is found

### Required Variables

The following variables must be set in your `config.sh`:

- `DATA_HOME` - Input data directory
- `OUTPUT_BASE` - Output directory
- `REF_DIR` - Reference files directory
- `JET2` - JET singularity image path
- `TEProf2` - TEProf2 singularity image path
- `refDir` - Reference directory for JET (can be same as `REF_DIR`)

Additional variables are documented in `config_template.sh`.

## Batch Job Configuration (SLURM)

### SBATCH Settings Overview

The scripts support three ways to customize SLURM job settings (in order of precedence):

1. **Command-line override** (highest priority)
2. **batch-config.sh file**
3. **Built-in defaults** (lowest priority)

### Method 1: Command-Line Override (Recommended for one-time changes)

Override SBATCH settings at job submission time:

```bash
# Override time limit and memory
sbatch --export=ALL,BATCH_TIME=48:00:00,BATCH_MEM=64G process_array.sh

# Override multiple settings
sbatch --export=ALL,BATCH_TIME=72:00:00,BATCH_MEM=128G,BATCH_CPUS=16 process_samples.sh
```

**Available environment variables:**
- `BATCH_TIME` - Job time limit (e.g., `24:00:00`, `48:00:00`, `3-00:00:00`)
- `BATCH_MEM` - Memory allocation (e.g., `32G`, `64G`, `128G`)
- `BATCH_CPUS` - CPUs per task (e.g., `8`, `16`, `32`)
- `BATCH_ARRAY` - Array job range (e.g., `1-400`, `1-100`) - **process_array.sh only**
- `BATCH_ARRAY_THROTTLE` - Max concurrent array jobs (e.g., `20`, `50`) - **process_array.sh only**
- `BATCH_LOG_DIR` - Log directory (e.g., `logs`, `custom_logs`)

### Method 2: batch-config.sh File (Recommended for persistent settings)

1. Copy the template:
   ```bash
   cd scripts
   cp batch-config.sh.template batch-config.sh
   ```

2. Edit `batch-config.sh` to set your preferred defaults:
   ```bash
   # Example: Increase resources for large datasets
   export BATCH_TIME="72:00:00"
   export BATCH_MEM="128G"
   export BATCH_CPUS="16"
   export BATCH_ARRAY="1-400"
   export BATCH_ARRAY_THROTTLE="30"
   ```

3. Submit jobs normally:
   ```bash
   sbatch process_array.sh
   ```

The scripts will automatically source `batch-config.sh` if it exists.

### Method 3: Built-in Defaults

If neither command-line overrides nor `batch-config.sh` are provided, the scripts use these defaults:

**process_array.sh:**
- Time: `24:00:00`
- Memory: `32G`
- CPUs: `8`
- Array: `1-400%20` (max 20 concurrent)

**process_samples.sh:**
- Time: `48:00:00`
- Memory: `32G`
- CPUs: `8`

## Usage Examples

### Example 1: First-time Setup

```bash
cd scripts

# Setup pipeline configuration
cp config_template.sh config.sh
vim config.sh  # Edit with your paths

# Use default SBATCH settings
sbatch process_array.sh
```

### Example 2: Custom SBATCH Settings for Specific Job

```bash
# Run with more memory and longer time limit
sbatch --export=ALL,BATCH_TIME=96:00:00,BATCH_MEM=256G process_samples.sh
```

### Example 3: Persistent Custom SBATCH Settings

```bash
cd scripts

# Create batch config
cp batch-config.sh.template batch-config.sh

# Edit batch-config.sh
cat >> batch-config.sh << 'EOF'
export BATCH_TIME="48:00:00"
export BATCH_MEM="64G"
export BATCH_CPUS="16"
export BATCH_PARTITION="large-mem"
export BATCH_ACCOUNT="myproject"
EOF

# All future submissions will use these settings
sbatch process_array.sh
```

### Example 4: Combining batch-config.sh and Command-Line

```bash
# batch-config.sh sets: BATCH_MEM=64G, BATCH_CPUS=16
# Override only the time limit at submission
sbatch --export=ALL,BATCH_TIME=96:00:00 process_array.sh

# Result: Time=96h (from command-line), Memory=64G (from batch-config.sh), CPUs=16 (from batch-config.sh)
```

## Important Notes

### SLURM #SBATCH Directive Limitations

The `#SBATCH` directives in the script headers are parsed **before** the script runs, so environment variables cannot be used directly in these lines. The current implementation:

1. Sets **safe defaults** in `#SBATCH` directives
2. Documents how to override via environment variables or `batch-config.sh`
3. Users can override at submission time using `sbatch --export=ALL,VAR=value`

To fully customize `#SBATCH` directives with variables, you would need to either:
- Use the command-line override method (recommended)
- Edit the `#SBATCH` lines in the scripts directly (not recommended)

### Sample List for process_array.sh

The `process_array.sh` script requires a `sample_list.txt` file with format:
```
sample_name fq1_path fq2_path rel_path
```

Generate this file using `generate_sample_list.sh` or create it manually.

### Error Handling

Both scripts now use `set -euo pipefail` for strict error handling:
- Scripts will exit immediately on any error
- Unset variables will cause errors
- Pipeline failures will propagate

## Troubleshooting

### Config file not found

```
ERROR: Neither config.sh nor config_template.sh found
```

**Solution:** Make sure you're running the script from the correct directory or that `config_template.sh` exists in the `scripts/` directory.

### Missing required variables

```
ERROR: Required configuration variables not set:
  - DATA_HOME
  - JET2
```

**Solution:** Edit your `config.sh` and set all required variables listed in the error message.

### batch-config.sh not loaded

```
INFO: batch-config.sh not found, using default SBATCH settings
```

**Note:** This is informational only. The scripts will work fine with default settings. To customize, create `batch-config.sh` from the template.

## Migration from Old Scripts

If you're migrating from the old hardcoded configuration scripts:

1. Your old configuration values should be moved to `config.sh`
2. Any custom SBATCH settings should be moved to `batch-config.sh` or passed via command-line
3. The scripts maintain the same functionality and calling conventions
4. No changes needed to `function.sh` or other supporting scripts

## Additional Resources

- See `config_template.sh` for complete list of configuration variables
- See `batch-config.sh.template` for SBATCH customization options
- See `test_functions.sh` for examples of configuration loading
