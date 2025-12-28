# Test Script for scripts/function.sh

## Overview

`test_functions.sh` is a comprehensive test script that verifies all functions in `scripts/function.sh` by executing them with a test configuration file and capturing their output for analysis.

## Features

- ✅ **Configuration Loading**: Automatically loads config.sh or falls back to config_template.sh
- ✅ **Environment Validation**: Checks for Singularity, configuration variables, and required files
- ✅ **Test Environment Setup**: Creates isolated test data directories in /tmp
- ✅ **Function Testing**: Tests run_jet_step1, run_jet_step2, and run_teprof2
- ✅ **Output Capture**: Captures both stdout and stderr for each function
- ✅ **Error Detection**: Identifies Singularity errors and execution failures
- ✅ **Detailed Logging**: Timestamped logs with all execution details
- ✅ **Summary Report**: Clear pass/fail status for each function

## Usage

### Basic Usage

```bash
cd scripts
./test_functions.sh
```

### With Custom Configuration

```bash
# First, create your configuration file
cp config_template.sh config.sh
# Edit config.sh with your actual paths

# Then run the test
./test_functions.sh
```

## Output

The script generates:

1. **Console Output**: Color-coded test results with:
   - Configuration loading status
   - Environment validation warnings
   - Test execution progress
   - Function output summaries
   - Final test report

2. **Log Directory**: `/tmp/function_test_<timestamp>/` containing:
   - `test_functions.log` - Main test log with timestamps
   - `<function_name>_output.log` - Standard output from each function
   - `<function_name>_error.log` - Standard error from each function
   - `test_data/` - Test data directory structure

## Test Functions

### 1. run_jet_step1
Tests JET Step 1 - STAR Alignment
- Checks JET2 singularity image availability
- Executes alignment with test metadata
- Verifies execution without Singularity errors

### 2. run_jet_step2  
Tests JET Step 2 - R Analysis
- Checks JET2 singularity image availability
- Executes R analysis with test parameters
- Validates R processing completion

### 3. run_teprof2
Tests TEProf2 Analysis
- Checks TEProf2 singularity image availability
- Executes TEProf2 analysis workflow
- Confirms analysis execution status

## Exit Codes

- `0`: All tests passed
- `1`: One or more tests failed (expected when Singularity is not available)

## Example Output

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Function Test Suite for scripts/function.sh
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[2025-12-28 04:30:19] Starting test suite

...

╔══════════════════════════════════════════════════════════════════════╗
║ Testing: run_jet_step1
╚══════════════════════════════════════════════════════════════════════╝

ℹ Testing JET Step 1 - STAR Alignment
ℹ Checking prerequisites...
ℹ Executing function: run_jet_step1
ℹ Capturing output to: /tmp/function_test_20251228_043019/run_jet_step1_output.log
ℹ Analyzing execution results...
✗ run_jet_step1: Function failed (exit code: 1)

...

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Test Summary Report
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Test Execution Summary:
  Total Tests Run: 3
  Passed: 0
  Failed: 3

Detailed Results:
  ✗ run_teprof2: Function failed (exit code: 1)
  ✗ run_jet_step2: Function failed (exit code: 1)
  ✗ run_jet_step1: Function failed (exit code: 1)

ℹ Test logs saved to: /tmp/function_test_20251228_043019
ℹ Main log file: /tmp/function_test_20251228_043019/test_functions.log
```

## Expected Behavior

### In Testing Environment (without Singularity)
- Functions will execute but fail with "command not found" errors
- Error detection will identify the missing Singularity installation
- All test logs will be captured for review
- Exit code will be 1 (expected failure)

### In Production Environment (with Singularity)
- Functions will execute with actual Singularity containers
- Success/failure depends on data availability and configuration
- Singularity errors will be detected and reported
- Full output captured for debugging

## Configuration Requirements

The test script requires these configuration variables:
- `JET2` - Path to JET singularity image
- `TEProf2` - Path to TEProf2 singularity image
- Reference file paths (fastaFile, gtfGeneFile, etc.)
- Tool paths within containers (samtoolsBinDir, starBinDir, etc.)

See `config_template.sh` for complete list.

## Troubleshooting

### Tests Fail with "Singularity not found"
**Solution**: This is expected in environments without Singularity. The test script still validates that functions execute and captures their output correctly.

### Tests Fail with "Container not found"
**Solution**: Update paths in config.sh to point to actual container image files.

### Tests Fail with "Reference files not found"
**Solution**: Update reference file paths in config.sh or ensure files exist at specified locations.

## Integration with CI/CD

This script can be integrated into continuous integration pipelines:

```bash
# In CI pipeline
cd scripts
./test_functions.sh
# Check exit code for pass/fail status
```

## Maintenance

When adding new functions to `function.sh`:
1. Add a corresponding `test_<function_name>()` function in `test_functions.sh`
2. Call the new test function in the `main()` function's test execution section
3. Update this README with the new function's test description

## Related Files

- `scripts/function.sh` - Functions being tested
- `scripts/config_template.sh` - Configuration template
- `scripts/config.sh` - User configuration (create from template)

## Support

For issues or questions about the test script, please refer to the main repository README or open an issue.
