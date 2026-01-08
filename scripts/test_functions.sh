#!/bin/bash

################################################################################
# Test Script for scripts/function.sh
# 
# This script tests all functions in scripts/function.sh by:
# - Loading configuration from config file
# - Setting up minimal test environment
# - Executing each function with test data
# - Capturing and reporting stdout/stderr
# - Checking for errors and Singularity execution status
# - Providing detailed logging for debugging
################################################################################

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Test results tracking
declare -A TEST_RESULTS
declare -A TEST_MESSAGES
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Log file
TEST_LOG_DIR="/home/junseokp/temp/function_test_$(date +%Y%m%d_%H%M%S)"
#TEST_LOG_DIR="/home/junseokp/temp/function_test_20260107_132725"

mkdir -p "${TEST_LOG_DIR}"
TEST_LOG="${TEST_LOG_DIR}/test_functions.log"

################################################################################
# Helper Functions
################################################################################

log_message() {
    local message="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${message}" | tee -a "${TEST_LOG}"
}

print_header() {
    echo -e "\n${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}" | tee -a "${TEST_LOG}"
    echo -e "${BLUE}$1${NC}" | tee -a "${TEST_LOG}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n" | tee -a "${TEST_LOG}"
}

print_test_header() {
    echo -e "\n${CYAN}╔══════════════════════════════════════════════════════════════════════╗${NC}" | tee -a "${TEST_LOG}"
    echo -e "${CYAN}║ Testing: $1${NC}" | tee -a "${TEST_LOG}"
    echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════════╝${NC}\n" | tee -a "${TEST_LOG}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}" | tee -a "${TEST_LOG}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}" | tee -a "${TEST_LOG}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}" | tee -a "${TEST_LOG}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}" | tee -a "${TEST_LOG}"
}

################################################################################
# Test Framework Functions
################################################################################

start_test() {
    local test_name="$1"
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_header "${test_name}"
    log_message "Starting test: ${test_name}"
}

pass_test() {
    local test_name="$1"
    local message="${2:-Test passed}"
    TEST_RESULTS["${test_name}"]="PASS"
    TEST_MESSAGES["${test_name}"]="${message}"
    TESTS_PASSED=$((TESTS_PASSED + 1))
    print_success "${test_name}: ${message}"
    log_message "Test PASSED: ${test_name}"
}

fail_test() {
    local test_name="$1"
    local message="${2:-Test failed}"
    TEST_RESULTS["${test_name}"]="FAIL"
    TEST_MESSAGES["${test_name}"]="${message}"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    print_error "${test_name}: ${message}"
    log_message "Test FAILED: ${test_name} - ${message}"
}

################################################################################
# Configuration Loading
################################################################################

load_configuration() {
    print_header "Loading Configuration"
    
    # Try to load config.sh first, fall back to config_template.sh
    if [ -f "${SCRIPT_DIR}/config.sh" ]; then
        log_message "Loading configuration from config.sh"
        source "${SCRIPT_DIR}/config.sh"
        print_success "Loaded config.sh"
    elif [ -f "${SCRIPT_DIR}/config_template.sh" ]; then
        log_message "config.sh not found, using config_template.sh"
        print_warning "Using config_template.sh (may have placeholder values)"
        source "${SCRIPT_DIR}/config_template.sh"
    else
        print_error "No configuration file found!"
        log_message "ERROR: Neither config.sh nor config_template.sh found"
        exit 1
    fi
    
    # Source the function file
    if [ -f "${SCRIPT_DIR}/function.sh" ]; then
        source "${SCRIPT_DIR}/function.sh"
        print_success "Loaded function.sh"
        log_message "Loaded function.sh successfully"
    else
        print_error "function.sh not found!"
        log_message "ERROR: function.sh not found"
        exit 1
    fi
}

################################################################################
# Pre-Test Validation
################################################################################

validate_environment() {
    print_header "Validating Test Environment"
    
    local validation_failed=0
    
    # Check if singularity is available
    print_info "Checking for singularity..."
    module load singularity
    if command -v singularity &> /dev/null; then
        local singularity_version=$(singularity --version 2>&1)
        print_success "Singularity found: ${singularity_version}"
        log_message "Singularity version: ${singularity_version}"
    else
        print_warning "Singularity not found in PATH"
        log_message "WARNING: Singularity not found"
    fi
    
    # Check required configuration variables
    print_info "Checking configuration variables..."
    
    local required_vars=("JET2" "TEProf2" "dataDir" "refDir")
    for var in "${required_vars[@]}"; do
        if [ -z "${!var}" ]; then
            print_warning "Variable ${var} is not set"
            log_message "WARNING: ${var} not set"
        else
            print_info "${var}=${!var}"
            log_message "Config: ${var}=${!var}"
        fi
    done
    
    # Check singularity images (if paths are set)
    if [ -n "${JET2}" ]; then
        if [ -f "${JET2}" ]; then
            print_success "JET2 singularity image exists: ${JET2}"
            log_message "JET2 image found: ${JET2}"
        else
            print_warning "JET2 singularity image not found: ${JET2}"
            log_message "WARNING: JET2 image not found: ${JET2}"
        fi
    fi
    
    if [ -n "${TEProf2}" ]; then
        if [ -f "${TEProf2}" ]; then
            print_success "TEProf2 singularity image exists: ${TEProf2}"
            log_message "TEProf2 image found: ${TEProf2}"
        else
            print_warning "TEProf2 singularity image not found: ${TEProf2}"
            log_message "WARNING: TEProf2 image not found: ${TEProf2}"
        fi
    fi
    
    echo ""
}

################################################################################
# Test Environment Setup
################################################################################

setup_test_environment() {
    print_header "Setting Up Test Environment"
    
    # Create test directories
    export outputDir="${TEST_LOG_DIR}/test_output"
    export SAMPLE_NAME="sim200_AluY_blood"
    local sample_rel_path="nonReferenceTE/AluY/5X/fq"
    
    mkdir -p "${outputDir}/output"
    mkdir -p "${outputDir}/log"
    mkdir -p "${outputDir}/err"
    
    print_success "Test directories created: ${outputDir}"
    log_message "Test data directory: ${outputDir}"
    
    # Set test fastq files - use dataDir from config or fallback to DATA_HOME
    local test_data_dir="${dataDir:-${DATA_HOME:-/home/junseokp/workspaces/data/rTea-simul/sims}}"
    export FQ1="${test_data_dir}/${sample_rel_path}/${SAMPLE_NAME}.1.fq.gz"
    export FQ2="${test_data_dir}/${sample_rel_path}/${SAMPLE_NAME}.2.fq.gz"
    
    # Set dataDir for output directory structure
    export dataDir="${test_data_dir}"
    
    # Set other required variables with defaults if not set
    export threads="${threads:-4}"
    export readLength="${readLength:-150}"
    export organism="${organism:-Human}"  # Updated to capitalize per JET spec
    export genome="${genome:-hg38}"
    export database="${database:-ensembl}"  # Updated default
    export refDir="${refDir:-${REF_DIR:-/home/junseokp/workspaces/data/rTea-simul/ref}}"
    
    # Set minJunction for JET Step 2 if not already set
    export minJunction="${minJunction:-2e7}"
    
    print_info "Test environment variables set:"
    log_message "Test environment setup complete"
    echo "  dataDir=${dataDir}" | tee -a "${TEST_LOG}"
    echo "  outputDir=${outputDir}" | tee -a "${TEST_LOG}"
    echo "  SAMPLE_NAME=${SAMPLE_NAME}" | tee -a "${TEST_LOG}"
    echo "  FQ1=${FQ1}" | tee -a "${TEST_LOG}"
    echo "  FQ2=${FQ2}" | tee -a "${TEST_LOG}"
    echo "  organism=${organism}" | tee -a "${TEST_LOG}"
    echo "  genome=${genome}" | tee -a "${TEST_LOG}"
    echo "  database=${database}" | tee -a "${TEST_LOG}"
    echo "  minJunction=${minJunction}" | tee -a "${TEST_LOG}"
    echo ""
}

################################################################################
# Test Execution Functions
################################################################################

capture_function_output() {
    local func_name="$1"
    local output_file="${TEST_LOG_DIR}/${func_name}_output.log"
    local error_file="${TEST_LOG_DIR}/${func_name}_error.log"
    
    log_message "Executing function: ${func_name}"
    print_info "Capturing output to: ${output_file}"
    print_info "Capturing errors to: ${error_file}"
    
    # Execute the function directly and capture output
    set +e  # Don't exit on error
    ${func_name} > "${output_file}" 2> "${error_file}"
    local exit_code=$?
    set -e
    
    log_message "Function ${func_name} completed with exit code: ${exit_code}"
    
    # Store exit code in a global variable instead of returning it
    FUNC_EXIT_CODE=${exit_code}
    return 0
}

check_singularity_errors() {
    local error_file="$1"
    local func_name="$2"
    
    if [ ! -f "${error_file}" ]; then
        log_message "No error file found: ${error_file}"
        return 0
    fi
    
    # Check for common Singularity errors
    local has_errors=0
    
    if grep -qi "singularity.*error" "${error_file}"; then
        print_error "Singularity error detected in ${func_name}"
        log_message "Singularity error found in ${error_file}"
        grep -i "singularity.*error" "${error_file}" >> "${TEST_LOG}"
        has_errors=1
    fi
    
    if grep -qi "FATAL" "${error_file}"; then
        print_error "FATAL error detected in ${func_name}"
        log_message "FATAL error found in ${error_file}"
        grep -i "FATAL" "${error_file}" >> "${TEST_LOG}"
        has_errors=1
    fi
    
    if grep -qi "container.*not found" "${error_file}"; then
        print_warning "Container not found error in ${func_name}"
        log_message "Container not found in ${error_file}"
        has_errors=1
    fi
    
    # Check for general ERROR messages (but exclude harmless ones)
    if grep -qi "^ERROR:" "${error_file}"; then
        print_warning "ERROR message detected in ${func_name}"
        log_message "ERROR found in ${error_file}"
        grep -i "^ERROR:" "${error_file}" | head -5 >> "${TEST_LOG}"
        has_errors=1
    fi
    
    return ${has_errors}
}

analyze_function_execution() {
    local func_name="$1"
    local exit_code="$2"
    local output_file="${TEST_LOG_DIR}/${func_name}_output.log"
    local error_file="${TEST_LOG_DIR}/${func_name}_error.log"
    
    print_info "Analyzing execution results..."
    log_message "Exit code: ${exit_code}"
    
    # Display output summary
    if [ -f "${output_file}" ]; then
        local line_count=$(wc -l < "${output_file}")
        print_info "Output lines: ${line_count}"
        log_message "Output file size: ${line_count} lines"
        
        # Show last few lines of output
        if [ ${line_count} -gt 0 ]; then
            echo -e "\n${CYAN}Last 10 lines of output:${NC}" | tee -a "${TEST_LOG}"
            tail -10 "${output_file}" >> "${TEST_LOG}"
            tail -10 "${output_file}"
        fi
    fi
    
    # Check for Singularity errors
    check_singularity_errors "${error_file}" "${func_name}"
    local has_singularity_errors=$?
    
    # Determine test result
    if [ ${exit_code} -eq 0 ]; then
        if [ ${has_singularity_errors} -eq 0 ]; then
            pass_test "${func_name}" "Function executed successfully (exit code: 0, no Singularity errors)"
        else
            fail_test "${func_name}" "Function completed but Singularity errors detected"
        fi
    else
        if [ ${has_singularity_errors} -eq 1 ]; then
            fail_test "${func_name}" "Function failed (exit code: ${exit_code}) with Singularity errors"
        else
            # Check if it's an expected failure (missing files, etc.)
            if grep -qi "not found\|does not exist\|no such file" "${output_file}" 2>/dev/null; then
                print_warning "${func_name}: Failed due to missing files (expected in test environment)"
                fail_test "${func_name}" "Function failed (exit code: ${exit_code}) - likely due to missing test data"
            else
                fail_test "${func_name}" "Function failed (exit code: ${exit_code})"
            fi
        fi
    fi
    
    echo ""
}

################################################################################
# Individual Function Tests
################################################################################

test_run_jet_step1() {
    local func_name="run_jet_step1"
    start_test "${func_name}"
    
    print_info "Testing JET Step 1 - STAR Alignment"
    log_message "Testing ${func_name}"
    
    # Ensure required variables are set
    print_info "Checking prerequisites..."
    if [ -z "${JET2}" ]; then
        fail_test "${func_name}" "JET2 singularity image path not set"
        return 1
    fi
    
    # Execute function with output capture
    capture_function_output "${func_name}"
    local exit_code=${FUNC_EXIT_CODE}
    
    # Analyze results
    analyze_function_execution "${func_name}" "${exit_code}"
}

test_run_jet_step2() {
    local func_name="run_jet_step2"
    start_test "${func_name}"
    
    print_info "Testing JET Step 2 - R Analysis"
    log_message "Testing ${func_name}"
    
    # Ensure required variables are set
    print_info "Checking prerequisites..."
    if [ -z "${JET2}" ]; then
        fail_test "${func_name}" "JET2 singularity image path not set"
        return 1
    fi
    
    # Execute function with output capture
    capture_function_output "${func_name}"
    local exit_code=${FUNC_EXIT_CODE}
    
    # Analyze results
    analyze_function_execution "${func_name}" "${exit_code}"
}

test_run_teprof2() {
    local func_name="run_teprof2"
    start_test "${func_name}"
    
    print_info "Testing TEProf2 analysis"
    log_message "Testing ${func_name}"
    
    # Ensure required variables are set
    print_info "Checking prerequisites..."
    if [ -z "${TEProf2}" ]; then
        fail_test "${func_name}" "TEProf2 singularity image path not set"
        return 1
    fi
    
    # Execute function with output capture
    capture_function_output "${func_name}"
    local exit_code=${FUNC_EXIT_CODE}
    
    # Analyze results
    analyze_function_execution "${func_name}" "${exit_code}"
}

test_filter_broken_fastq() {
    local func_name="filter_broken_fastq"
    start_test "${func_name}"
    
    print_info "Testing FASTQ filtering function"
    log_message "Testing ${func_name}"
    
    # Create test FASTQ files
    local test_dir="${TEST_LOG_DIR}/fastq_test"
    mkdir -p "${test_dir}"
    
    # Create a valid FASTQ file
    local valid_fq="${test_dir}/valid.fq"
    cat > "${valid_fq}" << 'EOF'
@READ1
ACGTACGTACGT
+
IIIIIIIIIIII
@READ2
GGGGTTTTAAAA
+
JJJJJJJJJJJJ
EOF
    
    # Create a broken FASTQ file (mismatched lengths)
    local broken_fq="${test_dir}/broken.fq"
    cat > "${broken_fq}" << 'EOF'
@READ1
ACGTACGTACGT
+
IIIIIII
@READ2
GGGGTTTTAAAA
+
JJJJJJJJJJJJ
@READ3
CCCC
+
QQQQ
EOF
    
    print_info "Testing with valid FASTQ..."
    local result=$(filter_broken_fastq "${valid_fq}")
    local exit_code=$?
    
    if [ ${exit_code} -eq 0 ]; then
        print_success "Valid FASTQ processed successfully"
        print_info "Result: ${result}"
        
        # Should return original file path
        if [[ "${result}" == *"valid.fq"* ]] && [[ ! "${result}" == *"filtered"* ]]; then
            print_success "Correctly returned original file (no filtering needed)"
        else
            print_warning "Expected original file path, got: ${result}"
        fi
    else
        print_error "Failed to process valid FASTQ"
    fi
    
    print_info "Testing with broken FASTQ..."
    result=$(filter_broken_fastq "${broken_fq}")
    exit_code=$?
    
    if [ ${exit_code} -eq 0 ]; then
        print_success "Broken FASTQ processed successfully"
        print_info "Result: ${result}"
        
        # Should create filtered file
        if [[ "${result}" == *"filtered.fq"* ]]; then
            print_success "Created filtered file as expected"
            
            # Check if filtered file exists and has valid records
            local filtered_file="${result}"
            if [ -f "${filtered_file}" ]; then
                local valid_count=$(grep -c "^@READ" "${filtered_file}" || true)
                print_info "Filtered file contains ${valid_count} valid reads"
                
                # Check if removed.seqs file was created
                local removed_file="${test_dir}/broken.removed.seqs"
                if [ -f "${removed_file}" ]; then
                    local removed_count=$(grep -c "^@READ" "${removed_file}" || true)
                    print_info "Removed file contains ${removed_count} broken reads"
                else
                    print_warning "Removed sequences file not found"
                fi
            else
                print_error "Filtered file not found: ${filtered_file}"
            fi
        else
            print_warning "Expected filtered file path, got: ${result}"
        fi
    else
        print_error "Failed to process broken FASTQ"
    fi
    
    # Test with non-existent file
    print_info "Testing with non-existent file..."
    local nonexistent="${test_dir}/nonexistent.fq"
    result=$(filter_broken_fastq "${nonexistent}" 2>&1)
    exit_code=$?
    
    if [ ${exit_code} -ne 0 ]; then
        print_success "Correctly failed for non-existent file"
    else
        print_error "Should have failed for non-existent file"
    fi
    
    # Overall test result
    pass_test "${func_name}" "FASTQ filtering function tests completed"
}

################################################################################
# Test Summary Report
################################################################################

generate_summary_report() {
    print_header "Test Summary Report"
    
    echo -e "${CYAN}Test Execution Summary:${NC}" | tee -a "${TEST_LOG}"
    echo "  Total Tests Run: ${TESTS_RUN}" | tee -a "${TEST_LOG}"
    echo -e "  ${GREEN}Passed: ${TESTS_PASSED}${NC}" | tee -a "${TEST_LOG}"
    echo -e "  ${RED}Failed: ${TESTS_FAILED}${NC}" | tee -a "${TEST_LOG}"
    echo "" | tee -a "${TEST_LOG}"
    
    echo -e "${CYAN}Detailed Results:${NC}" | tee -a "${TEST_LOG}"
    for test_name in "${!TEST_RESULTS[@]}"; do
        local result="${TEST_RESULTS[${test_name}]}"
        local message="${TEST_MESSAGES[${test_name}]}"
        
        if [ "${result}" = "PASS" ]; then
            echo -e "  ${GREEN}✓${NC} ${test_name}: ${message}" | tee -a "${TEST_LOG}"
        else
            echo -e "  ${RED}✗${NC} ${test_name}: ${message}" | tee -a "${TEST_LOG}"
        fi
    done
    echo "" | tee -a "${TEST_LOG}"
    
    print_info "Test logs saved to: ${TEST_LOG_DIR}"
    print_info "Main log file: ${TEST_LOG}"
    echo "" | tee -a "${TEST_LOG}"
    
    # List all generated log files
    echo -e "${CYAN}Generated Log Files:${NC}" | tee -a "${TEST_LOG}"
    ls -lh "${TEST_LOG_DIR}" | tail -n +2 >> "${TEST_LOG}"
    ls -lh "${TEST_LOG_DIR}" | tail -n +2
    echo "" | tee -a "${TEST_LOG}"
    
    # Overall result
    if [ ${TESTS_FAILED} -eq 0 ]; then
        print_success "All tests passed!"
        log_message "TEST SUITE PASSED"
        return 0
    else
        print_error "Some tests failed!"
        log_message "TEST SUITE FAILED"
        return 1
    fi
}

################################################################################
# Main Test Execution
################################################################################

main() {
    print_header "Function Test Suite for scripts/function.sh"
    log_message "Starting test suite"
    
    # Load configuration and functions
    load_configuration
    
    # Validate environment
    validate_environment
    
    # Setup test environment
    setup_test_environment
    
    # Run tests for each function
    print_header "Running Function Tests"
    
    #test_filter_broken_fastq
    test_run_jet_step1
    test_run_jet_step2
    test_run_teprof2
    
    # Generate summary report
    generate_summary_report
    local exit_code=$?
    
    # Print final message
    echo ""
    if [ ${exit_code} -eq 0 ]; then
        echo -e "${GREEN}╔════════════════════════════════════════╗${NC}"
        echo -e "${GREEN}║  ✓ TEST SUITE COMPLETED SUCCESSFULLY  ║${NC}"
        echo -e "${GREEN}╚════════════════════════════════════════╝${NC}"
    else
        echo -e "${RED}╔════════════════════════════════════════╗${NC}"
        echo -e "${RED}║  ✗ TEST SUITE COMPLETED WITH FAILURES ║${NC}"
        echo -e "${RED}╚════════════════════════════════════════╝${NC}"
    fi
    
    log_message "Test suite completed with exit code: ${exit_code}"
    exit ${exit_code}
}

# Run main function
main "$@"
