#!/bin/bash
set -euo pipefail

# ============================================
# Master Control Script for TE Analysis Pipeline
# ============================================
# Orchestrates the batch processing workflow
#
# Usage: ./pipeline.sh [command] [OPTIONS]
# ============================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper functions
print_header() {
    echo -e "\n${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

# Command functions
cmd_help() {
    cat << EOF
TE Analysis Pipeline - Master Control Script

USAGE:
    ./pipeline.sh [command] [OPTIONS]

COMMANDS:
    help            Show this help message
    
    setup           Initial setup (make scripts executable)
    
    list            Generate sample lists with batching
                    Options: --te_type, --coverage, --batch_size, --outdir, --data_home
    
    submit          Submit batch jobs for processing
                    Options: --sample_list_dir, --partition
    
    status          Check job/processing status
                    Options: --sample_list_dir or --job_list
    
    run             Full workflow: generate lists, submit jobs
                    Options: --te_type, --coverage, --batch_size, --outdir, --partition, --data_home
    
    check-config    Verify configuration and requirements
    
    clean-logs      Remove log files
    
    summary         Generate processing summary (legacy mode)

EXAMPLES:
    # Full workflow for L1 at 30x coverage
    ./pipeline.sh run --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists --partition compute

    # Generate sample lists only
    ./pipeline.sh list --te_type L1 --coverage 30x --batch_size 20 --outdir sample_lists

    # Submit jobs for existing batch files
    ./pipeline.sh submit --sample_list_dir sample_lists/L1/30x --partition compute

    # Check status
    ./pipeline.sh status --sample_list_dir sample_lists/L1/30x

For detailed information, see: README.md or QUICKSTART.txt
EOF
}

cmd_setup() {
    print_header "Setting up pipeline scripts"
    
    chmod +x generate_sample_list.sh
    chmod +x process_array.sh
    chmod +x process_batch_worker.sh
    chmod +x job_status_report.sh
    chmod +x check_status.sh
    chmod +x resubmit_failed.sh
    
    mkdir -p logs
    
    print_success "Scripts made executable"
    print_success "Logs directory created"
    print_info "Setup complete!"
}

cmd_check_config() {
    print_header "Checking configuration"
    
    # Load config to check
    if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
        source "${SCRIPT_DIR}/config.sh"
    elif [[ -f "${SCRIPT_DIR}/config_template.sh" ]]; then
        source "${SCRIPT_DIR}/config_template.sh"
    else
        print_error "No configuration file found"
        return 1
    fi
    
    errors=0
    
    # Check required variables
    required_vars=("DATA_HOME" "OUTPUT_BASE" "JET2" "TEProf2" "refDir")
    for var in "${required_vars[@]}"; do
        if [[ -z "${!var:-}" ]]; then
            print_error "Variable ${var} is not set"
            errors=$((errors + 1))
        else
            print_success "${var} is set"
        fi
    done
    
    # Check paths exist
    if [[ -n "${DATA_HOME:-}" && -d "${DATA_HOME}" ]]; then
        print_success "Data directory exists: ${DATA_HOME}"
    else
        print_error "Data directory not found: ${DATA_HOME:-NOT_SET}"
        errors=$((errors + 1))
    fi
    
    if [[ -n "${JET2:-}" && -f "${JET2}" ]]; then
        print_success "JET container found: ${JET2}"
    else
        print_error "JET container not found: ${JET2:-NOT_SET}"
        errors=$((errors + 1))
    fi
    
    if [[ -n "${TEProf2:-}" && -f "${TEProf2}" ]]; then
        print_success "TEProf2 container found: ${TEProf2}"
    else
        print_error "TEProf2 container not found: ${TEProf2:-NOT_SET}"
        errors=$((errors + 1))
    fi
    
    echo ""
    if [[ $errors -eq 0 ]]; then
        print_success "Configuration check passed!"
        return 0
    else
        print_error "Configuration check failed with $errors error(s)"
        return 1
    fi
}

cmd_list() {
    print_header "Generating sample lists"
    
    if [[ ! -x "generate_sample_list.sh" ]]; then
        print_error "generate_sample_list.sh not executable. Run: ./pipeline.sh setup"
        return 1
    fi
    
    # Pass all remaining arguments to generate_sample_list.sh
    ./generate_sample_list.sh "$@"
}

cmd_submit() {
    print_header "Submitting batch jobs"
    
    if [[ ! -x "process_array.sh" ]]; then
        print_error "process_array.sh not executable. Run: ./pipeline.sh setup"
        return 1
    fi
    
    # Pass all remaining arguments to process_array.sh
    ./process_array.sh "$@"
}

cmd_status() {
    print_header "Checking status"
    
    if [[ ! -x "check_status.sh" ]]; then
        print_error "check_status.sh not executable. Run: ./pipeline.sh setup"
        return 1
    fi
    
    # Pass all remaining arguments to check_status.sh
    ./check_status.sh "$@"
}

cmd_run() {
    print_header "Running full workflow"
    
    # Parse arguments for workflow
    local te_type="" coverage="" batch_size="" outdir="" partition="" data_home=""
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --te_type)
                te_type="$2"
                shift 2
                ;;
            --coverage)
                coverage="$2"
                shift 2
                ;;
            --batch_size)
                batch_size="$2"
                shift 2
                ;;
            --outdir)
                outdir="$2"
                shift 2
                ;;
            --partition)
                partition="$2"
                shift 2
                ;;
            --data_home)
                data_home="$2"
                shift 2
                ;;
            *)
                print_error "Unknown option for run command: $1"
                return 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [[ -z "${batch_size}" ]]; then
        print_error "--batch_size is required for run command"
        return 1
    fi
    
    if [[ -z "${outdir}" ]]; then
        outdir="sample_lists"
        print_info "Using default outdir: ${outdir}"
    fi
    
    # Build generate_sample_list.sh command
    print_info "Step 1: Generating sample lists..."
    local gen_cmd="./generate_sample_list.sh --batch_size ${batch_size} --outdir ${outdir}"
    [[ -n "${te_type}" ]] && gen_cmd="${gen_cmd} --te_type ${te_type}"
    [[ -n "${coverage}" ]] && gen_cmd="${gen_cmd} --coverage ${coverage}"
    [[ -n "${data_home}" ]] && gen_cmd="${gen_cmd} --data_home ${data_home}"
    
    eval "${gen_cmd}"
    
    if [[ $? -ne 0 ]]; then
        print_error "Failed to generate sample lists"
        return 1
    fi
    
    # Determine sample_list_dir from the output
    local sample_list_dir
    if [[ -n "${te_type}" && -n "${coverage}" ]]; then
        sample_list_dir="${outdir}/${te_type}/${coverage}"
    elif [[ -n "${te_type}" ]]; then
        sample_list_dir="${outdir}/${te_type}/all_coverages"
    elif [[ -n "${coverage}" ]]; then
        sample_list_dir="${outdir}/all_te_types/${coverage}"
    else
        sample_list_dir="${outdir}/all_samples"
    fi
    
    print_success "Sample lists created in: ${sample_list_dir}"
    echo ""
    
    # Submit jobs
    print_info "Step 2: Submitting jobs..."
    local submit_cmd="./process_array.sh --sample_list_dir ${sample_list_dir}"
    [[ -n "${partition}" ]] && submit_cmd="${submit_cmd} --partition ${partition}"
    
    eval "${submit_cmd}"
    
    if [[ $? -ne 0 ]]; then
        print_error "Failed to submit jobs"
        return 1
    fi
    
    print_success "Jobs submitted successfully"
    echo ""
    print_info "Monitor jobs with:"
    print_info "  squeue -u \${USER}"
    print_info "  ./pipeline.sh status --sample_list_dir ${sample_list_dir}"
}

cmd_clean_logs() {
    print_header "Cleaning log files"
    
    if [[ -d "logs" ]]; then
        log_count=$(find logs -type f 2>/dev/null | wc -l)
        print_info "Found $log_count log files"
        
        read -p "Delete all log files? [y/N] " -n 1 -r
        echo
        
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            rm -f logs/*
            print_success "Log files deleted"
        else
            print_info "Cancelled"
        fi
    else
        print_info "No logs directory found"
    fi
}

cmd_summary() {
    print_header "Generating processing summary (legacy mode)"
    
    # Load config
    if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
        source "${SCRIPT_DIR}/config.sh"
    elif [[ -f "${SCRIPT_DIR}/config_template.sh" ]]; then
        source "${SCRIPT_DIR}/config_template.sh"
    fi
    
    OUTPUT_BASE="${OUTPUT_BASE:-/home/junseokp/workspaces/data/rTea-simul/output}"
    
    if [[ ! -d "$OUTPUT_BASE" ]]; then
        print_error "Output directory not found: $OUTPUT_BASE"
        return 1
    fi
    
    print_info "Counting results..."
    
    jet_bams=$(find "$OUTPUT_BASE" -name "*.sorted.bam" 2>/dev/null | wc -l)
    teprof2_dirs=$(find "$OUTPUT_BASE" -type d -name "TEProf2" 2>/dev/null | wc -l)
    
    echo ""
    print_info "JET BAM files: $jet_bams"
    print_info "TEProf2 output directories: $teprof2_dirs"
    
    disk_usage=$(du -sh "$OUTPUT_BASE" 2>/dev/null | cut -f1)
    print_info "Total disk usage: $disk_usage"
}

# Main script logic
main() {
    if [[ $# -eq 0 ]]; then
        cmd_help
        exit 0
    fi
    
    # Change to script directory
    cd "${SCRIPT_DIR}"
    
    command=$1
    shift
    
    case $command in
        help|--help|-h)
            cmd_help
            ;;
        setup)
            cmd_setup
            ;;
        check-config|config)
            cmd_check_config
            ;;
        list)
            cmd_list "$@"
            ;;
        submit)
            cmd_submit "$@"
            ;;
        status)
            cmd_status "$@"
            ;;
        run)
            cmd_run "$@"
            ;;
        clean-logs)
            cmd_clean_logs
            ;;
        summary)
            cmd_summary
            ;;
        *)
            print_error "Unknown command: $command"
            echo ""
            cmd_help
            exit 1
            ;;
    esac
}

main "$@"


cmd_setup() {
    print_header "Setting up pipeline scripts"
    
    chmod +x generate_sample_list.sh
    chmod +x process_array.sh
    chmod +x process_samples.sh
    chmod +x check_status.sh
    chmod +x resubmit_failed.sh
    
    mkdir -p logs
    
    print_success "Scripts made executable"
    print_success "Logs directory created"
    print_info "Setup complete!"
}

cmd_check_config() {
    print_header "Checking configuration"
    
    DATA_HOME=/home/junseokp/workspaces/data/rTea-simul
    REF=${DATA_HOME}/ref
    JET=/home/sasidharp/jet_docker/jet.sif
    TEProf2=/home/sasidharp/jet_docker/teprof2.sif
    
    errors=0
    
    # Check DATA_HOME
    if [ -d "$DATA_HOME" ]; then
        print_success "Data directory exists: $DATA_HOME"
    else
        print_error "Data directory not found: $DATA_HOME"
        errors=$((errors + 1))
    fi
    
    # Check REF directory
    if [ -d "$REF" ]; then
        print_success "Reference directory exists: $REF"
    else
        print_warning "Reference directory not found: $REF"
        print_info "You may need to create this directory and add reference files"
    fi
    
    # Check singularity containers
    if [ -f "$JET" ]; then
        print_success "JET container found: $JET"
    else
        print_error "JET container not found: $JET"
        errors=$((errors + 1))
    fi
    
    if [ -f "$TEProf2" ]; then
        print_success "TEProf2 container found: $TEProf2"
    else
        print_error "TEProf2 container not found: $TEProf2"
        errors=$((errors + 1))
    fi
    
    # Check for sample directories
    if [ -d "$DATA_HOME/nonReferenceTE" ]; then
        print_success "nonReferenceTE directory found"
    else
        print_error "nonReferenceTE directory not found"
        errors=$((errors + 1))
    fi
    
    if [ -d "$DATA_HOME/referenceTE" ]; then
        print_success "referenceTE directory found"
    else
        print_error "referenceTE directory not found"
        errors=$((errors + 1))
    fi
    
    # Check module availability
    if command -v module &> /dev/null; then
        print_success "Module system available"
    else
        print_warning "Module system not found (may not be needed)"
    fi
    
    echo ""
    if [ $errors -eq 0 ]; then
        print_success "Configuration check passed!"
        return 0
    else
        print_error "Configuration check failed with $errors error(s)"
        return 1
    fi
}

cmd_list() {
    print_header "Generating sample list"
    
    if [ ! -x "generate_sample_list.sh" ]; then
        print_error "generate_sample_list.sh not executable. Run: ./pipeline.sh setup"
        return 1
    fi
    
    ./generate_sample_list.sh
    
    if [ -f "sample_list.txt" ]; then
        sample_count=$(grep -v "^#" sample_list.txt | wc -l)
        print_success "Sample list created with $sample_count samples"
        print_info "File: sample_list.txt"
    else
        print_error "Failed to create sample list"
        return 1
    fi
}

cmd_submit() {
    print_header "Submitting array job"
    
    if [ ! -f "sample_list.txt" ]; then
        print_error "sample_list.txt not found. Run: ./pipeline.sh list"
        return 1
    fi
    
    sample_count=$(grep -v "^#" sample_list.txt | wc -l)
    print_info "Total samples to process: $sample_count"
    
    # Update array size in script
    sed -i "s/--array=1-[0-9]*%20/--array=1-${sample_count}%20/" process_array.sh
    
    print_info "Submitting job..."
    sbatch process_array.sh
    
    if [ $? -eq 0 ]; then
        print_success "Job submitted successfully"
        print_info "Monitor with: squeue -u \$USER"
        print_info "Check logs in: logs/"
    else
        print_error "Job submission failed"
        return 1
    fi
}

cmd_submit_seq() {
    print_header "Submitting sequential processing job"
    
    print_warning "Sequential processing will take much longer than array job"
    read -p "Continue? [y/N] " -n 1 -r
    echo
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Cancelled"
        return 0
    fi
    
    print_info "Submitting job..."
    sbatch process_samples.sh
    
    if [ $? -eq 0 ]; then
        print_success "Job submitted successfully"
    else
        print_error "Job submission failed"
        return 1
    fi
}

cmd_status() {
    print_header "Checking processing status"
    
    if [ ! -f "sample_list.txt" ]; then
        print_error "sample_list.txt not found. Run: ./pipeline.sh list"
        return 1
    fi
    
    ./check_status.sh
}

cmd_resubmit() {
    print_header "Resubmitting failed samples"
    
    if [ ! -f "sample_list.txt" ]; then
        print_error "sample_list.txt not found"
        return 1
    fi
    
    ./resubmit_failed.sh
    
    # Check if any resubmission script was created
    resubmit_script=$(ls -t resubmit_failed_*.sh 2>/dev/null | head -1)
    
    if [ -n "$resubmit_script" ]; then
        print_info "Resubmission script created: $resubmit_script"
        read -p "Submit now? [y/N] " -n 1 -r
        echo
        
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            sbatch $resubmit_script
            print_success "Resubmission job submitted"
        fi
    else
        print_success "No failed samples found!"
    fi
}

cmd_clean() {
    print_header "Cleaning output directory"
    
    OUTPUT_BASE=/home/junseokp/workspaces/data/rTea-simul/output
    
    if [ ! -d "$OUTPUT_BASE" ]; then
        print_info "Output directory does not exist: $OUTPUT_BASE"
        return 0
    fi
    
    print_warning "This will DELETE all output files in: $OUTPUT_BASE"
    print_warning "This action CANNOT be undone!"
    echo ""
    read -p "Are you sure? Type 'DELETE' to confirm: " confirm
    
    if [ "$confirm" = "DELETE" ]; then
        print_info "Removing output directory..."
        rm -rf "$OUTPUT_BASE"
        print_success "Output directory cleaned"
    else
        print_info "Cancelled"
    fi
}

cmd_clean_logs() {
    print_header "Cleaning log files"
    
    if [ -d "logs" ]; then
        log_count=$(find logs -type f | wc -l)
        print_info "Found $log_count log files"
        
        read -p "Delete all log files? [y/N] " -n 1 -r
        echo
        
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            rm -f logs/*
            print_success "Log files deleted"
        else
            print_info "Cancelled"
        fi
    else
        print_info "No logs directory found"
    fi
}

cmd_summary() {
    print_header "Generating processing summary"
    
    OUTPUT_BASE=/home/junseokp/workspaces/data/rTea-simul/output
    
    if [ ! -d "$OUTPUT_BASE" ]; then
        print_error "Output directory not found: $OUTPUT_BASE"
        return 1
    fi
    
    print_info "Counting results..."
    
    jet_bams=$(find $OUTPUT_BASE -name "*.sorted.bam" 2>/dev/null | wc -l)
    teprof2_dirs=$(find $OUTPUT_BASE -type d -name "TEProf2" 2>/dev/null | wc -l)
    
    echo ""
    print_info "JET BAM files: $jet_bams"
    print_info "TEProf2 output directories: $teprof2_dirs"
    
    disk_usage=$(du -sh $OUTPUT_BASE 2>/dev/null | cut -f1)
    print_info "Total disk usage: $disk_usage"
}

cmd_validate() {
    print_header "Validating outputs"
    
    print_info "Running status check..."
    ./check_status.sh
    
    print_info "Validating file integrity..."
    
    OUTPUT_BASE=/home/junseokp/workspaces/data/rTea-simul/output
    
    # Check for corrupted BAM files
    corrupted=0
    for bam in $(find $OUTPUT_BASE -name "*.sorted.bam" 2>/dev/null); do
        if ! samtools quickcheck "$bam" 2>/dev/null; then
            print_error "Corrupted BAM: $bam"
            corrupted=$((corrupted + 1))
        fi
    done
    
    if [ $corrupted -eq 0 ]; then
        print_success "All BAM files valid"
    else
        print_warning "Found $corrupted corrupted BAM files"
    fi
}

# Main script logic
main() {
    if [ $# -eq 0 ]; then
        cmd_help
        exit 0
    fi
    
    command=$1
    shift
    
    case $command in
        help|--help|-h)
            cmd_help
            ;;
        setup)
            cmd_setup
            ;;
        check-config|config)
            cmd_check_config
            ;;
        list)
            cmd_list
            ;;
        submit)
            cmd_submit
            ;;
        submit-seq)
            cmd_submit_seq
            ;;
        status)
            cmd_status
            ;;
        resubmit)
            cmd_resubmit
            ;;
        clean)
            cmd_clean
            ;;
        clean-logs)
            cmd_clean_logs
            ;;
        summary)
            cmd_summary
            ;;
        validate)
            cmd_validate
            ;;
        *)
            print_error "Unknown command: $command"
            echo ""
            cmd_help
            exit 1
            ;;
    esac
}

main "$@"
