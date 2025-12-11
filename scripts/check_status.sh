#!/bin/bash

# Check processing status for all samples
DATA_HOME=/home/junseokp/workspaces/data/rTea-simul
OUTPUT_BASE=${DATA_HOME}/output
SAMPLE_LIST="sample_list.txt"

REPORT_FILE="processing_status_$(date +%Y%m%d_%H%M%S).txt"

echo "TE Analysis Pipeline - Processing Status Report" > ${REPORT_FILE}
echo "Generated: $(date)" >> ${REPORT_FILE}
echo "=============================================" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

# Counters
total_samples=0
jet_step1_complete=0
jet_step1_failed=0
jet_step1_missing=0
jet_step2_complete=0
jet_step2_failed=0
jet_step2_missing=0
teprof2_complete=0
teprof2_failed=0
teprof2_missing=0

# Function to check if JET Step 1 completed successfully
check_jet_step1() {
    local sample_dir=$1
    
    # Check for STAR alignment output
    if [ -f "${sample_dir}/output/Aligned.sortedByCoord.out.bam" ] || \
       [ -f "${sample_dir}/output/"*"Aligned.sortedByCoord.out.bam" ]; then
        echo "COMPLETE"
    elif [ -d "${sample_dir}/output" ] && [ -f "${sample_dir}/log/step1_multisample_running_"*".log" ]; then
        echo "FAILED"
    else
        echo "NOT_STARTED"
    fi
}

# Function to check if JET Step 2 completed successfully
check_jet_step2() {
    local sample_dir=$1
    
    # Check for R analysis output
    if [ -d "${sample_dir}/output" ] && \
       [ -f "${sample_dir}/log/step2_multisample_running_"*".log" ] && \
       [ -f "${sample_dir}/output/"*"_results.txt" -o -f "${sample_dir}/output/"*"_TE_insertions.bed" ]; then
        echo "COMPLETE"
    elif [ -f "${sample_dir}/log/step2_multisample_running_"*".log" ]; then
        echo "FAILED"
    else
        echo "NOT_STARTED"
    fi
}

# Function to check if TEProf2 completed successfully
check_teprof2() {
    local sample_dir=$1
    
    # Check for expected output directory and files
    if [ -d "${sample_dir}/TEProf2" ] && \
       [ "$(ls -A ${sample_dir}/TEProf2 2>/dev/null)" ]; then
        echo "COMPLETE"
    elif [ -d "${sample_dir}/TEProf2" ]; then
        echo "FAILED"
    else
        echo "NOT_STARTED"
    fi
}

echo "Checking samples from: ${SAMPLE_LIST}" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

# Process each sample
while IFS=' ' read -r sample_name fq1 fq2 rel_path; do
    # Skip header
    if [[ "$sample_name" == "#"* ]]; then
        continue
    fi
    
    total_samples=$((total_samples + 1))
    sample_dir="${OUTPUT_BASE}/${rel_path}/${sample_name}"
    
    # Check JET Step 1 status
    jet_step1_status=$(check_jet_step1 "$sample_dir")
    case $jet_step1_status in
        "COMPLETE")
            jet_step1_complete=$((jet_step1_complete + 1))
            ;;
        "FAILED")
            jet_step1_failed=$((jet_step1_failed + 1))
            echo "JET STEP1 FAILED: ${sample_name} (${rel_path})" >> ${REPORT_FILE}
            ;;
        "NOT_STARTED")
            jet_step1_missing=$((jet_step1_missing + 1))
            ;;
    esac
    
    # Check JET Step 2 status
    jet_step2_status=$(check_jet_step2 "$sample_dir")
    case $jet_step2_status in
        "COMPLETE")
            jet_step2_complete=$((jet_step2_complete + 1))
            ;;
        "FAILED")
            jet_step2_failed=$((jet_step2_failed + 1))
            echo "JET STEP2 FAILED: ${sample_name} (${rel_path})" >> ${REPORT_FILE}
            ;;
        "NOT_STARTED")
            jet_step2_missing=$((jet_step2_missing + 1))
            ;;
    esac
    
    # Check TEProf2 status
    teprof2_status=$(check_teprof2 "$sample_dir")
    case $teprof2_status in
        "COMPLETE")
            teprof2_complete=$((teprof2_complete + 1))
            ;;
        "FAILED")
            teprof2_failed=$((teprof2_failed + 1))
            echo "TEProf2 FAILED: ${sample_name} (${rel_path})" >> ${REPORT_FILE}
            ;;
        "NOT_STARTED")
            teprof2_missing=$((teprof2_missing + 1))
            ;;
    esac
    
done < ${SAMPLE_LIST}

# Write summary
echo "" >> ${REPORT_FILE}
echo "=============================================" >> ${REPORT_FILE}
echo "SUMMARY" >> ${REPORT_FILE}
echo "=============================================" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}
echo "Total Samples: ${total_samples}" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}
echo "JET Step 1 Status:" >> ${REPORT_FILE}
echo "  Complete: ${jet_step1_complete} ($(awk "BEGIN {printf \"%.1f\", ${jet_step1_complete}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "  Failed: ${jet_step1_failed} ($(awk "BEGIN {printf \"%.1f\", ${jet_step1_failed}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "  Not Started: ${jet_step1_missing} ($(awk "BEGIN {printf \"%.1f\", ${jet_step1_missing}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}
echo "JET Step 2 Status:" >> ${REPORT_FILE}
echo "  Complete: ${jet_step2_complete} ($(awk "BEGIN {printf \"%.1f\", ${jet_step2_complete}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "  Failed: ${jet_step2_failed} ($(awk "BEGIN {printf \"%.1f\", ${jet_step2_failed}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "  Not Started: ${jet_step2_missing} ($(awk "BEGIN {printf \"%.1f\", ${jet_step2_missing}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}
echo "TEProf2 Status:" >> ${REPORT_FILE}
echo "  Complete: ${teprof2_complete} ($(awk "BEGIN {printf \"%.1f\", ${teprof2_complete}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "  Failed: ${teprof2_failed} ($(awk "BEGIN {printf \"%.1f\", ${teprof2_failed}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "  Not Started: ${teprof2_missing} ($(awk "BEGIN {printf \"%.1f\", ${teprof2_missing}/${total_samples}*100}")%)" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

# Breakdown by category
echo "=============================================" >> ${REPORT_FILE}
echo "BREAKDOWN BY CATEGORY" >> ${REPORT_FILE}
echo "=============================================" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

# nonReferenceTE breakdown
echo "nonReferenceTE:" >> ${REPORT_FILE}
for te_type in AluY L1HS LTR5 SVA_F; do
    count=$(grep -c "nonReferenceTE/${te_type}/" ${SAMPLE_LIST} 2>/dev/null || echo 0)
    if [ $count -gt 0 ]; then
        echo "  ${te_type}: ${count} samples" >> ${REPORT_FILE}
    fi
done
echo "" >> ${REPORT_FILE}

# referenceTE breakdown
echo "referenceTE/intron:" >> ${REPORT_FILE}
intron_count=$(grep -c "referenceTE/intron/" ${SAMPLE_LIST} 2>/dev/null || echo 0)
echo "  Total: ${intron_count} samples" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

echo "referenceTE/TSS:" >> ${REPORT_FILE}
tss_count=$(grep -c "referenceTE/TSS/" ${SAMPLE_LIST} 2>/dev/null || echo 0)
echo "  Total: ${tss_count} samples" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

# Coverage breakdown
echo "=============================================" >> ${REPORT_FILE}
echo "BREAKDOWN BY COVERAGE" >> ${REPORT_FILE}
echo "=============================================" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

for coverage in 5X 10X 50X 100X 200X; do
    count=$(grep -c "${coverage}/" ${SAMPLE_LIST} 2>/dev/null || echo 0)
    if [ $count -gt 0 ]; then
        echo "${coverage}: ${count} samples" >> ${REPORT_FILE}
    fi
done

# Display report
cat ${REPORT_FILE}

echo ""
echo "Report saved to: ${REPORT_FILE}"

# Check for failed samples
if [ $jet_step1_failed -gt 0 ] || [ $jet_step2_failed -gt 0 ] || [ $teprof2_failed -gt 0 ]; then
    echo ""
    echo "WARNING: Some samples failed processing!"
    echo "Review ${REPORT_FILE} for details"
    exit 1
fi

# Check if all complete
if [ $jet_step1_complete -eq $total_samples ] && \
   [ $jet_step2_complete -eq $total_samples ] && \
   [ $teprof2_complete -eq $total_samples ]; then
    echo ""
    echo "SUCCESS: All samples processed successfully!"
    exit 0
else
    echo ""
    echo "INFO: Processing incomplete"
    echo "  JET Step 1: ${jet_step1_complete}/${total_samples} complete"
    echo "  JET Step 2: ${jet_step2_complete}/${total_samples} complete"
    echo "  TEProf2: ${teprof2_complete}/${total_samples} complete"
    exit 2
fi
