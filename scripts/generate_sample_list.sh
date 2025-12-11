#!/bin/bash

# Generate sample list for array job processing
DATA_HOME=/home/junseokp/workspaces/data/rTea-simul
OUTPUT_FILE="sample_list.txt"

echo "Generating sample list..."
echo "# Sample_Name FQ1 FQ2 Output_Path" > ${OUTPUT_FILE}

# Function to add sample to list
add_sample() {
    local fq1=$1
    local fq2=$2
    local rel_path=$3
    local sample_name=$(basename "$fq1" | sed 's/.1.fq.gz$//')
    
    echo "${sample_name} ${fq1} ${fq2} ${rel_path}" >> ${OUTPUT_FILE}
}

# Process nonReferenceTE samples
for te_type in AluY L1HS LTR5 SVA_F; do
    for coverage in 5X 10X 50X 100X 200X; do
        fq_dir="${DATA_HOME}/nonReferenceTE/${te_type}/${coverage}/fq"
        
        if [ -d "$fq_dir" ]; then
            for fq1 in ${fq_dir}/*.1.fq.gz; do
                if [ -f "$fq1" ]; then
                    fq2="${fq1%.1.fq.gz}.2.fq.gz"
                    if [ -f "$fq2" ]; then
                        rel_path="nonReferenceTE/${te_type}/${coverage}"
                        add_sample "$fq1" "$fq2" "$rel_path"
                    fi
                fi
            done
        fi
    done
done

# Process referenceTE/intron samples
for coverage in 5X 10X 50X 100X 200X; do
    # Regular samples
    fq_dir="${DATA_HOME}/referenceTE/intron/${coverage}/fq"
    if [ -d "$fq_dir" ]; then
        for fq1 in ${fq_dir}/reftefu403_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                if [ -f "$fq2" ]; then
                    rel_path="referenceTE/intron/${coverage}/fq"
                    add_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
    
    # Mutated samples
    fq_mut_dir="${DATA_HOME}/referenceTE/intron/${coverage}/fq_mut"
    if [ -d "$fq_mut_dir" ]; then
        for fq1 in ${fq_mut_dir}/reftefu403_mut_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                if [ -f "$fq2" ]; then
                    rel_path="referenceTE/intron/${coverage}/fq_mut"
                    add_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
done

# Process referenceTE/TSS samples
for coverage in 5X 10X 50X 100X 200X; do
    # Regular samples
    fq_dir="${DATA_HOME}/referenceTE/TSS/${coverage}/fq"
    if [ -d "$fq_dir" ]; then
        for fq1 in ${fq_dir}/reftetss270_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                if [ -f "$fq2" ]; then
                    rel_path="referenceTE/TSS/${coverage}/fq"
                    add_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
    
    # Mutated samples
    fq_mut_dir="${DATA_HOME}/referenceTE/TSS/${coverage}/fq_mut"
    if [ -d "$fq_mut_dir" ]; then
        for fq1 in ${fq_mut_dir}/reftetss270_mut_*_${coverage}.1.fq.gz; do
            if [ -f "$fq1" ]; then
                fq2="${fq1%.1.fq.gz}.2.fq.gz"
                if [ -f "$fq2" ]; then
                    rel_path="referenceTE/TSS/${coverage}/fq_mut"
                    add_sample "$fq1" "$fq2" "$rel_path"
                fi
            fi
        done
    fi
done

# Count samples
total_samples=$(grep -v "^#" ${OUTPUT_FILE} | wc -l)
echo "Generated sample list with ${total_samples} samples"
echo "Output: ${OUTPUT_FILE}"

# Show summary
echo ""
echo "Summary by category:"
echo "===================="
echo "nonReferenceTE samples:"
grep "nonReferenceTE" ${OUTPUT_FILE} | wc -l
echo ""
echo "referenceTE/intron samples:"
grep "referenceTE/intron" ${OUTPUT_FILE} | wc -l
echo ""
echo "referenceTE/TSS samples:"
grep "referenceTE/TSS" ${OUTPUT_FILE} | wc -l
