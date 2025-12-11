# UPDATE NOTES - JET Pipeline Integration

## What Changed

I've updated the pipeline scripts to match the **actual JET workflow** you provided, which uses:
- **STAR aligner** (not BWA-MEM)
- **Metadata file format** for input
- **Two-step process**: Step1 (STAR alignment) + Step2 (R analysis)

## 📋 Required Configuration

Before running the pipeline, you **MUST** configure these paths:

### 1. Create Your Configuration File

```bash
cp config_template.sh config.sh
nano config.sh  # or vim, emacs, etc.
```

### 2. Fill in These REQUIRED Paths:

```bash
# JET Paths
JETProjectDir="/path/to/JET"              # Where JET is installed
samtoolsBinDir="/path/to/samtools/bin"    # Samtools binary directory
starBinDir="/path/to/STAR/bin"            # STAR binary directory
RlibDir="/path/to/R/library"              # R library for JET Step 2

# Reference Files
starIndexesDir="${REF_DIR}/STAR_indexes"  # STAR genome indexes
repeatsFile="${REF_DIR}/repeats.txt"      # Repeat elements file
gffFile="${REF_DIR}/TE_annotation.gff"    # TE annotation in GFF format

# Optional: Update if different
readLength=150        # Your sequencing read length
organism="human"      # Your organism
genome="hg38"         # Your genome version
database="db_name"    # Database name for JET
```

### 3. Validate Configuration

```bash
source config.sh
validate_config
```

## 📁 Updated Files

### New Files (Updated for actual JET):
1. **process_array_updated.sh** - SLURM array job with correct JET commands
2. **check_status_updated.sh** - Status checker for JET Step 1 & 2 outputs
3. **config_template.sh** - Configuration template with all required paths

### Original Files (Still Valid):
- generate_sample_list.sh - Sample discovery (no changes needed)
- pipeline.sh - Master control script (minor updates needed)

## 🔧 How JET Works Now

### Metadata File Format
For each sample, the script creates a metadata file:
```
sample          fastq1                                  fastq2
sim200_AluY_blood  /path/to/sim200_AluY_blood.1.fq.gz  /path/to/sim200_AluY_blood.2.fq.gz
```

### Directory Structure Per Sample
```
output/{rel_path}/{sample_name}/
├── metadata.txt                    # Created automatically
├── output/                         # JET results
│   ├── Aligned.sortedByCoord.out.bam  # Step 1 output
│   └── *_TE_insertions.bed            # Step 2 output
├── log/                            # JET logs
│   ├── step1_multisample_running_YYYYMMDD.log
│   └── step2_multisample_running_YYYYMMDD.log
├── err/                            # Error logs
└── TEProf2/                        # TEProf2 results
```

## ⚠️ Important Notes

### 1. STAR Index Must Exist
Before running, you need STAR genome indexes:
```bash
STAR --runMode genomeGenerate \
     --genomeDir ${REF_DIR}/STAR_indexes \
     --genomeFastaFiles ${REF_DIR}/reference.fa \
     --sjdbGTFfile ${REF_DIR}/gene_annotation.gtf \
     --runThreadN 16
```

### 2. TEProf2 Command Structure
The TEProf2 command is still a **placeholder**. Please provide the actual TEProf2 command format, similar to what you provided for JET.

Current placeholder in script:
```bash
singularity exec ${TEProf2} teprof2 \
    --fq1 ${FQ1} \
    --fq2 ${FQ2} \
    --ref ${refDir}/reference.fa \
    --te-annot ${refDir}/TE_annotation.gtf \
    ...
```

### 3. File Format Requirements
- **repeats.txt**: Format required by JET Step 2
- **TE_annotation.gff**: Must be GFF format (not BED or GTF)
- Gene annotation can be GTF format

## 🚀 Quick Start (Updated)

```bash
# 1. Configure
cp config_template.sh config.sh
nano config.sh  # Fill in all REQUIRED paths
source config.sh
validate_config

# 2. Generate sample list
./generate_sample_list.sh

# 3. Update process_array_updated.sh to source config.sh
# Add this at the top of the script:
# source /path/to/config.sh

# 4. Submit jobs
chmod +x process_array_updated.sh
sbatch process_array_updated.sh

# 5. Monitor
watch -n 60 'squeue -u $USER'
./check_status_updated.sh
```

## 🐛 Debugging Tips

### Check JET Step 1 Logs
```bash
tail -f output/{path}/{sample}/log/step1_multisample_running_*.log
```

### Check JET Step 2 Logs
```bash
tail -f output/{path}/{sample}/log/step2_multisample_running_*.log
```

### Check SLURM Logs
```bash
tail -f logs/TE_array_*.err
```

### Verify STAR Alignment Output
```bash
samtools view output/{path}/{sample}/output/Aligned.sortedByCoord.out.bam | head
```

## 📊 Expected Output Sizes

Per sample (rough estimates):
- **Step 1 (STAR)**: 2-10 GB BAM file
- **Step 2 (R)**: 10-100 MB result files
- **TEProf2**: 50-500 MB
- **Logs**: 1-10 MB

Total for 200 samples: **500 GB - 2 TB**

## ❓ What I Still Need From You

1. **Actual TEProf2 command structure**
   - Similar to how you provided JET Step 1 & 2
   - What are the exact parameters and flags?
   - What's the expected output structure?

2. **Confirm these file formats:**
   - What format is your repeats.txt file?
   - Is TE_annotation.gff in standard GFF3 format?

3. **Any JET-specific requirements:**
   - Special environment variables needed?
   - Required R packages for Step 2?
   - Any pre-processing steps?

## 🔄 Migration from Old Scripts

If you already ran the old scripts:
1. The **sample_list.txt** is still valid - no need to regenerate
2. The **output structure is different** - old outputs won't be found by new check_status
3. You may want to **clean and restart** or keep old results separate

## ✅ Checklist Before Running

- [ ] Copied config_template.sh to config.sh
- [ ] Filled in all REQUIRED paths in config.sh
- [ ] Ran validate_config successfully
- [ ] STAR genome indexes exist
- [ ] All reference files exist and are readable
- [ ] Generated sample_list.txt
- [ ] Updated process_array_updated.sh to source config.sh
- [ ] Created logs directory: `mkdir -p logs`
- [ ] Provided actual TEProf2 command (if different from placeholder)

---

## 📞 Next Steps

Once you provide:
1. **Actual TEProf2 command structure**
2. **Confirmation of file formats**

I can finalize the scripts and create the final package for download.
