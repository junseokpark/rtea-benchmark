# 🎯 QUICK ACTION ITEMS - What You Need to Do

## ✅ Step 1: Provide TEProf2 Command (REQUIRED)

I need the **actual TEProf2 command structure** just like you provided for JET.

Please share:
- The exact command you use to run TEProf2
- All parameters and flags
- Expected input/output file structure

Example format (like you did for JET):
```bash
# Your actual TEProf2 command
singularity exec ${TEProf2} /path/to/script.sh \
    --param1 value1 \
    --param2 value2 \
    ...
```

## ✅ Step 2: Fill in Configuration (REQUIRED)

Edit `config_template.sh` → Save as `config.sh`

**Must Update These Paths:**
```bash
# Line 15: JET installation
JETProjectDir="/path/to/JET"  # ← YOUR PATH HERE

# Line 19: Samtools binary
samtoolsBinDir="/path/to/samtools/bin"  # ← YOUR PATH HERE

# Line 20: STAR binary
starBinDir="/path/to/STAR/bin"  # ← YOUR PATH HERE

# Line 37: R library
RlibDir="/path/to/R/library"  # ← YOUR PATH HERE

# Lines 24-26: Verify these are correct
readLength=150  # Your read length?
organism="human"  # Correct?
genome="hg38"  # Correct?
database="database_name"  # What's the actual name?
```

## ✅ Step 3: Verify Reference Files Exist

Check these files are present:
```bash
ls -lh $DATA_HOME/ref/reference.fa
ls -lh $DATA_HOME/ref/gene_annotation.gtf
ls -lh $DATA_HOME/ref/TE_annotation.gff
ls -lh $DATA_HOME/ref/repeats.txt
ls -d $DATA_HOME/ref/STAR_indexes/
```

If STAR indexes don't exist, create them:
```bash
STAR --runMode genomeGenerate \
     --genomeDir $DATA_HOME/ref/STAR_indexes \
     --genomeFastaFiles $DATA_HOME/ref/reference.fa \
     --sjdbGTFfile $DATA_HOME/ref/gene_annotation.gtf \
     --runThreadN 16
```

## 📋 Reference File Format Questions

Please confirm:

1. **repeats.txt format** - What does it look like?
   - One repeat per line?
   - Tab-delimited with coordinates?
   - Example of first few lines?

2. **TE_annotation.gff** - Is this standard GFF3 format?
   - Or do you need to convert from GTF/BED?

## 🔧 Files to Use

**For the updated JET workflow:**
- ✅ `config_template.sh` → Configure and save as `config.sh`
- ✅ `process_array_updated.sh` → Main processing script
- ✅ `check_status_updated.sh` → Status monitoring
- ✅ `generate_sample_list.sh` → Sample discovery (unchanged)

**Original files (for reference only):**
- ⚠️ `process_array.sh` - Old version (BWA-based)
- ⚠️ `check_status.sh` - Old version

## 📝 Integration Checklist

Before running:
- [ ] Filled in config.sh with all paths
- [ ] Validated with: `source config.sh && validate_config`
- [ ] STAR genome indexes exist
- [ ] Provided actual TEProf2 command
- [ ] Created logs directory: `mkdir -p logs`
- [ ] Tested JET Step 1 manually on one sample
- [ ] Tested JET Step 2 manually on one sample

## 🚀 Once Everything is Ready

```bash
# 1. Generate sample list
./generate_sample_list.sh

# 2. Edit process_array_updated.sh
# Add near the top (line ~20):
source /path/to/your/config.sh

# 3. Update array size based on sample count
# Edit line: #SBATCH --array=1-N%20
# Where N = number from: wc -l sample_list.txt

# 4. Submit
chmod +x process_array_updated.sh
sbatch process_array_updated.sh

# 5. Monitor
./check_status_updated.sh
```

## ❓ Questions for You

1. **What's the actual TEProf2 command?**
2. **What format is repeats.txt?**
3. **Are all reference files in the expected format?**
4. **Do you need help creating STAR indexes?**
5. **Any special R packages needed for JET Step 2?**

## 📞 Once You Provide This Info

I'll create the **final, production-ready** package with:
- ✅ Correct TEProf2 integration
- ✅ Pre-configured scripts
- ✅ Complete documentation
- ✅ Testing checklist
- ✅ Troubleshooting guide

---

**Most Important:** Send me the actual TEProf2 command structure!
