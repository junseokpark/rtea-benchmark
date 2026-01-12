RepeatMasker="/home/junseokp/workspaces/data/rTea-simul/ref/hg38/hg38.RepeatMasker-4.0.6-Dfam-2.0.fa.out"


# awk 'BEGIN{OFS="\t"}
# $1 ~ /^[0-9]+$/ {
#   rep=$10
#   split($11,a,"/")
#   cls=a[1]
#   fam=(a[2] ? a[2] : a[1])
#   print rep,cls,fam
# }' "$RepeatMasker" \
# | sort -u \
# > repeatmasker_hg38_description_uniq.lst


awk 'BEGIN{
        OFS="\t";
        print "#genoName","genoStart","genoEnd","strand","repName","repClass","repFamily"
     }
     !/^#/ && NF>=15 {
        chr=$5;
        start=$6;
        end=$7;
        strand=($9=="C" ? "-" : "+");
        repName=$10;
        split($11, a, "/");
        repClass=a[1];
        repFamily=(a[2] ? a[2] : a[1]);
        print repName,repFamily,repClass
     }' $RepeatMasker | sort | uniq >  repeatmasker_hg38_description_uniq.lst