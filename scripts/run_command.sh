
#!/bin/sh

./aggregateProcessedAnnotation.R -a ./arguments.txt
mkdir filterreadstats
python ./commandsmax_speed.py filter_combined_candidates.tsv ./
bash filterreadcommands.sh

find ./filterreadstats -name "*.stats" -type f -maxdepth 1 -print0 | xargs -0 -n128 -P1 grep e > resultgrep_filterreadstatsdone.txt
cat resultgrep_filterreadstatsdone.txt | sed 's/\:/\t/g' > filter_read_stats.txt
./filterReadCandidates.R
gffread -E candidate_transcripts.gff3 -T -o candidate_transcripts.gtf
echo candidate_transcripts.gtf > cuffmergegtf.list
cuffmerge -o ./merged_asm_full -g ../genome_46/gencode.v46.basic.annotation.sorted.gtf cuffmergegtf.list
mv ./merged_asm_full/merged.gtf reference_merged_candidates.gtf
gffread -E reference_merged_candidates.gtf -o- > reference_merged_candidates.gff3
python ./rmskhg38_annotate_gtf_update_test_tpm_cuff.py reference_merged_candidates.gff3 ./arguments.txt
find ./ -name "*bam" | while read file ; do xbase=${file##*/}; echo "samtools view -q 255 -h "$file" | stringtie - -o "${xbase%.*}".gtf -e -b "${xbase%.*}"_stats -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf" >> quantificationCommands.txt ; done ;
bash quantificationCommands.txt

./mergeAnnotationProcess.R

find . -name "*i_data.ctab" > ctab_i.txt

cat ctab_i.txt | while read ID ; do fileid=$(echo "$ID" | awk -F "/" '{print $2}'); cat <(printf 'chr\tstrand\tstart\tend\t'${fileid/_stats/}'\n') <(grep -F -f candidate_introns.txt $ID | awk -F'\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6 }') > ${ID}_cand ; done ;

cat <(find . -name "*i_data.ctab_cand" | head -1 | while read file ; do cat $file | awk '{print $1"\t"$2"\t"$3"\t"$4}' ; done;) > table_i_all

find . -name "*i_data.ctab_cand" | while read file ; do paste -d'\t' <(cat table_i_all) <(cat $file | awk '{print $5}') > table_i_all_temp; mv table_i_all_temp table_i_all; done ;

ls ./*stats/t_data.ctab > ctablist.txt

cat ctablist.txt | while read file ; do echo "python ./stringtieExpressionFrac.py $file" >> stringtieExpressionFracCommands.txt ; done;

ls ./*stats/t_data.ctab_frac_tot > ctab_frac_tot_files.txt
ls ./*stats/t_data.ctab_tpm > ctab_tpm_files.txt

cat <(echo "TranscriptID") <(find . -name "*ctab_frac_tot" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_frac_tot
cat ctab_frac_tot_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_frac_tot) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_frac_tot_temp; mv table_frac_tot_temp table_frac_tot; done ;

cat <(echo "TranscriptID") <(find . -name "*ctab_tpm" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_tpm
cat ctab_tpm_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_tpm) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_tpm_temp; mv table_tpm_temp table_tpm; done ;

cat <(head -1 table_frac_tot) <(grep -Ff candidate_names.txt table_frac_tot) > table_frac_tot_cand
cat <(head -1 table_tpm) <(grep -Ff candidate_names.txt table_tpm) > table_tpm_cand
./finalStatisticsOutput.R -a ./arguments.txt
