#intersect with coords with full exon coords (no overlaps dups etc removed)
# intersect coords without flanks- to not mess up alignment, then add flansk back eg
tail -n +2 /home/maria/cactus/to_human_ref/coord_filt_no_ins.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 + 1, $3 - 1, $4, $5, $6}' | \
awk -F'\t' 'BEGIN{OFS="\t"} ($3>$2) {print}' | \
bedtools intersect -wa -wb -a - -b /home/maria/filter_transcripts/output/exon_merged_ids_sort.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $5, $6, $7, $8 - 1, $9 + 1, $10, $12}' \
> "/home/maria/cactus_target_size/auxillary/extracted_df2.bed"

#redoing with primates
tail -n +2 /home/maria/cactus/to_human_ref/primates.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 + 1, $3 - 1, $4, $5, $6, $7, $8}' | \
awk -F'\t' 'BEGIN{OFS="\t"} ($3>$2) {print}' | \
bedtools intersect -wa -wb -a - -b /home/maria/filter_transcripts/output/exon_merged_ids_sort.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $5, $6, $7, $8, $9, $10-1, $11+1, $12, $14}' \
> "/home/maria/cactus/to_human_ref/primates_intersected.bed"

#to redo withoutcpg
tail -n +2 /home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 + 1, $3 - 1, $4, $5, $6}' | \
awk -F'\t' 'BEGIN{OFS="\t"} ($3>$2) {print}' | \
bedtools intersect -wa -wb -a - -b /home/maria/filter_transcripts/output/exon_merged_ids_sort.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $5, $6, $7, $8 - 1, $9 + 1, $10, $12}' \
> "/home/maria/cactus_target_size/auxillary/extracted_df2_nocpg.bed"


#old

#add flanks to coords, and instersect with befile 
# rm top line of file
#no longer use, now use overlapping 
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $6}' /home/maria/filter_transcripts/output/exon_cords_ids.bed \
> /home/maria/cactus_target_size/auxillary/coord_flank.bed


tail -n +2 /home/maria/cactus/to_human_ref/full_sing_filtgaps.maf_no_ins.bed | \
bedtools subtract -a - -b /home/maria/cactus_target_size/auxillary/overlaps.bed | \
bedtools intersect -wa -wb -a - -b /home/maria/cactus_target_size/auxillary/coords_no_overlaps.bed \
> /home/maria/cactus_target_size/auxillary/extracted_df.bed




#because i extracted the hclca coords using the inputted exon coords in theory should be fine, 
#always should have hclca inside coord_flank
#if not then may need to identify acual overlap and cut

#intersect without flanks, to avoid artifical overlaps?
#plus remove overlapping segments from gene sequences 
#next time will prepocess to remove before extractinge seqs
tail -n +2 /home/maria/cactus/to_human_ref/full_sing_filtgaps.maf_no_ins.bed | \
bedtools subtract -a - -b /home/maria/cactus_target_size/auxillary/overlaps.bed | \
bedtools intersect -wa -wb -a - -b /home/maria/cactus_target_size/auxillary/coords_no_overlaps.bed \
> /home/maria/cactus_target_size/auxillary/extracted_df.bed

#rm /home/maria/cactus_target_size/auxillary/coord_flank.bed