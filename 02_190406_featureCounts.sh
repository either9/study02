#!/bin/sh

. ./config.txt

mkdir ${results_folder}${results_featureCounts_STAR_100}
mkdir ${results_folder}${results_featureCounts_STAR_50}

## main
for seqlib in ${seq_list[@]}
do

featureCounts \
-T 24 \
-t exon \
-g gene_id \
-a ${ref_gtf} \
-o ${results_folder}${results_featureCounts_STAR_100}/counts_${seqlib}.txt \
${results_folder}${results_STAR_100}/${seqlib}Aligned.out.sam


featureCounts \
-T 24 \
-t exon \
-g gene_id \
-a ${ref_gtf} \
-o ${results_folder}${results_featureCounts_STAR_50}/counts_${seqlib}.txt \
${results_folder}${results_STAR_50}/${seqlib}Aligned.out.sam

done
