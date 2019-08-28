#!/bin/sh

. ./config.txt

## activating conda environment
#conda activate general

## STAR index building sjdbOverhang 100
#mkdir ${ref_STAR_index_100}
#STAR \
#--runThreadN 12 \
#--runMode genomeGenerate \
#--genomeDir ${ref_STAR_index_100} \
#--genomeFastaFiles ${ref_genome_fa} \
#--sjdbGTFfile ${ref_gtf} \
#--sjdbOverhang 100

mkdir ${results_folder}${results_miRNA}

for seqlib in ${miRNA_seq_list[@]}
do

STAR \
--runThreadN 12 \
--genomeDir ${ref_STAR_index_100} \
--readFilesIn ${miRNA_seq}/${seqlib}_trimmed.fq.gz \
--readFilesCommand gunzip -c \
--alignIntronMin 1 \
--outFileNamePrefix ${results_folder}${results_miRNA}/${seqlib}_

featureCounts \
-T 24 \
-t "miRNA" \
-g Name \
-a ${ref_miRNA_gtf} \
-o ${results_folder}${results_miRNA}/counts_${seqlib}.txt \
-s 1 \
-R SAM \
-M \
${results_folder}${results_miRNA}/${seqlib}_Aligned.out.sam

done
