#!/bin/sh

. ./config.txt

#conda activate general

for seqlib in ${miRNA_seq_list[@]}
do

#gzip -d -k ${miRNA_seq}/${seqlib}.fastq.gz
trim_galore -o ${miRNA_seq} ${miRNA_seq}/${seqlib}.fastq.gz

done
