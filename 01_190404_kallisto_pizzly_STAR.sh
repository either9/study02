#!/bin/bash

. ./config.txt

## activating conda environment
#conda activate gneral

## kallisto-index building
kallisto index -i ${ref_kallisto_index} ${ref_cdna_fa}

## STAR index building sjdbOverhang 100
mkdir ${ref_STAR_index_100}
STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ${ref_STAR_index_100} \
--genomeFastaFiles ${ref_genome_fa} \
--sjdbGTFfile ${ref_gtf} \
--sjdbOverhang 100

## STAR index building sjdbOverhang 50
mkdir ${ref_STAR_index_50}
STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ${ref_STAR_index_50} \
--genomeFastaFiles ${ref_genome_fa} \
--sjdbGTFfile ${ref_gtf} \
--sjdbOverhang 50


## mkdir
mkdir ${results_folder}${results_STAR_100}
mkdir ${results_folder}${results_STAR_50}
mkdir ${results_folder}${results_kallisto}
mkdir ${results_folder}${results_pizzly}

## main
for seqlib in ${seq_list[@]}
do

## mkdir for kallsito
mkdir ${results_folder}${results_kallisto}/${seqlib}

## kallisto quant
kallisto quant \
-i ${ref_kallisto_index} \
-o ${results_folder}${results_kallisto}/${seqlib} \
-b 100 \
--fusion \
-t 12 \
${seq}/${seqlib}_1.fastq.gz \
${seq}/${seqlib}_2.fastq.gz


## pizly
pizzly \
-k 31 \
--gtf ${ref_gtf} \
--align-score 2 \
--insert-size 400 \
--fasta ${ref_genome_fa} \
--output ${results_pizzly}/${seqlib} \
${results_folder}${results_kallisto}/${seqlib}/fusion.txt


## STAR with 100
STAR --runThreadN 12 \
--genomeDir ${ref_STAR_index_100} \
--sjdbGTFfile ${ref_gtf} \
--readFilesCommand gunzip -c \
--readFilesIn \
${seq}/${seqlib}_1.fastq.gz \
${seq}/${seqlib}_2.fastq.gz \
--outFileNamePrefix ${results_folder}${results_STAR_100}/${seqlib}

## STAR with 50
STAR --runThreadN 12 \
--genomeDir ${ref_STAR_index_50} \
--sjdbGTFfile ${ref_gtf} \
--readFilesCommand gunzip -c \
--readFilesIn \
${seq}/${seqlib}_1.fastq.gz \
${seq}/${seqlib}_2.fastq.gz \
--outFileNamePrefix ${results_folder}${results_STAR_50}/${seqlib}

done
