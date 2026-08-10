#!/bin/bash
module load strainest
module load bowtie2/2.5.1
source /apps/profiles/modules_asax.sh.dyn
module load bowtie2/2.5.1
module load samtools/1.11

strainest mapgenomes Xp91-118.fna Xeu_85-10.fna Xg_USA01-FM3.fna Xv_LMG911.fna Xp91-118.fna MR.fna
strainest map2snp Xp91-118.fna MR.fna snp.dgrp
strainest snpdist snp.dgrp snp_dist.txt hist.pdf
strainest snpclust snp.dgrp snp_dist.txt snp_clust.dgrp clusters.txt
strainest mapgenomes Xp91-118.fna Xeu_85-10.fna Xg_USA01-FM3.fna Xv_LMG911.fna Xp91-118.fna MA.fna

bowtie2-build -f MA.fna MA

bowtie2 --very-fast --no-unal -x MA -1 sample_1.fq -2 Sample_2.fq  -S Sample.sam 
samtools view -b Sample.sam > Sample.bam
samtools sort Sample.bam -o Sample.sorted.bam
samtools index Sample.sorted.bam 

strainest est ./snp_clust.dgrp sample.sorted.bam ./Sample_Inter