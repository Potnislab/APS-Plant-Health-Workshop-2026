#!/bin/bash

module load anaconda/3-2020.02
prokka \
  --prefix 91-118 \
  --locustag 91-118 \
  --increment 10 \
  --outdir 91-118_Prokka \
  --force \
  --addgenes \
  --genus Xanthomonas \
  --gcode 11 \
  GCF_000192045.2_ASM19204v3_genomic.fna

module load parsnp/1.5.6
parsnp -r ./91-118_prokka/91-118.fna -g ./91-118_prokka/91-118.gbk  -d ../fna/*.fna -o Xp_parsnp_output -p 30

wget https://github.com/marbl/harvest-tools/releases/download/v1.2/harvesttools-Linux64-v1.2.tar.gz
tar -xvf harvesttools-Linux64-v1.2.tar.gz
cd harvesttools-Linux64-v1.2/
./harvesttools \
  -x ../Xp_parsnp_output/parsnp.xmfa \
  -i ../Xp_parsnp_output/parsnp.ggr \
  -M ../SNPs.fa
cd ..
