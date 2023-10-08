#!/bin/bash

set -e

datapath=/home/zhepan/Project/PanCancerAtlas/data/NSCLC_28
path=$datapath/sra
out=$datapath/fastq


for i in $(ls $path);
do
  fasterq-dump --split-files -p -e 12 --include-technical --outdir $out/${i:0:11} $path/$i 
  cd $out/${i:0:11}
  pigz *
  rm -f $path/$i
  
  #cd $datapath/count
  
  #mv -u $datapath/fastq/${i:0:11}/${i:0:11}_1.fastq.gz $datapath/fastq/${i:0:11}/${i:0:11}'_S1_L001_R1_001.fastq.gz'
  #mv -u $datapath/fastq/${i:0:11}/${i:0:11}_2.fastq.gz $datapath/fastq/${i:0:11}/${i:0:11}'_S1_L001_R2_001.fastq.gz'
  
  #fastq=$datapath/fastq/${i:0:11}
  #cellranger count --id ${i:0:11} \
  #                 --localcores=12 \
  #                 --localmem=200 \
  #                 --transcriptome=$ref \
  #                 --no-bam --nosecondary --fastqs=$fastq  #--include-introns=false \
  #rm -rf $fastq
done
