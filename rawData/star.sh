#!/bin/bash
# STAR env
set -e # Exit immediately if a command exits with a non-zero status

datapath=/home/zhepan/Project/PanCancerAtlas/data/ESCC_1

refer=/home/zhepan/Reference/STAR_index
whitelist=/home/zhepan/Reference/STAR_index/3M-february-2018.txt

for sample in $(echo $(diff <(cd $datapath/star && ls | sort) <(cd $datapath/count && ls | sort)) | grep 'SRR........' -o | uniq);
do
  bam=$datapath/count/$sample/outs/*.bam
  mkdir -p $datapath/star/$sample
  STAR --runThreadN 12 \
       --genomeDir $refer \
       --soloType CB_UMI_Simple \
       --readFilesType SAM SE --readFilesIn $bam  \
       --readFilesCommand samtools view -F 0x100 \
       --soloInputSAMattrBarcodeSeq CR UR \
       --soloInputSAMattrBarcodeQual CY UY \
       --soloCBwhitelist $whitelist \
       --soloCellFilter EmptyDrops_CR \
       --soloFeatures Gene SJ Velocyto \
       --outFileNamePrefix $datapath/star/$sample/ \
       --outSAMtype None
  rm -rf $bam
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=star
