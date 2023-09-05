set -e # Exit immediately if a command exits with a non-zero status

datapath=/home/zhepan/Project/PanCancerAtlas/data/BRCA_13

ref=/home/zhepan/Reference/Cellranger

for sample in $(ls $datapath/HRA002051);
do
  cd $datapath/count
  #if [ $(ls $datapath/fastq/$sample | wc -w) == "3" ]; then
  mv -u $datapath/HRA002051/$sample/${sample}_f1.fastq.gz $datapath/HRA002051/$sample/$sample'_S1_L001_R1_001.fastq.gz'
  mv -u $datapath/HRA002051/$sample/${sample}_r2.fastq.gz $datapath/HRA002051/$sample/$sample'_S1_L001_R2_001.fastq.gz'
  #mv -u $datapath/fastq/$sample/${sample}*_I1.fastq.gz $datapath/fastq/$sample/$sample'_S1_L001_I1_001.fastq.gz'
  fastq=$datapath/HRA002051/$sample
  #else 
  #mv -u $datapath/fastq/$sample/${sample}_1*_R1.fastq.gz $datapath/fastq/$sample/$sample'_S1_L001_R1_001.fastq.gz'
  #mv -u $datapath/fastq/$sample/${sample}_1*_R2.fastq.gz $datapath/fastq/$sample/$sample'_S1_L001_R2_001.fastq.gz'
  #mv -u $datapath/fastq/$sample/${sample}_1*_I1.fastq.gz $datapath/fastq/$sample/$sample'_S1_L001_I1_001.fastq.gz'
  #mv -u $datapath/fastq/$sample/${sample}_2*_R1.fastq.gz $datapath/fastq/$sample/$sample'_S1_L002_R1_001.fastq.gz'
  #mv -u $datapath/fastq/$sample/${sample}_2*_R2.fastq.gz $datapath/fastq/$sample/$sample'_S1_L002_R2_001.fastq.gz'
  #mv -u $datapath/fastq/$sample/${sample}_2*_I1.fastq.gz $datapath/fastq/$sample/$sample'_S1_L002_I1_001.fastq.gz'
  #fastq=$datapath/fastq/$sample
  #fi
  cellranger count --id $sample \
                   --localcores=12 \
                   --localmem=150 \
                   --transcriptome=$ref \
                   --chemistry=fiveprime \
                   --no-bam --nosecondary --fastqs=$fastq  #--include-introns=false \
  rm -rf $fastq
done

curl -X POST https://sctapi.ftqq.com/SCT167047TyIQz75wID4jMUfJ4EGGbELEN.send?title=Done