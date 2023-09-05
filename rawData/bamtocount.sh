set -e # Exit immediately if a command exits with a non-zero status

datapath=/home/zhepan/Project/PanCancerAtlas/data/LIHC_9
fq_path=$datapath/fastq

ref=/home/zhepan/Reference/Cellranger

cd $datapath/count

for sample in $(ls $datapath/bam);
do
  bam=$datapath/bam/$sample
  fq_name=$(echo $sample | sed 's/_possorted_genome_bam.bam.1//')
  cellranger bamtofastq --nthreads=12 --reads-per-fastq=5000000000 $bam $fq_path/$fq_name
  
  fq_file=$(ls $fq_path/$fq_name)
  
  cellranger count --id $fq_name \
                   --localcores=12 \
                   --localmem=150 \
                   --transcriptome=$ref \
                   --no-bam --nosecondary --fastqs=$fq_path/$fq_name/$fq_file  #--include-introns=false \
  rm -rf $fq_path/$fq_name
  rm -rf $bam
done


curl -X POST https://sctapi.ftqq.com/SCT167047TyIQz75wID4jMUfJ4EGGbELEN.send?title=Done