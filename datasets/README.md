# Project: Pan-Cancer Atlas  

`python` `R` `single-cell` `multi-omics`

## Info  

Git repo: [kiklata/PanCancerAtlas](https://github.com/kiklata/PanCancerAtlas)  
Data: [xiyoucloud](https://xiyoucloud.pro:12001)  
Cloud storage: [aliyun](https://www.aliyundrive.com/drive/file/backup/64c8c567f53b025d56db4531ab96e3292f7f0903), [baiduyun](https://pan.baidu.com/disk/main#/index?category=all&path=%2Fapps%2Fbypy%2FProject%2FPanCancerAtlas)

## Analysis Process  

 ```ditaa {cmd=true args=["-E"]}
                                                                   part
  part                              all           +=------------------------------------+ 
/------\  10X: cellranger  +-------------------+  |  +------+  +---------+  +---------+ |
| fastq|------------------>| raw/filter matrix |  |  | ATAC |  | TCR/BCR |  | spatial | |
| bam  |  other: STAR  |   +-------------------+  |  +------+  +---------+  +---------+ |
\------/               |                          +-------------------------------------+
                       |                           part
                       |  +=----------------------------------------------------+
                       |  |  /---------\  /--------\  /----------\  /--------\  |
                       |  |  |  cnv    |  |  velo  |  | splicing |  | Mixcr  |  |
                       +->|  |  numbat |  |  STAR  |  | STAR     |  | TRUST4 |  |
                          |  \---------/  \--------/  \----------/  \--------/  |
                     |    +-----------------------------------------------------+        |          
                     +-------------------------------------------------------------------+
                                         |
                                         v
 ------------------------------------------------------------------------------------------------------     
           +=----------------------------------+   
           |    +------+         +------+      |       /---------\
           |    | seu  |         | h5ad |      |       | metadata|
           |    +------+         +------+      |       \---------/
           |      | Seurat        | Scanpy     |          | study
           |      | Signac        | mudata     |          | date
           |      | immuneRep     | scirpy     |          | platform
           |      +----------     | scverse    |          | type: cancer/normal
           |                      +---------   |          | treatment
           +-----------------------------------+          | timepoint
                                                          | patient
                                                          | sample
                                                          | celltype (by original paper)
                                                          +------------------------------
   ----------------------------------------------------------------------------------------------
                                         |
                                         v
                                   quality control
        +-------------------------------+         /-- nCount  
        | filter out low quality cells  |---------|   nFeature
        +-------------------------------+         \-- mt_percent
             manually OR check MAD = median(abs(diff(Xi,Xm))), 3 or 5 MAD?
                                           
        +--------------+         /-- Seu ---> SoupX (need raw & filter matrix) DecontX           
        | ambient RNA  |---------|         
        +--------------+         \-- h5ad -->   CellBender (GPU)
                      need Mod to get INT count
        +--------------+    
        |   doublet    |--------> scDbiFinder (benchmark)
        +--------------+    
     -----------------------------------------------------------------------------------------
                                         |              
                                         v
                               individual annotation
                     manually + automatic methods (see benchmark)
                 pick subset to integrate, use scib try multiple methods
                     use metacell-like methods to downsample
                                         |
                                         |
                                         v
+---------+   +----------------+   +---------------+   +---+   +-----------+   +----------+   /-- metastasis
|EPI cells|---|inferCNV/copykat|---|malignant cells|---|NMF|---|MetaProgram|---|phenotypes|---|   response
+---------+   +----------------+   +---------------+   +---+   +-----------+   +----------+   \-- adverse events ...
                                                                                   select from metadata
      Immune cells
      Stromal cells
      ...etc
   -----------------------------------------------------------------------------------------------------
                                         |
                                         v
                                     multiomics
 TBD
```

## Reference

### protocols

[sc-best-practice](https://www.sc-best-practices.org/preamble.html)
[Annotation](https://www.nature.com/articles/s41596-021-00534-0)

### software

[STAR](https://github.com/alexdobin/STAR)
[numbat](https://github.com/kharchenkolab/numbat)
[velocyto](http://velocyto.org/)
[Mixcr](https://mixcr.com/)
[TRUST4](https://github.com/liulab-dfci/TRUST4)
[scRepertoire](https://github.com/ncborcherding/scRepertoire)
[scverse](https://github.com/scverse)
[SoupX](https://github.com/constantAmateur/SoupX)
[DecontX](http://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html)
[CellBender](https://github.com/broadinstitute/CellBender)
[scDbiFinder](https://github.com/plger/scDblFinder)
[scib](https://github.com/theislab/scib)

### some papers

[PaperPDF](/ref/)