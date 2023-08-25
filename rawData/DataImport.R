library(Seurat)
library(dplyr)

# basic--------------
count = Read10X('~/Project/PanCancerAtlas/data/AML_7/raw/HTO_R306_R309')
count = Read10X_h5('~/Project/tempData/SKCM_10/GSE174401_filtered_feature_bc_matrix.h5')

meta = data.table::fread('~/Project/tempData/STAD_9/GSE234129_meta.tsv.gz') %>% as.data.frame()
rownames(meta) = meta$cell_barcodes
seu = CreateSeuratObject(count,min.cells = 3,min.features = 200)
seu = AddMetaData(seu,meta)
saveRDS(seu,file = '~/Project/tempData/STAD_9/seu/STAD_9.rds')

# named based samples------------
#subdir = 'F'
wd = paste0('~/Project/PanCancerAtlas/data/GCTB_2/')
all.sample = list.dirs(wd,recursive = F,full.names = F)
all.sample = list.files(wd,recursive = F,full.names = F,pattern = 'filter')

samples.name = stringr::str_sub(all.sample,12,-31)
samples.name = all.sample
samples.name = c('UM1','UM23')
#barcodes = list.files(wd,pattern = '.barcode',recursive = T)
#features = list.files(wd,pattern = '.gene',recursive = T)
#matrixs = list.files(wd,pattern = '.matrix',recursive = T)

seu.list=list()
for (i in 1:length(all.sample)) {

  count = Read10X(paste0(wd,all.sample[i]),)
  #count = count[,!colnames(count) %>% duplicated()]
  #count = count[rownames(count)[substring(rownames(count),1,6)=='GRCh38'],]
  #rownames(count) = gsub('GRCh38_','',x = rownames(count))
  #count = Read10X_h5(paste0(wd,all.sample[i]))
  
  seu = CreateSeuratObject(count,min.features = 200)
  #seu = CreateSeuratObject(count$`Gene Expression`,min.features = 200)
  #seu[['CMO']] = CreateAssayObject(counts = count$`Multiplexing Capture`[,colnames(seu)])
  
  #meta = read.delim(paste0(wd,i,'/metadata.csv.gz'))
  #rownames(meta) = meta[,1]
  #seu = AddMetaData(seu,meta)
  seu.list[[samples.name[i]]] = seu
}
saveRDS(seu.list,file = paste0(wd,'seu/seu.list.rds'))

# count basic----------------
count = data.table::fread('~/Project/PanCancerAtlas/data/HNSC_5/raw/GSM4546857_LSCC01_DBEC_UMI.csv.gz') %>% as.data.frame()
#count = Matrix::readMM('~/Project/tempData/SKCM_11/GSE200218_sc_sn_counts.mtx.gz')
count[1:10,1:10]
meta = read.delim("~/Project/tempData/NSCLC_27/GSE193531_cell-level-metadata.csv.gz",sep = ',' )
rownames(meta) = meta[,1]
meta = meta[,-1]
#meta =t(meta) %>% as.data.frame()
#meta_c = count[1:2,]
#meta[,1] = NULL
#meta = t(meta) %>% as.data.frame()
#count = count[-c(1:3),]
#meta <- read.delim("~/Project/PanCancerAtlas/data/LIHC_7/GSE98638.txt")
#count = count[!is.na(count[,3]),]
#count = count[nchar(count[,3]) != 0,]
#count = count[!duplicated(count[,1]),]

rownames(count) = count[,1]

#count[,colnames(count) %>% grep('Gene',.)] = NULL
count[,c(1)] = NULL
#count = t(count) %>% as.data.frame()
#rownames(meta) = meta[,1]
seu = CreateSeuratObject(count)
seu = AddMetaData(seu,metadata = meta)
saveRDS(seu,file ='~/Project/tempData/NSCLC_27/NSCLC_27.rds')


# 3CA ------------
count = Matrix::readMM('~/Project/PanCancerAtlas/data/NSCLC_12/3CA/Exp_data_counts.mtx') %>% as.data.frame()
cell = read.delim("~/Project/PanCancerAtlas/data/NSCLC_12/3CA/Cells.csv",sep = ',')
genes = read.delim("~/Project/PanCancerAtlas/data/NSCLC_12/3CA/Genes.txt",header = F) 

count[1:5,1:5]


colnames(count) = cell[,1]
count[,'gene'] = genes$V1
count = count[!duplicated(count[,'gene']),]

rownames(count) = count[,'gene']
count[,'gene'] = NULL

rownames(cell) = cell[,1]
seu = CreateSeuratObject(count,min.cells=0,min.features=200)
seu = AddMetaData(seu,cell)
saveRDS(seu,file ='~/Project/PanCancerAtlas/data/NSCLC_12/NSCLC_12.rds')
saveRDS(cell,file = '~/Project/PanCancerAtlas/data/NSCLC_12/cell_meta.rds')

# count samples--------------

dataset = c('~/Project/PanCancerAtlas/data/MM_15/')
filea = list.files(dataset,recursive = F,full.names = F,pattern = 'gz')
#tags = list.files(dataset,recursive = F,full.names = F,pattern = 'Sample')
samples = stringr::str_sub(filea,12,-26)
seu.list=list()
count = data.table::fread(paste0(dataset,filea[2])) %>% as.data.frame()
#tag = data.table::fread(paste0(dataset,tags[1])) %>% as.data.frame()
count[1:5,1:5]
#meta = data.table::fread('~/Project/tempData/SKCM_11/GSE200218_sc_sn_metadata.csv.gz') %>% as.data.frame()
#meta$V1 = NULL
#meta = meta[!duplicated(meta[,1]),]
#rownames(meta) = meta[,1]
#meta[,1] =NULL
#colnames(meta)[1] = 'CellID'
#symbol = bitr(count[,1],'ENSEMBL','SYMBOL','org.Hs.eg.db',drop = T)
#symbol = symbol[!duplicated(symbol$ENSEMBL),] 
#colnames(symbol)[1] = 'V1'
#count = left_join(count,symbol,'V1')

#meta = read.csv("~/Project/PanCancerAtlas/data/PRAD_6/meta/GSE143791_cell.annotation.human.csv.gz")
#rownames(meta) = meta[,1]
for (i in 1:length(filea)) {
#count = data.table::fread(paste0(dataset,filea[i])) %>%
# as.data.frame() %>% reshape2::dcast(.,Cell_Index~Gene,value.var = 'RSEC_Adjusted_Molecules') %>% data.table::setnafill(.,fill = 0)
count = data.table::fread(paste0(dataset,filea[i])) %>% as.data.frame() 
#tag = data.table::fread(paste0(dataset,tags[i])) %>% as.data.frame()
#colnames(count) = paste0('cell_',count[1,])
#count = count[-1,]
#count = t(count)
#symbol = bitr(count[,1],'ENSEMBL','SYMBOL','org.Hs.eg.db',drop = T)
#symbol = symbol[!duplicated(symbol$ENSEMBL),] 
#colnames(symbol)[1] = 'V1'
#count = left_join(count,symbol,'V1')

#count = count[!duplicated(count[,2]),]
#count = count[!is.na(count[,2]),]
rownames(count) = paste0(count[,1])
#rownames(tag) = paste0('cell_',tag[,1])

#for (k in 1:nrow(count)) {  rownames(count)[k] = strsplit(count[k,1],'__')[[1]][1]}
#meta = count[,2] %>% as.data.frame()
#colnames(meta) = 'Cluster_orig'
count[,c(1)] = NULL
#for (k in 1:ncol(count)){colnames(count)[k] = strsplit(colnames(count)[k],'\\|')[[1]][1]}
#colnames(count) = gsub('_membrane','',colnames(count)) %>% gsub('_secreted','',.) %>% gsub('_refseq','',.)
#count = t(count)
#rownames(meta) = colnames(count)
#rownames(count) = count[,'SYMBOL']
#colnames(count) = paste0('cell_',count[1,])
#count = count[-1,]
#count[, c('V1','SYMBOL')] = NULL
#count = count[,colnames(count)[!duplicated(colnames(count))]]
#count = t(count)
#colnames(count) = paste0('cell_',colnames(count))
seu = CreateSeuratObject(count)
#seu = AddMetaData(seu,meta)
#seu = AddMetaData(seu,meta)
seu.list[[samples[i]]] = seu
print(samples[i])
}

saveRDS(seu.list,file = paste0(dataset,'seu/','seu.list.rds'))

rm(list = ls())



# get samplename ----------------------------------------------------------
dataset = c('~/Project/PanCancerAtlas/data/AML_2/')
filea = list.files(dataset,recursive = F,full.names = F,pattern = 'GSM') %>% grep('barcodes',.,value = T)
samples =  stringr::str_sub(filea,12,-17) %>% as.data.frame(.)
write.table(samples,file = '~/samples.txt',col.names = F,row.names = F,quote = F,sep = '\t')

