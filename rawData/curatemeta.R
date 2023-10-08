library(Seurat)
library(dplyr)

dataset = 'HNSC_11'

meta.paper = read.delim(paste0("~/Project/PanCancerAtlas/cache/",dataset,"/meta/","metadata.txt.gz"))

meta.paper = fread('~/Project/PanCancerAtlas/cache/NPC_1/meta/GSE150430_npc_scRNA_EBV_genome_processed_data.txt.gz',)
#meta.paper = meta.paper.list[[2]]
colnames(meta.paper)

paper.col = c('Source.Name','Comment.ENA_SAMPLE.','Characteristics.age.',
              'Characteristics.sex.','Characteristics.individual.',
              'Characteristics.disease.','Characteristics.organism.part.','Factor.Value.disease.')          

meta.paper = meta.paper[,colnames(meta.paper) %in% paper.col]

saveRDS(meta.paper,file = paste0("~/Project/PanCancerAtlas/cache/",dataset,'/sample.meta.paper.rds'))

gs <- as.list(sample(10:100, size=100, replace=TRUE))
## sample gene sets
gs <- lapply(gs, function(n, p)
  paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))



meta.list = list()
all.meta = list.files(path = paste0("~/Project/PanCancerAtlas/cache/",dataset,"/meta/"),pattern = "")
paper.col = c('ID_paper','Cell_barcode','Technology','Cell_type_granular_mouse_correlations',
              'Cell_type_mouse_correlations', 'Cell_type_consensus_Jessa2022','Malignant_normal_consensus_Jessa2022')          

for (i in all.meta[1:11]) {
  
meta.paper = read.delim(paste0("~/Project/PanCancerAtlas/cache/",dataset,"/meta/",i),skip = 7,sep = ',')
sam.name = stringr::str_sub(i,11,-25)
#meta.paper = meta.paper[,colnames(meta.paper) %in% paper.col]
meta.paper$sample = sam.name

meta.list[[sam.name]] = meta.paper
}

meta = rbindlist(meta.list,use.names = T)
meta.paper = meta


meta.paper = read.delim(paste0("~/Project/PanCancerAtlas/cache/",dataset,"/meta/",all.meta[12]),sep = ',')
