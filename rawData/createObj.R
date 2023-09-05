source("~/Project/PanCancerAtlas/munge/rawData/genecorrect.R")
new.path = '~/Project/PanCancerAtlas/cache/'

human.map = HGNChelper::getCurrentHumanMap()

library(dplyr)

# correct gene symbol ---------------------------------------------

# seu --------------

data = 'BRCA_25'

seu.list@assays[["RNA"]]@data@Dimnames[[2]] %>% stringr::str_sub(.,18) %>% table
seu.list$sample_name = seu.list@assays[["RNA"]]@data@Dimnames[[2]] %>% stringr::str_sub(.,1,-18)

seu.list$sample_name = paste0(seu.list$sample_id)
seu.list$sample_name = Idents(seu.list)
seu.list$sample_name = data

seu.list$sample_name = stringr::str_sub(seu.list$Samples,1,-17)

# seu.list------------

for (i in 1:length(seu.list)) {
  #names(seu.list)[i] = seu.list[[i]]$orig.ident %>% unique() %>% as.character()
  seu.list[[i]]$sample_name = paste0(names(seu.list)[i])
  
  #seu.list[[i]]$sample_name = seu.list[[i]]$CellID %>% stringr::str_sub(.,18)
  #seu.list[[i]]$sample_name = paste0(names(seu.list)[i],'_',seu.list[[i]]@assays[["RNA"]]@data@Dimnames[[2]] %>% substring(.,18))
}

#sample_name = character()
#for (i in 1:ncol(seu.list)) {
#  sample_name[i] = paste0(strsplit(seu.list@assays[["RNA"]]@data@Dimnames[[2]][i],split = '\\.')[[1]][1],"_",
#                                     strsplit(seu.list@assays[["RNA"]]@data@Dimnames[[2]][i],split = '\\.')[[1]][2])}

seu.list = merge(seu.list[[1]],seu.list[2:length(seu.list)],add.cell.ids = names(seu.list)) 

seu.list = seu.list %>% gene.correct(., map = human.map)

saveRDS(seu.list,file = paste0(new.path,data,"/seu.list.rds"))


# ENSG --------------------------------------------------------------------
convertENSG = function(count, meta) {
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  
  symbols = bitr(rownames(count),
                 fromType = 'ENSEMBL',
                 'SYMBOL',
                 OrgDb = 'org.Hs.eg.db')
  
  count = as.data.frame(count %>% as.matrix())
  
  count$ENSEMBL = rownames(count)
  
  library(data.table)
  
  count = left_join(count, symbols, by = 'ENSEMBL')
  
  count = count[!is.na(count[, 'SYMBOL']), ]
  
  count[, 'ENSEMBL'] = NULL
  
  dup.gene = count[, 'SYMBOL'] %>% duplicated(.) %>% count[, 'SYMBOL'][.] %>% unique(.)
  
  new.mat = count[count[, 'SYMBOL'] %in% dup.gene, ] %>% as.data.table()
  
  count = count[!count[, 'SYMBOL'] %in% dup.gene, ]
  rownames(count) = count[, 'SYMBOL']
  count[, 'SYMBOL'] = NULL
  
  sub.sum.list = list()
  
  for (i in 1:length(dup.gene)) {
    sub.new.mat = new.mat[SYMBOL == dup.gene[i]] %>% .[, SYMBOL := NULL] %>% transpose()
    sub.sum.list[[dup.gene[i]]] = sub.new.mat[, sums := rowSums(.SD)] %>% .[, 'sums'] %>% transpose()
  }
  sub.new.mat.sum = rbindlist(sub.sum.list, use.names = T)
  
  colnames(sub.new.mat.sum) = colnames(count)
  sub.new.mat.sum[, 'SYMBOL'] = dup.gene
  sub.new.mat.sum = as.matrix(sub.new.mat.sum, rownames = 'SYMBOL')
  count = rbind(count, sub.new.mat.sum)
  
  library(Seurat)
  seu.list = CreateSeuratObject(counts = count, meta.data = meta)
  
  return(seu.list)
}

seu.list.new = convertENSG(count = seu.list@assays$RNA@counts, meta = seu.list@meta.data)



# RMS_1 -------------------------------------------------------------------

path = ('~/Project/PanCancerAtlas/data/RMS_1/seu')
all.run = list.dirs(path,recursive = F)

convert.seu.list = list()

for (i in all.run) {
  sub.files = list.files(path = i,pattern = '_filtered')
  for (k in sub.files) {
    obj = readRDS(file = paste0(i,'/',k))
    
    count = assay(obj,'counts')
    splice = assay(obj,'spliced')
    meta = colData(obj)@listData %>% as.data.frame()
    rownames(meta) = meta$barcodes
    
    seu = as.Seurat(obj,counts = 'counts',data = 'counts')
    seu[['spliced']] = CreateAssayObject(counts = splice)
    seu = AddMetaData(seu,metadata = meta)
    
    convert.seu.list[[paste0(i %>% gsub('/home/zhepan/Project/PanCancerAtlas/data/RMS_1/seu/','',.) ,'_',k %>% gsub('_filtered.rds','',.))]] = convertENSG(count = seu@assays$originalexp@counts, meta = seu@meta.data)
    
  }
}

seu.list = merge(convert.seu.list[[1]],convert.seu.list[2:length(convert.seu.list)],add.cell.ids = names(convert.seu.list)) 
seu.list = seu.list %>% gene.correct(., map = human.map)

