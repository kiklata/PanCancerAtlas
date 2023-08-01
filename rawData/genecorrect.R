gene.correct = function(obj) {
  library(dplyr)
  library(Seurat)
  library(HGNChelper)
  library(cli)
  library(data.table)
  
  exp.gene = rownames(obj)
  
  hgnc.check = checkGeneSymbols(exp.gene, species = 'human')
  trans.gene = filter(hgnc.check, Approved == 'FALSE') %>% filter(., Suggested.Symbol != '') # all symbols need to be converted
  
  for (i in 1:nrow(trans.gene)) {
    trans.gene$Suggested.Symbol[i] = strsplit(trans.gene$Suggested.Symbol[i], ' /// ')[[1]][1]
    
  }
  
  previous.exist.gene = trans.gene$x[trans.gene$Suggested.Symbol %in% exp.gene] # neo-symbols is existed in the matrix, extract the org-symbols
  previous.non.gene = trans.gene$x[!trans.gene$Suggested.Symbol %in% exp.gene] # neo-symbols is not existed in the matrix, extract the org-symbols
  
  previous.exist.gene.new = trans.gene$Suggested.Symbol[trans.gene$Suggested.Symbol %in% exp.gene] # converted pre-existed symbols
  previous.non.gene.new = trans.gene$Suggested.Symbol[!trans.gene$Suggested.Symbol %in% exp.gene] # converted non-existed symbols
  
  if (length(previous.exist.gene) == 0) {
    NULL
  }else {
  # previous.exist.gene.new have duplicated symbol, multi old symbol to one symbol---------------
  count = obj@assays$RNA@counts
  
  #cli_progress_bar(
  #  format = paste0(
  #    "{pb_spin} Converting existing symbols [{pb_current}/{pb_total}]  ETA:{pb_eta}"
  #  ),
  #  total = length(previous.exist.gene),
  #  type = 'tasks',
  #  clear = T
  #)
  time.start = Sys.time()
  #tmp_mat_list = list()
  
  #value.new = count[previous.exist.gene.new %>% unique(), ] %>% as.data.table()
  #value.new[,'symbol'] = previous.exist.gene.new %>% unique()
  
  all.exist.gene = c(previous.exist.gene,previous.exist.gene.new %>% unique())
  all.exist.gene.correct = c(previous.exist.gene.new, previous.exist.gene.new %>% unique())
  
  new.mat = count[all.exist.gene,] %>% as.data.table()
  
  # extract pre-existed symbol matrix
  #for (i in  1:length(all.exist.gene)) {
    #value = count[all.exist.gene[i], ]
    #value = value.old + value.new
    #names(value) = NULL
    
    #count = count[rownames(count) != previous.exist.gene[i], ]
    #count = count[rownames(count) != previous.exist.gene.new[i], ] 
    
    #tmp_mat = matrix(
    #  value,
    #  nrow = 1,
    #  ncol = length(value),
    #  dimnames = list(all.exist.gene.correct[i], names(value))
    #)
    #tmp_mat_list[[i]] = tmp_mat
    #count = rbind(count, tmp_mat)
    #cli_progress_update()
    
  #}
  
  count = count[!rownames(count) %in% all.exist.gene, ]

  #new.mat = data.table::rbindlist(list(tmp_mat_list)) %>% transpose()
  new.mat[,'symbol'] = all.exist.gene.correct
  new.mat = aggregate(.~symbol,new.mat,sum)
  rownames(new.mat) = new.mat[,'symbol']
  new.mat[,'symbol'] = NULL
  
  count = rbind(count, new.mat %>% as.matrix())
  
  time.end = Sys.time()
  tbld = difftime(time.end,time.start,units = 'auto')[[1]] %>% round(.,1)
  tbld_unit = difftime(time.end,time.start,units = 'auto') %>% units
  
  cli_alert_success(paste0("Converted ", length(previous.exist.gene)," existing symbols in {tbld} {tbld_unit}"))

  }
  # previous.non.gene check multi to one gene-------------
  count.new = count
  
  #cli_progress_bar(
  #  format = paste0(
  #    "{pb_spin} Converting non-existing symbols [{pb_current}/{pb_total}]  ETA:{pb_eta}"
  #  ),
  #  total = length(previous.non.gene),
  #  type = 'tasks',
  #  clear = T
  #)
  
  time.start = Sys.time()
  
  all.non.gene = c(previous.non.gene)
  all.non.gene.correct = c(previous.non.gene.new)
  
  new.mat = count.new[all.non.gene,] %>% as.data.table()
  

  #for (i in 1:length(previous.non.gene)) {
  #  value.old = count.new[previous.non.gene[i], ]
  #  if (previous.non.gene.new[i] %in% rownames(count.new)) {
  #    value.new = count.new[previous.non.gene.new[i], ]
  #  } else{
  #    value.new = 0
  #  }
  #  value = value.old + value.new
    
    #count.new = count.new[rownames(count.new) != previous.non.gene[i], ]
    
    #if (previous.non.gene.new[i] %in% rownames(count.new)) {
    #  count.new = count.new[rownames(count.new) != previous.non.gene.new[i], ]
    #}
    
  #  tmp_mat = matrix(
  #    value,
  #    nrow = 1,
  #    ncol = length(value),
  #    dimnames = list(previous.non.gene.new[i], names(value))
  #  )
  #  tmp_mat_list[[i]] = tmp_mat
    
    
    #count.new = rbind(count.new, tmp_mat)

   # cli_progress_update()
    
  #}
  
  count.new = count.new[!rownames(count.new) %in% all.non.gene, ]
  
  #new.mat = data.table::rbindlist(list(tmp_mat_list)) %>% transpose()
  new.mat[,'symbol'] = all.non.gene.correct
  new.mat = aggregate(.~symbol,new.mat,sum)
  rownames(new.mat) = new.mat[,'symbol']
  new.mat[,'symbol'] = NULL
  
  count.new = rbind(count.new, new.mat %>% as.matrix())
  
  time.end = Sys.time()
  tbld = difftime(time.end,time.start,units = 'auto')[[1]] %>% round(.,1)
  tbld_unit = difftime(time.end,time.start,units = 'auto') %>% units
  
  cli_alert_success(paste0("Converted ", length(all.non.gene)," non-existing symbols in {tbld} {tbld_unit}"))
  
  correct.obj = CreateSeuratObject(count.new,
                                   meta.data = obj@meta.data,
                                   min.cells = 3)
  
  exp.gene = rownames(correct.obj)
  
  hgnc.check = checkGeneSymbols(exp.gene)
  left.gene = filter(hgnc.check, Approved == 'TRUE')$x
  
  correct.obj = correct.obj[left.gene, ]
  
  if (length(obj@assays)==1) {
    NULL
  }else {
    for (n in 2:length(obj@assays)) {
      omics = obj[[names(obj@assays)[n]]]
      correct.obj[[names(obj@assays)[n]]] = subset(omics,cells = colnames(correct.obj))
    }
  }
  return(correct.obj)
  
}
