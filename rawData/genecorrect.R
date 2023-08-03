gene.correct = function(obj,check.dup = FALSE) {
  library(dplyr)
  library(Seurat)
  library(HGNChelper)
  library(cli)
  library(data.table)
  
  # This code checks the gene symbols in the 'exp.gene' variable against the HGNC database to ensure they are valid human gene symbols.
  # If any symbols are not valid, they are filtered out and their suggested symbols are used instead.
  # The resulting suggested symbols are then used to replace the original gene symbols in the 'exp.gene' variable.
  human.map = getCurrentHumanMap()
  
  exp.gene = rownames(obj)
  exp.gene = stringr::str_replace(exp.gene,
                                  pattern = c("(\\.)[0-9+]"),
                                  replacement = '')
  
  hgnc.check = checkGeneSymbols(exp.gene, species = 'human', map = human.map)
  
  remove.gene = hgnc.check %>% filter(., is.na(.$Suggested.Symbol) == T) %>% .$x
  
  trans.gene = filter(hgnc.check, Approved == 'FALSE') %>% filter(., Suggested.Symbol != '') # all symbols need to be converted
  
  for (i in 1:nrow(trans.gene)) {
    trans.gene$Suggested.Symbol[i] = strsplit(trans.gene$Suggested.Symbol[i], ' /// ')[[1]][1]
    
  }
  # This code block extracts the original gene symbols from a matrix based on whether they exist or not in a list of gene symbols.
  # It then converts the extracted symbols to their corresponding suggested symbols using a translation matrix.
  # The resulting converted symbols are stored in two separate vectors: previous.exist.gene.new and previous.non.gene.new.
  # The original symbols that were extracted are also stored in two separate vectors: previous.exist.gene and previous.non.gene.
  
  previous.exist.gene = trans.gene$x[trans.gene$Suggested.Symbol %in% exp.gene] # neo-symbols is existed in the matrix, extract the org-symbols
  previous.non.gene = trans.gene$x[!trans.gene$Suggested.Symbol %in% exp.gene] # neo-symbols is not existed in the matrix, extract the org-symbols
  
  previous.exist.gene.new = trans.gene$Suggested.Symbol[trans.gene$Suggested.Symbol %in% exp.gene] # converted pre-existed symbols
  previous.non.gene.new = trans.gene$Suggested.Symbol[!trans.gene$Suggested.Symbol %in% exp.gene] # converted non-existed symbols
  
  # This function corrects gene symbols in the RNA-seq count matrix. It takes in a list of previous gene symbols and a list of new gene symbols,
  # and replaces the old symbols with the new ones in the count matrix. If there are duplicated symbols, it combines the counts for those symbols.
  # The function returns a new count matrix with corrected gene symbols.
  count = obj@assays$RNA@counts
  
  row.n = nrow(count)
  col.n = ncol(count)
  
  cli_alert_info("Detect a {col.n} cells by {row.n} genes matrix")
  
  rownames(count) = exp.gene
  count = count[!rownames(count) %in% remove.gene, ]
  
  # deal with original duplicated symbols-----------------------
  if(check.dup == TRUE){
  all.dup.gene = rownames(count)[rownames(count) %>% duplicated()] %>% unique()
  
  if (length(all.dup.gene) > 0) {
    time.start = Sys.time()
    
    new.mat = count[all.dup.gene, ] %>% as.matrix() %>% as.data.table()
    new.mat[, 'symbol'] = all.dup.gene
    count = count[!rownames(count) %in% all.dup.gene,]
    
    sub.new.mat.sum.list = list()
    for (i in 1:length(all.dup.gene)) {
      sub.new.mat = new.mat[symbol == all.dup.gene[i]] %>% .[, symbol := NULL] %>% transpose() %>% as.matrix()
      #sub.new.mat = sub.new.mat[, lapply(.SD,sum), by = symbol]
      sub.new.mat.sum.list[[all.dup.gene[i]]] = apply(sub.new.mat, 1, sum) %>% as.data.table() %>% transpose()
      
    }
    #new.mat = as.matrix(new.mat[, lapply(.SD,sum), by = symbol])
    sub.new.mat.sum = rbindlist(sub.new.mat.sum.list, use.names = T)
    colnames(sub.new.mat.sum) = colnames(count)
    sub.new.mat.sum[, 'symbol'] = all.dup.gene
    #sub.new.mat.sum = as.matrix(sub.new.mat.sum)
    sub.new.mat.sum = as.matrix(sub.new.mat.sum, rownames = 'symbol')
    #rownames(sub.new.mat.sum) = all.dup.gene
    #rownames(new.mat) = new.mat[,'symbol']
    #sub.new.mat.sum[,'symbol'] = NULL
    sub.new.mat.sum.predup = sub.new.mat.sum
    
    count = rbind(count, sub.new.mat.sum.predup)
    
    time.end = Sys.time()
    tbld = difftime(time.end, time.start, units = 'auto')[[1]] %>% round(., 1)
    tbld_unit = difftime(time.end, time.start, units = 'auto') %>% units
    
    cli_alert_success(paste0(
      "Converted ",
      length(all.dup.gene),
      " original duplicated symbols in {tbld} {tbld_unit}"
    ))
    
    
  } else {
    NULL
  }
  }else{NULL}
  
  if (length(previous.exist.gene) == 0) {
    NULL
  } else {
    # previous.exist.gene.new have duplicated symbol, multi old symbol to one symbol---------------
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
    
    all.exist.gene = c(previous.exist.gene, previous.exist.gene.new %>% unique())
    all.exist.gene.correct = c(previous.exist.gene.new,
                               previous.exist.gene.new %>% unique())
    
    new.mat = count[all.exist.gene, ] %>% as.matrix() %>% as.data.table()
    
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
    
    count = count[!rownames(count) %in% all.exist.gene,]
    
    #new.mat = data.table::rbindlist(list(tmp_mat_list)) %>% transpose()
    new.mat[, 'symbol'] = all.exist.gene.correct
    all.dup.gene = all.exist.gene.correct[all.exist.gene.correct %>% duplicated()] %>% unique()
    sub.new.mat.sum.list = list()
    for (i in 1:length(all.dup.gene)) {
      sub.new.mat = new.mat[symbol == all.dup.gene[i]] %>% .[, symbol := NULL] %>% transpose() %>% as.matrix()
      #sub.new.mat = sub.new.mat[, lapply(.SD,sum), by = symbol]
      sub.new.mat.sum.list[[all.dup.gene[i]]] = apply(sub.new.mat, 1, sum) %>% as.data.table() %>% transpose()
      
    }
    #new.mat = as.matrix(new.mat[, lapply(.SD,sum), by = symbol])
    sub.new.mat.sum = rbindlist(sub.new.mat.sum.list, use.names = T)
    colnames(sub.new.mat.sum) = colnames(count)
    sub.new.mat.sum[, 'symbol'] = all.dup.gene
    #sub.new.mat.sum = as.matrix(sub.new.mat.sum)
    sub.new.mat.sum = as.matrix(sub.new.mat.sum, rownames = 'symbol')
    #rownames(sub.new.mat.sum) = all.dup.gene
    #rownames(new.mat) = new.mat[,'symbol']
    #sub.new.mat.sum[,'symbol'] = NULL
    sub.new.mat.sum.exist = sub.new.mat.sum
    
    count = rbind(count, sub.new.mat.sum.exist)
    
    time.end = Sys.time()
    tbld = difftime(time.end, time.start, units = 'auto')[[1]] %>% round(., 1)
    tbld_unit = difftime(time.end, time.start, units = 'auto') %>% units
    
    cli_alert_success(paste0(
      "Converted ",
      length(previous.exist.gene),
      " existing symbols in {tbld} {tbld_unit}"
    ))
    
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
  new.mat = count.new[all.non.gene, ] %>% as.matrix() %>% as.data.table()
  
  count.new = count.new[!rownames(count.new) %in% all.non.gene,]
  
  
  new.mat[, 'symbol'] = all.non.gene.correct
  all.dup.gene = all.non.gene.correct[all.non.gene.correct %>% duplicated()] %>% unique()
  sub.new.mat.sum.list = list()
  for (i in 1:length(all.dup.gene)) {
    sub.new.mat = new.mat[symbol == all.dup.gene[i]] %>% .[, symbol := NULL] %>% transpose() %>% as.matrix()
    #sub.new.mat = sub.new.mat[, lapply(.SD,sum), by = symbol]
    sub.new.mat.sum.list[[all.dup.gene[i]]] = apply(sub.new.mat, 1, sum) %>% as.data.table() %>% transpose()
    
  }
  #new.mat = as.matrix(new.mat[, lapply(.SD,sum), by = symbol])
  sub.new.mat.sum = rbindlist(sub.new.mat.sum.list, use.names = T)
  colnames(sub.new.mat.sum) = colnames(count)
  sub.new.mat.sum[, 'symbol'] = all.dup.gene
  #rownames(sub.new.mat.sum) = all.dup.gene
  #sub.new.mat.sum = as.matrix(sub.new.mat.sum)
  #rownames(new.mat) = new.mat[,'symbol']
  #new.mat[,'symbol'] = NULL
  new.mat = new.mat[!symbol %in% all.dup.gene]
  #new.mat = data.table::rbindlist(list(tmp_mat_list)) %>% transpose()
  #new.mat = aggregate(.~symbol,new.mat,sum)
  new.mat = rbind(new.mat, sub.new.mat.sum) %>% as.matrix(., rownames = 'symbol')
  sub.new.mat.sum.nonexist = new.mat
  
  #rownames(new.mat) = all.non.gene.correct[!(all.non.gene.correct %in% unique(all.non.gene.correct[duplicated(all.non.gene.correct)]))]
  #new.mat[,'symbol'] = NULL
  count.new = rbind(count.new, sub.new.mat.sum.nonexist)
  
  time.end = Sys.time()
  tbld = difftime(time.end, time.start, units = 'auto')[[1]] %>% round(., 1)
  tbld_unit = difftime(time.end, time.start, units = 'auto') %>% units
  
  cli_alert_success(paste0(
    "Converted ",
    length(all.non.gene),
    " non-existing symbols in {tbld} {tbld_unit}"
  ))
  
  correct.obj = CreateSeuratObject(count.new,
                                   meta.data = obj@meta.data,
                                   min.cells = 3)
  
  if (length(obj@assays) == 1) {
    NULL
  } else {
    for (n in 2:length(obj@assays)) {
      omics = obj[[names(obj@assays)[n]]]
      correct.obj[[names(obj@assays)[n]]] = subset(omics, cells = colnames(correct.obj))
    }
  }
  return(correct.obj)
  
}
