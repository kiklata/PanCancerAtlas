## general.sample.meta
setwd('Project/PanCancerAtlas/')

info.list = readRDS('docs/info.list.rds')

meta.path = c('docs/samplemeta')
obj.path = c('data')

data = 'AML_7'
info.list[[data]]

# obj_only ----------------------------------------------------------------

obj = readRDS(file.path(obj.path,data,'seu.list.rds'))
#meta.paper = readRDS(file.path(obj.path,data,'meta','meta.paper.rds'))
meta.generate = read.delim(file.path(meta.path,paste0(data,'.txt')))

samples = obj$sample_name %>% table %>% names
samples
n_sample = samples %>% length

meta.generate[2:n_sample,] = NA

meta.generate$SampleID = paste(data,samples,sep = '.')
meta.generate$StudyID = data
meta.generate$PatientID = paste(data,c('UPN29','UPN23'),sep = '.')
#meta.generate$PatientID = stringr::str_sub(meta.generate$SampleID,1,-2)
meta.generate$Kit = 'sc'
meta.generate$Platform = '10x3v3'
meta.generate$Format = 'count'
meta.generate$Type = c('cancer')
meta.generate$Site = c('PBMC')
meta.generate$Metastasis = 'N'
meta.generate$Treatment = c('naive')
#meta.generate$Treatment = if_else(meta.generate$SampleID %>% do::right(.,1) == 'D','naive','treatment' )
meta.generate$Timepoint = c('pre')
meta.generate$Timepoint = if_else(meta.generate$Treatment == 'naive','pre','post')
meta.generate$Sort = 'none'

write.table(meta.generate,file = file.path(obj.path,data,'meta.generate.txt'),
            col.names = T,row.names = F,quote = F)



meta.read = data.frame(SampleID = paste(data,samples,sep = '.'), Age = NA, Sex = NA)
meta.read$Age = c(47,48)
meta.read$Sex = c('F','F')  

write.table(meta.read,file = file.path(obj.path,data,'meta.read.txt'),
            col.names = T,row.names = F,quote = F)
