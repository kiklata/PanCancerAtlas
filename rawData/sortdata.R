library(rlist)

path = '~/Project/PanCancerAtlas/cache'
dataset = list.dirs(path,full.names = F,recursive = F)
info.list = list()

for (i in dataset) {
    tmp = list.files(path = file.path(path,i),include.dirs = T,recursive = T)
    #print(tmp)
    info.list[[i]] = list(obj = tmp %>% grep('seu',.,value = T),
                          meta = tmp %>% grep('meta',.,value = T),
                          vdj = tmp %>% grep('VDJ',.,value = T))
}

list.filter(info.list,length(meta) != 0) %>% names
list.filter(info.list,length(vdj) != 0) %>% names
