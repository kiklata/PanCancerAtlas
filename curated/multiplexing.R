#filter cell with no RNA or CMO count!!
seu.list = subset(seu.list, nCount_RNA >0 | nCount_CMO >0)

# ref: 10xgenomics.com/resources/analysis-guides/bioinformatics-tools-for-sample-demultiplexing
#      Nature Methods doi: 10.1038/s41592-019-0433-8
#DefaultAssay(seu.list) = 'CMO'
seu.list = NormalizeData(seu.list,assay = 'CMO',normalization.method = 'CLR')
seu.list = HTODemux(seu.list,assay = 'CMO')
seu.list = MULTIseqDemux(seu.list,assay = 'CMO')
# de-multiplexing parameters need refine?
#read CMO condition file, and subset exist combination of CMO+library
