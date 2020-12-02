source('/home/zhangfeng/project/EVOL/source/MultiTools.R')

setwd('/home/zhangfeng/project/hTSC/data/hTSC_paper_RNA_seq/htseq_count')

MAT=read.table('MAT_FULL.NAME.txt',header=T,row.names=1,sep='\t')

MAT[1:5,1:5]



GRP=as.character(read.table('https://raw.githubusercontent.com/jumphone/BEER/master/SUP/KEGG_Ribosome.txt',sep='\t')[,1])
GCC=as.character(read.table('https://raw.githubusercontent.com/jumphone/BEER/master/SUP/KEGG_CellCycle.txt',sep='\t')[,1])




library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()



pbmc <- CreateSeuratObject(counts = MAT,  project = "hTSC", min.cells = 0, min.features = 0)
###################
COUNT=as.matrix(pbmc@assays$RNA@counts)
SUM=apply(COUNT, 2,sum)
SRP=apply(COUNT[which(rownames(COUNT) %in% GRP),],2,sum)/ SUM
SCC=apply(COUNT[which(rownames(COUNT) %in% GCC),],2,sum)/ SUM
pbmc[['cc']]=SCC
pbmc[['rp']]=SRP
######################

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress=c('cc','rp'))

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#VariableFeaturePlot(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

DimPlot(pbmc, reduction = "pca")
saveRDS(pbmc, file='pbmc.rds')
###############################################################################


DimPlot(pbmc, reduction = "pca")+ NoLegend()







































































############
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress=c('cc','rp'))
#pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
#pbmc <- RunPCA(pbmc, features = all.genes, npcs = 50)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)
#ElbowPlot(pbmc, ndims=50)
#pbmc <- RunUMAP(pbmc, dims = 1:30,min.dist=0.3,n.neighbors=10)
#DimPlot(pbmc, reduction = "umap")
