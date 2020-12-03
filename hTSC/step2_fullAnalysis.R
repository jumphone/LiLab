
library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()
source('/home/zhangfeng/project/EVOL/source/MultiTools.R')
setwd('/home/zhangfeng/project/hTSC/data/hTSC_paper_RNA_seq/htseq_count')

###
#Load data
MAT=read.table('MAT_FULL.NAME.txt',header=T,row.names=1,sep='\t')
MAT[1:5,1:5]


###
#Select protein coding genes
PCGENE=read.table('/home/database/annotation/hg19/Homo_sapiens.GRCh37.75.chr.pc.gene.SYM.bed',sep='\t',header=F)
MAT=MAT[which(rownames(MAT) %in% PCGENE[,4]),]
MAT[1:5,1:5]

##############
GRP=as.character(read.table('https://raw.githubusercontent.com/jumphone/BEER/master/SUP/KEGG_Ribosome.txt',sep='\t')[,1])
GCC1=as.character(read.table('https://raw.githubusercontent.com/jumphone/BEER/master/SUP/G1S.txt',sep='\t')[,1])
GCC2=as.character(read.table('https://raw.githubusercontent.com/jumphone/BEER/master/SUP/G2M.txt',sep='\t')[,1])

##############
pbmc <- CreateSeuratObject(counts = MAT,  project = "hTSC", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

##
#Use All Genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features =all.genes)
pbmc <- RunPCA(pbmc, features = all.genes, npcs = 10)
DimPlot(pbmc, reduction = "pca", pt.size=3) #+ NoLegend()

##
#Use Variable Genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
VariableFeaturePlot(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 10)
DimPlot(pbmc, reduction = "pca", pt.size=3) #+ NoLegend()


write.table(pbmc@meta.data, file='META.txt',row.names=T,col.names=T,quote=F,sep='\t')





###################
COUNT=as.matrix(pbmc@assays$RNA@counts)
SUM=apply(COUNT, 2,sum)
SRP=apply(COUNT[which(rownames(COUNT) %in% GRP),],2,sum)/ SUM
SCC1=apply(COUNT[which(rownames(COUNT) %in% GCC1),],2,sum)/ SUM
SCC2=apply(COUNT[which(rownames(COUNT) %in% GCC2),],2,sum)/ SUM
pbmc[['cc1']]=SCC1
pbmc[['cc2']]=SCC2
pbmc[['rp']]=SRP
######################

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#VariableFeaturePlot(pbmc)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc), vars.to.regress=c('cc1','cc2','rp'))

#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
#VariableFeaturePlot(pbmc)
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

#ElbowPlot(pbmc, ndims=50)
#pbmc <- RunUMAP(pbmc, dims = 1:10,min.dist=0.3,n.neighbors=15)
#DimPlot(pbmc, reduction = "umap", pt.size=3)

DimPlot(pbmc, reduction = "pca",dims=c(1,2), pt.size=3)
saveRDS(pbmc, file='pbmc.rds')
###############################################################################


DimPlot(pbmc, reduction = "pca", pt.size=3) #+ NoLegend()







































































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
