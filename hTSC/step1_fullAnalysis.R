source('/home/zhangfeng/project/EVOL/source/MultiTools.R')

setwd('/home/zhangfeng/project/hTSC/data/hTSC_paper_RNA_seq/htseq_count')

MAT=read.table('MAT_FULL.NAME.txt',header=T,row.names=1,sep='\t')

MAT[1:5,1:5]


library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()



pbmc <- CreateSeuratObject(counts = MAT,  project = "HairSkin", min.cells = 0, min.features = 0)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 30)
#pbmc <- RunUMAP(pbmc, dims = 1:10,min.dist=0.3,n.neighbors=30)



DimPlot(pbmc, reduction = "pca")

