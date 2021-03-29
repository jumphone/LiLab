
/home/toolkit/tools/R4.0.3/bin/R




setwd('/home/database/data/QW_LAB_data')


library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()



library(dplyr)
library(Seurat)
library(patchwork)


DATA=read.table('ATAC_CHROMVAR_MAT.txt',sep='\t',row.names=1,header=T)


MAT=DATA-min(DATA)
MAT=100*MAT

BATCH=c(rep('YZ6038',2),rep('YZ6080',2), rep('YZ5736',1),rep('YZ5834',5),rep('YZ5736',1),rep('YZ5834',4))

sample_name=colnames(DATA)

sample_tag=c('D0','D3','D0_S2KO','D3_S2KO','D0','D1.5','D2.5','D3.5','D3.5','D5','D0_Foxh1KO','D1.5_Foxh1KO','D2.5_Foxh1KO','D3.5_Foxh1KO','D5_Foxh1KO')


pbmc <- CreateSeuratObject(counts = MAT, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)

pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = all.genes,npcs=14)

pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag
pbmc@meta.data$batch=BATCH

DimPlot(pbmc, reduction = "pca",group.by='batch',label=T,pt.size=3)+NoLegend()


DimPlot(pbmc, reduction = "pca",group.by='stag',label=T,pt.size=3,dims=c(1,2))+NoLegend()




##############################











######################
source('/home/database/data/QW_LAB_data/BEER.R')

DATA.COMBAT=.combat(DATA,BATCH)


MAT=DATA.COMBAT-min(DATA.COMBAT)
MAT=100*MAT

BATCH=c(rep('YZ6038',2),rep('YZ6080',2), rep('YZ5736',1),rep('YZ5834',5),rep('YZ5736',1),rep('YZ5834',4))

sample_name=colnames(DATA)

sample_tag=c('D0','D3','D0_S2KO','D3_S2KO','D0','D1.5','D2.5','D3.5','D3.5','D5','D0_Foxh1KO','D1.5_Foxh1KO','D2.5_Foxh1KO','D3.5_Foxh1KO','D5_Foxh1KO')



pbmc <- CreateSeuratObject(counts = MAT, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


pbmc <- RunPCA(pbmc, features = all.genes,npcs=10)


pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag
pbmc@meta.data$batch=BATCH

DimPlot(pbmc, reduction = "pca",group.by='batch',label=T,pt.size=3)+NoLegend()


DimPlot(pbmc, reduction = "pca",group.by='stag',label=T,pt.size=3,dims=c(1,2))+NoLegend()





##############################

USED_INDEX=c(1,5,6,7,2,8,9,10)


USED_NAME=c( "D0_1"   ,"D0_2"   ,"D1.5" ,"D2.5" ,"D3"  , "D3.5_1", "D3.5_2", "D5")
USED_TAG=sample_tag[USED_INDEX]
USED_DATA=DATA.COMBAT[,USED_INDEX]

colnames(USED_DATA)=USED_NAME


################################



SMAT=t(apply(USED_DATA,1,scale))
rownames(SMAT)=rownames(USED_DATA)
colnames(SMAT)=USED_NAME




DIFF=cbind(abs(SMAT[,1]-SMAT[,2]),abs(SMAT[,6]-SMAT[,7]))

MEAN.DIFF=apply(DIFF,1,mean)
SMAT=SMAT[which(MEAN.DIFF<0.5),]




D0=apply(SMAT[,1:2],1,mean)
D1.5=SMAT[,3]
D2.5=SMAT[,4]
D3=SMAT[,5]
D3.5=apply(SMAT[,6:7],1,mean)
D5=SMAT[,8]


MEAN.MAT=cbind(D0,D1.5,D2.5,D3,D3.5,D5)


D=dist(MEAN.MAT)
H=hclust(D)




ha = rowAnnotation(foo = anno_mark(at = c( which(rownames(SMAT)=='MA0141.3_ESRRB'),
                                           which(rownames(SMAT)=='MA0668.1_NEUROD2'),
                                           which(rownames(SMAT)=='MA0662.1_MIXL1')
                                           
                                            ), 
                         labels = c('MA0141.3_ESRRB','MA0668.1_NEUROD2','MA0662.1_MIXL1')))




library(ComplexHeatmap)
library(circlize)
color_fun =colorRamp2(c(-1, 0, 1), c('royalblue3','white','indianred3'))


CN=3

par(mfrow=c(1,1))
C=cutree(H,k=CN)

Heatmap(SMAT,row_title='',name="Z",
        cluster_columns=F, cluster_rows=F,
	      show_column_dend = F, show_row_dend = F, right_annotation = ha,
	      show_column_names=T, show_row_names=F,
	      col=color_fun, border = TRUE, 
	      column_names_gp = gpar(fontsize = 8),
	      row_names_gp = gpar(fontsize = 8),
	      heatmap_width = unit(0.7, "npc"),
	      row_split = C
	      )



MEAN=c()
UPPER=c()
LOWER=c()


i=1
while(i<=CN){
    this.mean=apply(MEAN.MAT[which(C==i),],2,mean)
    this.upper=apply(MEAN.MAT[which(C==i),],2,quantile, 0.975)
    this.lower=apply(MEAN.MAT[which(C==i),],2,quantile, 0.025)
    MEAN=cbind(MEAN,this.mean)
    UPPER=cbind(UPPER,this.upper)
    LOWER=cbind(LOWER,this.lower)
    i=i+1
    }



par(mfrow=c(CN,1))
par(mar=c(0.5,1,0.5,1))
i=1
while(i<=CN){
    plot(MEAN[,i],ylim=c(-2,2),type='b',pch=16,col='black',axes=FALSE, frame.plot=TRUE,cex=2)
    #arrows(x0=c(1:5),y0=UPPER[,i],x1=c(1:5),y1=LOWER[,i],angle = 90,code = 3,length = 0.1)
    i=i+1
    }

table(C)


OUT=cbind(C[order(C)], SMAT[order(C),])
colnames(OUT)[1]='Cluster'

OUT=cbind(row.names(OUT),OUT)
colnames(OUT)[1]='TF'


write.table(OUT,'ATAC_CHROMVAR_CLUSTER_MAT.txt',row.names=F,col.names=T,sep='\t',quote=F)












#################################################


/home/toolkit/tools/R4.0.3/bin/R


setwd('/home/database/data/QW_LAB_data')


library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()



DATA=read.table('ATAC_CHROMVAR_CLUSTER_MAT_StepMiner.pcl',sep='\t',row.names=1,header=T)
DATA=DATA[,3:ncol(DATA)]


SMAT=as.matrix(DATA)


library(ComplexHeatmap)
library(circlize)
color_fun =colorRamp2(c(-1,0,1), c('royalblue3','white','indianred3'))




ha = rowAnnotation(foo = anno_mark(at = c( which(rownames(SMAT)=='MA0141.3_ESRRB'),
                                           which(rownames(SMAT)=='MA0668.1_NEUROD2'),
                                           which(rownames(SMAT)=='MA0662.1_MIXL1')
                                           
                                            ), 
                         labels = c('MA0141.3_ESRRB','MA0668.1_NEUROD2','MA0662.1_MIXL1')))




par(mfrow=c(1,1))
C=cutree(H,k=CN)

Heatmap(SMAT,row_title='',name="Z",
        cluster_columns=F, cluster_rows=F,
          show_column_dend = F, show_row_dend = F, right_annotation = ha,
          show_column_names=T, show_row_names=F,
          col=color_fun, border = TRUE, 
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 8),
          heatmap_width = unit(0.7, "npc")
          )




















#COR=cor(t(MEAN.MAT),method='spearman')




library(ComplexHeatmap)
library(circlize)
color_fun =colorRamp2(c(-1, 0, 1), c('royalblue3','white','indianred3'))

Heatmap(SMAT,row_title='',name="Z",
        cluster_columns=F, cluster_rows=T,
	      show_column_dend = F, show_row_dend = T, 
	      show_column_names=T, show_row_names=F,
	      col=color_fun, border = TRUE, 
	      column_names_gp = gpar(fontsize = 8),
	      row_names_gp = gpar(fontsize = 8),
	      heatmap_width = unit(0.7, "npc")
	      )
























































































RPKM=read.table('ATAC_RPKM_MAT.txt',sep='\t',row.names=1,header=T)
BATCH=c(rep('YZ6038',2),rep('YZ6080',2), rep('YZ5736',1),rep('YZ5834',5),rep('YZ5736',1),rep('YZ5834',4))


sample_name=colnames(RPKM)

sample_tag=c('D0','D3','D0_S2KO','D3_S2KO','D0','D1.5','D2.5','D3.5','D3.5','D5','D0_Foxh1KO','D1.5_Foxh1KO','D2.5_Foxh1KO','D3.5_Foxh1KO','D5_Foxh1KO')



##
#Without batch correction


pbmc <- CreateSeuratObject(counts = RPKM, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=14)


pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag
pbmc@meta.data$batch=BATCH

DimPlot(pbmc, reduction = "pca",group.by='batch',label=T,pt.size=3)+NoLegend()


DimPlot(pbmc, reduction = "pca",group.by='stag',label=T,pt.size=3,dims=c(1,4))+NoLegend()




pbmc <- RunUMAP(pbmc, dims = 1:10, n.neighbors=3, min.dist=0.3,umap.method='umap-learn')


DimPlot(pbmc, reduction = "umap",group.by='stag',label=T,pt.size=3,dims=c(1,2))+NoLegend()







#With batch correction



source('/home/zhangfeng/project/BEER4/BEER.R')

RPKM.COMBAT=.combat(RPKM,BATCH)
RPKM.COMBAT[which(RPKM.COMBAT<0)]=0



pbmc <- CreateSeuratObject(counts = RPKM.COMBAT, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=14)



pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag

pbmc@meta.data$batch=BATCH
DimPlot(pbmc, reduction = "pca",group.by='batch',label=F,pt.size=3)+NoLegend()
DimPlot(pbmc, reduction = "pca",group.by='stag',label=T,pt.size=3,dims=c(1,2))+NoLegend()




pbmc <- RunUMAP(pbmc, dims = 1:10, n.neighbors=5, min.dist=0.5,umap.method='umap-learn')
DimPlot(pbmc, reduction = "umap",group.by='stag',label=T,pt.size=3,dims=c(1,2))+NoLegend()






















































































COUNT=read.table('ATAC_PEAK_INTER_MERGED_COUNT_MAT.txt',sep='\t',row.names=1,header=T)
BATCH=c(rep('YZ6038',2),rep('YZ6080',2), rep('YZ5736',1),rep('YZ5834',5),rep('YZ5736',1),rep('YZ5834',4))


























##
#Without batch correction


pbmc <- CreateSeuratObject(counts = COUNT, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=14)


sample_name=colnames(COUNT)

sample_tag=c('D0','D3','D0_S2KO','D3_S2KO','D0','D1.5','D2.5','D3.5','D3.5','D5','D0_Foxh1KO','D1.5_Foxh1KO','D2.5_Foxh1KO','D3.5_Foxh1KO','D5_Foxh1KO')

pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag
pbmc@meta.data$batch=BATCH

DimPlot(pbmc, reduction = "pca",group.by='batch',label=T,pt.size=3)+NoLegend()


DimPlot(pbmc, reduction = "pca",group.by='stag', label=T,pt.size=3,dims=c(4,3))+NoLegend()



#With batch correction



source('/home/zhangfeng/project/BEER4/BEER.R')


pbmc <- CreateSeuratObject(counts = COUNT, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

DATA=pbmc@assays$RNA@data
DATA.COMBAT=.combat(DATA,BATCH)
DATA.COMBAT[which(DATA.COMBAT<0)]=0
pbmc@assays$RNA@data=DATA.COMBAT

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 20000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=14)


sample_name=colnames(COUNT)

sample_tag=c('D0','D3','D0_S2KO','D3_S2KO','D0','D1.5','D2.5','D3.5','D3.5','D5','D0_Foxh1KO','D1.5_Foxh1KO','D2.5_Foxh1KO','D3.5_Foxh1KO','D5_Foxh1KO')



PCA=pbmc@reductions$pca@cell.embeddings

TIME=c(0,3,0,3,0,1.5,2.5,3.5,3.5,5,0,1.5,2.5,3.5,5)

TCOR=apply(PCA,2,cor,TIME, method='spearman')
TCOR


pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag

pbmc@meta.data$batch=BATCH
DimPlot(pbmc, reduction = "pca",group.by='batch',label=F,pt.size=3)+NoLegend()


DimPlot(pbmc, reduction = "pca",group.by='stag',label=T,pt.size=3,dims=c(9,2))+NoLegend()


pbmc <- RunUMAP(pbmc, dims = 1:10, n.neighbors = 3)


saveRDS(pbmc,'ATAC_seurat.rds')
