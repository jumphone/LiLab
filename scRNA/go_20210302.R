setwd('/home/database/data/backup_20201008_before/backup_20200818/scRNA_rdsFiles')

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()


library(Seurat)
pbmc=readRDS('d7_scRNA.rds')

FeaturePlot(pbmc, features=c('EPCAM','TACSTD2','DLK1'))

EXP=as.matrix(pbmc@assays$RNA@data)

EPCAM=rep(0,ncol(EXP))
TACSTD2=rep(0,ncol(EXP))
DLK1=rep(0,ncol(EXP))

EPCAM[which(EXP[which(rownames(EXP)=='EPCAM'),]>0)]=1
TACSTD2[which(EXP[which(rownames(EXP)=='TACSTD2'),]>0)]=1
DLK1[which(EXP[which(rownames(EXP)=='DLK1'),]>0)]=1

pbmc@meta.data$EPCAM=EPCAM
pbmc@meta.data$TACSTD2=TACSTD2
pbmc@meta.data$DLK1=DLK1

TACSTD2_GROUP=rep(0,ncol(EXP))
TACSTD2_GROUP[which(pbmc@meta.data$EPCAM>0 & pbmc@meta.data$TACSTD2>0 & pbmc@meta.data$DLK1==0)]=1

DLK1_GROUP=rep(0,ncol(EXP))
DLK1_GROUP[which(pbmc@meta.data$EPCAM>0 & pbmc@meta.data$TACSTD2==0 & pbmc@meta.data$DLK1>0)]=1

pbmc@meta.data$TACSTD2_GROUP=TACSTD2_GROUP
pbmc@meta.data$DLK1_GROUP=DLK1_GROUP


FeaturePlot(pbmc, c('TACSTD2_GROUP','DLK1_GROUP'))


CLUSTER=rep(0,ncol(pbmc))
CLUSTER[which(TACSTD2_GROUP>0)]=1
CLUSTER[which(DLK1_GROUP>0)]=2

pbmc@meta.data$cluster=CLUSTER

Idents(pbmc)=CLUSTER

TACSTD2.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = 2, only.pos = TRUE, min.pct = 0.25)
DLK1.markers <- FindMarkers(pbmc, ident.1 = 2, ident.2 = 1, only.pos = TRUE, min.pct = 0.25)

saveRDS(TACSTD2.markers,'TACSTD2.markers.rds')
saveRDS(DLK1.markers,'DLK1.markers.rds')

write.table(TACSTD2.markers,file='TACSTD2.markers.txt',sep='\t',row.names=T,col.names=T,quote=F)
write.table(DLK1.markers,file='DLK1.markers.txt',sep='\t',row.names=T,col.names=T,quote=F)



MAT1=EXP[,which(CLUSTER==1)]
MAT2=EXP[,which(CLUSTER==2)]

MAT=cbind(MAT1,MAT2)
TAG=c(rep('TACSTD2',ncol(MAT1)),rep('DLK1',ncol(MAT2)))












