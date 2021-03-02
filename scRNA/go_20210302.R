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


G1=as.character(read.table('geneset/SUB_GO_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT.txt',sep='\t',row.names=NULL,header=T)[,1])
G2=as.character(read.table('geneset/SUB_GO_CONNECTIVE_TISSUE_DEVELOPMENT.txt',sep='\t',row.names=NULL,header=T)[,1])
G3=as.character(read.table('geneset/SUB_GO_COLLAGEN_TRIMER.txt',sep='\t',row.names=NULL,header=T)[,1])
G4=as.character(read.table('geneset/SUB_GO_EPIDERMIS_DEVELOPMENT.txt',sep='\t',row.names=NULL,header=T)[,1])
G5=as.character(read.table('geneset/SUB_GO_CELL_CELL_JUNCTION.txt',sep='\t',row.names=NULL,header=T)[,1])
G6=as.character(read.table('geneset/SUB_GO_SKIN_DEVELOPMENT.txt',sep='\t',row.names=NULL,header=T)[,1])


#G1=as.character(read.table('geneset/GO_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT.txt',sep='\t',row.names=NULL,header=T)[,1])
#G2=as.character(read.table('geneset/GO_CONNECTIVE_TISSUE_DEVELOPMENT.txt',sep='\t',row.names=NULL,header=T)[,1])
#G3=as.character(read.table('geneset/GO_COLLAGEN_TRIMER.txt',sep='\t',row.names=NULL,header=T)[,1])
#G4=as.character(read.table('geneset/GO_EPIDERMIS_DEVELOPMENT.txt',sep='\t',row.names=NULL,header=T)[,1])
#G5=as.character(read.table('geneset/GO_CELL_CELL_JUNCTION.txt',sep='\t',row.names=NULL,header=T)[,1])
#G6=as.character(read.table('geneset/GO_SKIN_DEVELOPMENT.txt',sep='\t',row.names=NULL,header=T)[,1])


SM=function(x){
   return(smooth.spline(x,df=100)$y)
   }

###################################################


mat=MAT[which(rownames(MAT)%in% G1),]
o.mat=t(apply(mat,1,scale))
rownames(o.mat)=rownames(mat)
colnames(o.mat)=colnames(mat)
G1.S=apply(o.mat,2,mean)

mat=MAT[which(rownames(MAT)%in% G2),]
o.mat=t(apply(mat,1,scale))
rownames(o.mat)=rownames(mat)
colnames(o.mat)=colnames(mat)
G2.S=apply(o.mat,2,mean)

mat=MAT[which(rownames(MAT)%in% G3),]
o.mat=t(apply(mat,1,scale))
rownames(o.mat)=rownames(mat)
colnames(o.mat)=colnames(mat)
G3.S=apply(o.mat,2,mean)

mat=MAT[which(rownames(MAT)%in% G4),]
o.mat=t(apply(mat,1,scale))
rownames(o.mat)=rownames(mat)
colnames(o.mat)=colnames(mat)
G4.S=apply(o.mat,2,mean)

mat=MAT[which(rownames(MAT)%in% G5),]
o.mat=t(apply(mat,1,scale))
rownames(o.mat)=rownames(mat)
colnames(o.mat)=colnames(mat)
G5.S=apply(o.mat,2,mean)

mat=MAT[which(rownames(MAT)%in% G6),]
o.mat=t(apply(mat,1,scale))
rownames(o.mat)=rownames(mat)
colnames(o.mat)=colnames(mat)
G6.S=apply(o.mat,2,mean)




s.mat=t(cbind(G1.S,G2.S,G3.S,G4.S,G5.S,G6.S))
o.mat=t(apply(s.mat,1,scale))
o.mat=t(apply(s.mat,1,SM))
rownames(o.mat)=rownames(s.mat)
colnames(o.mat)=colnames(s.mat)

diff=apply(o.mat[,which(TAG=='TACSTD2')],1,mean)-apply(o.mat[,which(TAG=='DLK1')],1,mean)
o.mat=o.mat[order(diff),]


library('ComplexHeatmap')
library('circlize')
library('seriation')

color_fun_3 =colorRamp2(c(-0.06,-0.01,0,0.01,0.06 ), c('royalblue3','white','white','white','indianred3'))
#color_fun_3 =colorRamp2(c(-0.05,-0.005,0,0.005,0.05 ), c('royalblue3','white','white','white','indianred3'))

#color_fun_3 =colorRamp2(c(-1.5,-0.5,0,0.5,1.5 ), c('royalblue3','white','white','white','indianred3'))

############################################
ha_top = HeatmapAnnotation(  
     GROUP = TAG,
     col=list(GROUP=c('TACSTD2'='blue','DLK1'='red'))
     )

Heatmap(o.mat,row_title='',name="Score",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=FALSE, show_row_names=FALSE,
	      col=color_fun_3, border = TRUE,
        top_annotation = ha_top
	      )







#################################

mat=MAT[which(rownames(MAT)%in% G6),]
o.mat=t(apply(mat,1,scale))
o.mat=t(apply(o.mat,1,SM))
rownames(o.mat)=rownames(mat)
colnames(o.mat)=colnames(mat)

diff=apply(o.mat[,which(TAG=='TACSTD2')],1,mean)-apply(o.mat[,which(TAG=='DLK1')],1,mean)
o.mat=o.mat[order(diff),]

library('ComplexHeatmap')
library('circlize')
library('seriation')

color_fun_3 =colorRamp2(c(-0.5,-0.1,0,0.1,0.5 ), c('royalblue3','white','white','white','indianred3'))
#color_fun_3 =colorRamp2(c(-1.5,-0.5,0,0.5,1.5 ), c('royalblue3','white','white','white','indianred3'))

############################################
ha_top = HeatmapAnnotation(  
     GROUP = TAG,
     col=list(GROUP=c('TACSTD2'='blue','DLK1'='red'))
     )

Heatmap(o.mat,row_title='',name="C",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=FALSE, show_row_names=FALSE,
	      col=color_fun_3, border = TRUE,
        top_annotation = ha_top
	      )



























