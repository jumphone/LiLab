


/home/toolkit/tools/R4.0.3/bin/R


setwd('/home/database/data/QW_LAB_data')


library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()



library(dplyr)
library(Seurat)
library(patchwork)


RPKM=read.table('RNA_RPKM_MAT.txt',sep='\t',row.names=1,header=T)
BATCH=c(rep('P3733',10),rep('p08143',16))





##
#Without batch correction


pbmc <- CreateSeuratObject(counts = RPKM, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=25)


sample_name=colnames(RPKM)

sample_tag=c('D0','D0','D1','D1','D2','D2','D3','D3','D4','D4','D0','D4','D0','D4','D0_S2.3DKO','D4_S2.3DKO','D0_S2.3DKO','D4_S2.3DKO','D0_S4KO',
'D4_S4KO','D0_S4KO','D4_S4KO','D0_Foxh1KO','D4_Foxh1KO','D0_Foxh1KO','D4_Foxh1KO')


pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag
pbmc@meta.data$batch=BATCH

DimPlot(pbmc, reduction = "pca",group.by='batch',label=T,pt.size=3)+NoLegend()



#With batch correction



source('/home/zhangfeng/project/BEER4/BEER.R')

RPKM.COMBAT=.combat(RPKM,BATCH)
RPKM.COMBAT[which(RPKM.COMBAT<0)]=0


COL=rep('grey',26)
COL[c(15,16,17,18)]='red'
par(las=2)
par(mar=c(10,3,3,3))
barplot(t(RPKM.COMBAT)[,which(rownames(RPKM.COMBAT)=='Smad2')],col=COL)

COL=rep('grey',26)
COL[c(19:22)]='red'
par(las=2)
par(mar=c(10,3,3,3))
barplot(t(RPKM.COMBAT)[,which(rownames(RPKM.COMBAT)=='Smad4')],col=COL)

COL=rep('grey',26)
COL[c(23:26)]='red'
par(las=2)
par(mar=c(10,3,3,3))
barplot(t(RPKM.COMBAT)[,which(rownames(RPKM.COMBAT)=='Foxh1')],col=COL)




pbmc <- CreateSeuratObject(counts = RPKM.COMBAT, project = "rna", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=25)


sample_name=colnames(RPKM)

sample_tag=c('D0','D0','D1','D1','D2','D2','D3','D3','D4','D4','D0','D4','D0','D4','D0_S2.3DKO','D4_S2.3DKO','D0_S2.3DKO','D4_S2.3DKO','D0_S4KO',
'D4_S4KO','D0_S4KO','D4_S4KO','D0_Foxh1KO','D4_Foxh1KO','D0_Foxh1KO','D4_Foxh1KO')


pbmc@meta.data$sname=sample_name
pbmc@meta.data$stag=sample_tag

pbmc@meta.data$batch=BATCH
DimPlot(pbmc, reduction = "pca",group.by='batch',label=F,pt.size=3)+NoLegend()


DimPlot(pbmc, reduction = "pca",group.by='stag',label=T,pt.size=3)+NoLegend()

saveRDS(pbmc,'RNA_seurat.rds')





/home/toolkit/tools/R4.0.3/bin/R

#############################################################################
#Heatmap & Cluster

RPKM=read.table('RNA_RPKM_MAT.txt',sep='\t',row.names=1,header=T)
MAT=RPKM[,1:10]
ROW_MEAN_MAT=apply(MAT,1,mean)
MAT=MAT[which(ROW_MEAN_MAT>1),]

SMAT=t(apply(MAT,1,scale))
rownames(SMAT)=rownames(MAT)
colnames(SMAT)=c('D0a','D0b','D1a','D1b','D2a','D2b','D3a','D3b','D4a','D4b')

DIFF=cbind(abs(SMAT[,1]-SMAT[,2]),abs(SMAT[,3]-SMAT[,4]),abs(SMAT[,5]-SMAT[,6]),abs(SMAT[,7]-SMAT[,8]),abs(SMAT[,9]-SMAT[,10]))

MEAN.DIFF=apply(DIFF,1,mean)
SMAT=SMAT[which(MEAN.DIFF<1),]




D0=apply(SMAT[,1:2],1,mean)
D1=apply(SMAT[,3:4],1,mean)
D2=apply(SMAT[,5:6],1,mean)
D3=apply(SMAT[,7:8],1,mean)
D4=apply(SMAT[,9:10],1,mean)

MEAN.MAT=cbind(D0,D1,D2,D3,D4)


D=dist(MEAN.MAT)
H=hclust(D)

#

ha = rowAnnotation(foo = anno_mark(at = c( which(rownames(SMAT)=='Nanog'),
                                           which(rownames(SMAT)=='T'),
                                           which(rownames(SMAT)=='Eomes'),
                                           which(rownames(SMAT)=='Gsc'),
                                           which(rownames(SMAT)=='Mixl1'),
                                           which(rownames(SMAT)=='Foxa2'),
                                           which(rownames(SMAT)=='Foxh1'),
                                           which(rownames(SMAT)=='Pou5f1'),
                                           which(rownames(SMAT)=='Smad2'),
                                           which(rownames(SMAT)=='Smad3'),
                                           which(rownames(SMAT)=='Smad4')
                                            ), 
                         labels = c('Nanog','T','Eomes','Gsc','Mixl1','Foxa2','Foxh1','Pou5f1','Smad2','Smad3','Smad4')))




library(ComplexHeatmap)
library(circlize)
color_fun =colorRamp2(c(-1, 0, 1), c('royalblue3','white','indianred3'))


C=cutree(H,k=8)

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
while(i<=8){
    this.mean=apply(MEAN.MAT[which(C==i),],2,mean)
    this.upper=apply(MEAN.MAT[which(C==i),],2,quantile, 0.975)
    this.lower=apply(MEAN.MAT[which(C==i),],2,quantile, 0.025)
    MEAN=cbind(MEAN,this.mean)
    UPPER=cbind(UPPER,this.upper)
    LOWER=cbind(LOWER,this.lower)
    i=i+1
    }



par(mfrow=c(8,1))
par(mar=c(0.5,1,0.5,1))
i=1
while(i<=8){
    plot(MEAN[,i],ylim=c(-2,2),type='b',pch=16,col='black',axes=FALSE, frame.plot=TRUE,cex=2)
    #arrows(x0=c(1:5),y0=UPPER[,i],x1=c(1:5),y1=LOWER[,i],angle = 90,code = 3,length = 0.1)
    i=i+1
    }



OUT=cbind(C[order(C)], SMAT[order(C),])
colnames(OUT)[1]='Cluster'

OUT=cbind(row.names(OUT),OUT)
colnames(OUT)[1]='GENE'


write.table(OUT,'RNA_GENE_CLUSTER_MAT.txt',row.names=F,col.names=T,sep='\t',quote=F)





################################################










/home/toolkit/tools/R4.0.3/bin/R


setwd('/home/database/data/QW_LAB_data')


library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()



DATA=read.table('RNA_GENE_CLUSTER_MAT_StepMiner.pcl',sep='\t',row.names=1,header=T)
DATA=DATA[,3:ncol(DATA)]


SMAT=as.matrix(DATA)


library(ComplexHeatmap)
library(circlize)
color_fun =colorRamp2(c(-1.2,0,1.2), c('royalblue3','white','indianred3'))





ha = rowAnnotation(foo = anno_mark(at = c( which(rownames(SMAT)=='Nanog'),
                                           which(rownames(SMAT)=='T'),
                                           which(rownames(SMAT)=='Eomes'),
                                           which(rownames(SMAT)=='Gsc'),
                                           which(rownames(SMAT)=='Mixl1'),
                                           which(rownames(SMAT)=='Foxa2'),
                                           which(rownames(SMAT)=='Foxh1'),
                                           which(rownames(SMAT)=='Pou5f1'),
                                           which(rownames(SMAT)=='Smad2'),
                                           which(rownames(SMAT)=='Smad3'),
                                           which(rownames(SMAT)=='Smad4')
                                            ), 
                         labels = c('Nanog','T','Eomes','Gsc','Mixl1','Foxa2','Foxh1','Pou5f1','Smad2','Smad3','Smad4')))




Heatmap(SMAT,row_title='',name="Z",
        cluster_columns=F, cluster_rows=F,
        show_column_dend = F, show_row_dend = F, right_annotation = ha,
        show_column_names=T, show_row_names=F,
        col=color_fun, border = TRUE, 
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        heatmap_width = unit(0.7, "npc")
        )




D3C=read.table('RNA_D3C.txt',header=F,row.names=NULL,sep='\t')[,1]
D4C=read.table('RNA_D4C.txt',header=F,row.names=NULL,sep='\t')[,1]



D3C_EXP=apply(SMAT[which(rownames(SMAT) %in% D3C),],2,mean)
D4C_EXP=apply(SMAT[which(rownames(SMAT) %in% D4C),],2,mean)


plot(D4C_EXP,lwd=2,type='b',col='black',pch=16)
points(D3C_EXP,lwd=2,type='b',col='royalblue3',pch=16)



par(las=2)
X=barplot(t(cbind(D3C_EXP,D4C_EXP))+1.5,beside=T,width=1)




























######################
CCC=apply(SMAT, 1, cor, c(0,0,1,1,2,2,3,3,4,4),method='spearman')
VEC=cbind(1:length(CCC),CCC)


library(gatepoints)

plot(VEC)
selectedPoints <- fhs(VEC, pch=16,col='red3',cex=1,mark = TRUE)
max(which(names(CCC) %in% selectedPoints))




plot(CCC)
abline(v=2081,col='red')
abline(v=3300,col='red')
abline(v=4192,col='red')
abline(v=5240,col='red')
abline(v=5992,col='red')
abline(v=6542,col='red')
abline(v=7698,col='red')
abline(v=8589,col='red')
abline(v=9444,col='red')



SMC=rep(10,length(CCC))
SMC[1:9444]=9
SMC[1:8589]=8
SMC[1:7698]=7
SMC[1:6542]=6
SMC[1:5992]=5
SMC[1:5240]=4
SMC[1:4192]=3
SMC[1:3300]=2
SMC[1:2081]=1




ha = rowAnnotation(foo = anno_mark(at = c( which(rownames(SMAT)=='Nanog'),
                                           which(rownames(SMAT)=='T'),
                                           which(rownames(SMAT)=='Eomes'),
                                           which(rownames(SMAT)=='Gsc'),
                                           which(rownames(SMAT)=='Mixl1'),
                                           which(rownames(SMAT)=='Foxa2'),
                                           which(rownames(SMAT)=='Foxh1'),
                                           which(rownames(SMAT)=='Pou5f1'),
                                           which(rownames(SMAT)=='Smad2'),
                                           which(rownames(SMAT)=='Smad3'),
                                           which(rownames(SMAT)=='Smad4')
                                            ), 
                         labels = c('Nanog','T','Eomes','Gsc','Mixl1','Foxa2','Foxh1','Pou5f1','Smad2','Smad3','Smad4')))




library(ComplexHeatmap)
library(circlize)
color_fun =colorRamp2(c(-1.2, 0, 1.2), c('royalblue3','white','indianred3'))


Heatmap(SMAT,row_title='',name="Z",
        cluster_columns=F, cluster_rows=F,
        show_column_dend = F, show_row_dend = F, right_annotation = ha,
        show_column_names=T, show_row_names=F,
        col=color_fun, border = TRUE, 
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        heatmap_width = unit(0.7, "npc"),
        row_split = SMC
        )


















MEAN=c()
UPPER=c()
LOWER=c()


i=1
while(i<=8){
    this.mean=apply(MEAN.MAT[which(C==i),],2,mean)
    this.upper=apply(MEAN.MAT[which(C==i),],2,quantile, 0.975)
    this.lower=apply(MEAN.MAT[which(C==i),],2,quantile, 0.025)
    MEAN=cbind(MEAN,this.mean)
    UPPER=cbind(UPPER,this.upper)
    LOWER=cbind(LOWER,this.lower)
    i=i+1
    }



par(mfrow=c(8,1))
par(mar=c(0.5,1,0.5,1))
i=1
while(i<=8){
    plot(MEAN[,i],ylim=c(-2,2),type='b',pch=16,col='black',axes=FALSE, frame.plot=TRUE,cex=2)
    #arrows(x0=c(1:5),y0=UPPER[,i],x1=c(1:5),y1=LOWER[,i],angle = 90,code = 3,length = 0.1)
    i=i+1
    }



OUT=cbind(C[order(C)], SMAT[order(C),])
colnames(OUT)[1]='Cluster'

OUT=cbind(row.names(OUT),OUT)
colnames(OUT)[1]='GENE'


write.table(OUT,'RNA_GENE_CLUSTER_MAT.txt',row.names=F,col.names=T,sep='\t',quote=F)


















































































library(dendextend)
row_dend = as.dendrogram(H)

Heatmap(SMAT,row_title='',name="Z",
        cluster_columns=F, cluster_rows=row_dend,
	      show_column_dend = F, show_row_dend = T, 
	      show_column_names=T, show_row_names=F,
	      col=color_fun, border = TRUE, 
	      column_names_gp = gpar(fontsize = 8),
	      row_names_gp = gpar(fontsize = 8),
	      heatmap_width = unit(0.7, "npc"),
	      row_split = 8
	      )



CSMAT=t(.generate_mean(t(SMAT),C))

Heatmap(CSMAT[order(as.numeric(rownames(CSMAT))),],row_title='',name="C",
        cluster_columns=F, cluster_rows=F,
	      show_column_dend = F, show_row_dend = F, 
	      show_column_names=T, show_row_names=F,
	      col=color_fun, border = TRUE, 
	      column_names_gp = gpar(fontsize = 8),
	      row_names_gp = gpar(fontsize = 8),
	      heatmap_width = unit(0.7, "npc")
	      )





















































CSMAT=t(.generate_mean(t(SMAT),C))

COR.SCORE=apply(CSMAT,1,cor,c(0,0,1,1,2,2,3,3,4,4),method='spearman')
OC=order(COR.SCORE)

NEW.C=C
i=1
while(i<=length(OC)){
	NEW.C[which(C==OC[i])]=i
	i=i+1
}


C=NEW.C

























TREND=matrix(0,nrow=nrow(SMAT),ncol=4)
rownames(TREND)=rownames(SMAT)
colnames(TREND)=c('D0_D1','D1_D2','D2_D3','D3_D4')

CUT=0.5

TREND[,1][which(D1-D0> CUT)]=1
TREND[,1][which(D1-D0< -CUT)]=-1
TREND[,2][which(D2-D1> CUT)]=1
TREND[,2][which(D2-D1< -CUT)]=-1
TREND[,3][which(D3-D2> CUT)]=1
TREND[,3][which(D3-D2< -CUT)]=-1
TREND[,4][which(D4-D3> CUT)]=1
TREND[,4][which(D4-D3< -CUT)]=-1



D=dist(TREND)
