
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
# All samples
pbmc <- CreateSeuratObject(counts = MAT,  project = "hTSC", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features =all.genes)
pbmc <- RunPCA(pbmc, features = all.genes, npcs = 10)
DimPlot(pbmc, reduction = "pca", pt.size=3) #+ NoLegend()



#################################
# Only V1, V7, V10, WT

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()
source('/home/zhangfeng/project/EVOL/source/MultiTools.R')
setwd('/home/zhangfeng/project/hTSC/data/hTSC_paper_RNA_seq/htseq_count')

###
#Load data
MAT=read.table('MAT_FULL.NAME.txt',header=T,row.names=1,sep='\t')
MAT[1:5,1:5]
TYPE=readRDS('TYPE.rds')
MAT=MAT[,which(TYPE=='NA')]
PCGENE=read.table('/home/database/annotation/hg19/Homo_sapiens.GRCh37.75.chr.pc.gene.SYM.bed',sep='\t',header=F)
MAT=MAT[which(rownames(MAT) %in% PCGENE[,4]),]

pbmc <- CreateSeuratObject(counts = MAT,  project = "hTSC", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features =all.genes)
pbmc <- RunPCA(pbmc, features = all.genes, npcs = 10)
DimPlot(pbmc, reduction = "pca", pt.size=3) #+ NoLegend()

VEC=pbmc@reductions$pca@cell.embeddings
TEXT=str_replace(colnames(pbmc),pattern='NA_','')

plot(VEC, pch=16, cex=2, col='grey70',xlim=c(-80,170),ylim=c(-100,80))
set.seed(8)
text(TEXT, x=VEC[,1],y=VEC[,2], pos=sample(c(1,2,3,4,2,3,1,4),nrow(VEC),replace=T))
#text(TEXT, x=VEC[,1],y=VEC[,2], pos=c(1,2,3,4))

saveRDS(pbmc, file='pbmc_onlyVWT.rds')

DATA=as.matrix(pbmc@assays$RNA@data)
################
COR1=apply(DATA,1, cor, VEC[,1],method='spearman')
COR1[which(is.na(COR1))]=0
TOP.COR1=COR1[order(-COR1)][1:100]
write.table(TOP.COR1,file='TOP.COR1.txt',sep='\t',quote=F,col.names=F,row.names=T)

################
COR2=apply(DATA,1, cor, VEC[,2],method='spearman')
COR2[which(is.na(COR2))]=0
TOP.COR2=COR2[order(-COR2)][1:100]
write.table(TOP.COR2,file='TOP.COR2.txt',sep='\t',quote=F,col.names=F,row.names=T)



colnames(pbmc)
'''
 [1] "NA_WT.ES8P74.D3" "NA_V10.12h"      "NA_V10.18h"      "NA_V10.3h"
 [5] "NA_V10.6h"       "NA_V10.D1"       "NA_V10.D2"       "NA_V10.D3"
 [9] "NA_V1.D1"        "NA_V1.D2"        "NA_V1.D3"        "NA_V7.12h"
[13] "NA_V7.18h"       "NA_V7.3h"        "NA_V7.6h"        "NA_WT.3h"
[17] "NA_WT.6h"        "NA_WT.DMSO.D3"

'''
TAG=c('WT.D3','V.12h','V.18h','V.3h',
      'V.6h','V.D1','V.D2','V.D3',
      'V.D1','V.D2','V.D3','V.12h',
       'V.18h','V.3h','V.6h','WT.3h',
       'WT.6h','WT.D3')

Idents(pbmc)=TAG
DimPlot(pbmc, reduction = "pca", pt.size=3, label=T) + NoLegend()

DATA=as.matrix(pbmc@assays$RNA@data)

VAR=apply(DATA,1,var)

SORT_VAR=sort(VAR)


plot(x=1:length(SORT_VAR),y=SORT_VAR, pch=16,col='grey70', xlab='Index', ylab='VAR')
used_var=which(SORT_VAR >quantile(SORT_VAR, 0.95))
points(x=(1:length(SORT_VAR))[used_var], y=SORT_VAR[used_var] , pch=16,col='black')


##############################
S.DATA=as.matrix(pbmc@assays$RNA@scale.data)
mat=S.DATA[which(rownames(S.DATA) %in% names(used_var)),]
mat.colname=c('WT.D3','V.12h_r1','V.18h_r1','V.3h_r1',
      'V.6h_r1','V.D1_r1','V.D2_r1','V.D3_r1',
      'V.D1_r2','V.D2_r2','V.D3_r2','V.12h_r2',
       'V.18h_r2','V.3h_r2','V.6h_r2','WT.3h',
       'WT.6h','WT.D3')

colnames(mat)=mat.colname


library('ComplexHeatmap')
library('circlize')
library('seriation')

color_fun =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))
Heatmap(mat,row_title='',name="C",
        cluster_columns=TRUE, cluster_rows=TRUE,
	      show_column_dend = TRUE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE
        )
        
pdf('Heatmap.pdf',width=5,height=4)

Heatmap(mat,row_title='',name="C",
        cluster_columns=TRUE, cluster_rows=TRUE,
	      show_column_dend = TRUE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE
        )

dev.off()


#D=dist(mat)
#H=hclust(D)
#C=cutree(H,k=5)




HM=Heatmap(mat,row_title='',name="C",
        cluster_columns=TRUE, cluster_rows=TRUE,
	      show_column_dend = TRUE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:5))),
        row_km = 5
        )


pdf('Heatmap_anno.pdf',width=5,height=4)

HM=Heatmap(mat,row_title='',name="C",
        cluster_columns=TRUE, cluster_rows=TRUE,
	      show_column_dend = TRUE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:5))),
        row_km = 5
        )
HM  
dev.off()


ht = draw(HM)
tmp=row_dend(ht)
names(tmp)



GENE.OUT=c(labels(tmp$'1'),
    labels(tmp$'2'),
    labels(tmp$'3'),
    labels(tmp$'4'),
    labels(tmp$'5')
	   )
GENE.OUT.TYPE=c( rep('1',length(labels(tmp$'1'))),
                 rep('2',length(labels(tmp$'2'))),
                 rep('3',length(labels(tmp$'3'))),
                 rep('4',length(labels(tmp$'4'))),
                 rep('5',length(labels(tmp$'5')))
                     )
OUT=cbind(GENE.OUT,GENE.OUT.TYPE)
write.table(OUT,'GENE_5CLUSTER.txt',sep='\t',row.names=F,col.names=F,quote=F)
