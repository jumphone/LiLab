

reticulate::use_python('/home/toolkit/local/bin/python3',required=TRUE)

library(dplyr)
library(Seurat)
library(patchwork)

source('/home/toolkit/src/BEER.R')
source('bbknn.R')
library(harmony)

data=readRDS(file='data/mat_data_h9.rds')
batch=readRDS(file='data/mat_batch_h9.rds')

sobj=CreateSeuratObject(counts = data, min.cells = 10, min.features = 0, project = "ALL")
sobj$batch=batch

sobj <- subset(sobj, subset = nFeature_RNA > 200 )

DATA=sobj[['RNA']]@counts
BATCH=sobj$batch

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE )

saveRDS(mybeer, 'data/mybeer.rds')

#########################
sobj=mybeer$seurat
oriPCA=sobj@reductions$pca@cell.embeddings
PCA=oriPCA[,mybeer$select]

HAR=harmony::HarmonyMatrix(
  data_mat  = PCA,       # Matrix with coordinates for each cell (row) along many PCs (columns)
  meta_data = sobj$batch, # Dataframe with information for each cell (row)
  do_pca    = FALSE      # Since we are providing PCs, do not run PCA
  )

BBKNN=.bbknn(HAR, sobj$batch, NB=3, NT=10, DM=2)

UMAP=BBKNN
rownames(UMAP)=colnames(sobj)
colnames(UMAP)=paste0('UMAP_',c(1:ncol(UMAP)))
sobj$batch=sobj$batch
sobj@reductions$umap@cell.embeddings=UMAP

sobj@reductions$harmony=HAR
sobj@reductions$bbknn=BBKNN

saveRDS(sobj, 'data/sobj_chbb.rds')
saveRDS(BBKNN, 'data/sobj_chbb_bbknn.rds')




pdf('plot/p00_umap_batch.pdf',width=6,height=5)
DimPlot(sobj, group.by='batch')
dev.off()



pdf('plot/p01_combined_batch_umap.pdf',width=7,height=7)

VEC=BBKNN
plot(VEC,pch=16,col='grey70',cex=0.5)
points(VEC[which(BATCH=='D7'),],pch=16,col='red1',cex=0.5)

plot(VEC,pch=16,col='grey70',cex=0.5)
points(VEC[which(BATCH=='D21'),],pch=16,col='red1',cex=0.5)

plot(VEC,pch=16,col='grey70',cex=0.5)
points(VEC[which(BATCH=='D41'),],pch=16,col='red1',cex=0.5)

plot(VEC,pch=16,col='grey70',cex=0.5)
points(VEC[which(BATCH=='KC'),],pch=16,col='red1',cex=0.5)

dev.off()




source('/home/toolkit/src/fitdevo.R')
BGW=readRDS('/home/toolkit/src/BGW.rds')
#DP=fitdevo(sobj[['RNA']]@counts,BGW)
#saveRDS(DP,file='data/fitdevo_dp.rds')

DP=readRDS(file='data/fitdevo_dp.rds')

VEC=BBKNN
#FIELD=fitdevo.field(DP, VEC)



pdf('plot/p02_fitdevo.pdf',width=7,height=7)

V=DP
COL=.vcol(V, c(min(V),median(V),max(V)),c('blue','gold1','red'))
plot(BBKNN[order(V),], col=COL[order(V)],pch=16,cex=0.5)

plot(x=V,y=rep(1,length(V)),col=COL,type='h',lwd=2,ylim=c(0,1))

FIELD=fitdevo.field(DP, VEC)

dev.off()





pdf('plot/p03_umap_gene.pdf',width=5,height=5)
source('/home/toolkit/src/VISA.R')

GENES=c('TACSTD2','DLK1','PAX6','KRT18','KRT8','KRT5','KRT14','TP63')

for(GENE in GENES){
    GENE=GENE
    V=sobj[['RNA']]@data[GENE,]
    COL=visa.vcol(V, c(0, median(V[which(V>0)]) ,quantile(V[which(V>0)],0.99)),c('blue','grey95','red'))
    plot(BBKNN[order(V),], col=COL[order(V)], pch=16,cex=0.5,main=GENE)
    }

dev.off()






