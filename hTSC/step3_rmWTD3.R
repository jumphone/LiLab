
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
MAT=MAT[,2:(ncol(MAT)-1)]
mat.colname=c('V.12h_r1','V.18h_r1','V.3h_r1',
      'V.6h_r1','V.D1_r1','V.D2_r1','V.D3_r1',
      'V.D1_r2','V.D2_r2','V.D3_r2','V.12h_r2',
       'V.18h_r2','V.3h_r2','V.6h_r2','WT.3h',
       'WT.6h')
colnames(MAT)=mat.colname
TAG=c('V.12h','V.18h','V.3h',
      'V.6h','V.D1','V.D2','V.D3',
      'V.D1','V.D2','V.D3','V.12h',
       'V.18h','V.3h','V.6h','WT.3h',
       'WT.6h')
       
PCGENE=read.table('/home/database/annotation/hg19/Homo_sapiens.GRCh37.75.chr.pc.gene.SYM.bed',sep='\t',header=F)
MAT=MAT[which(rownames(MAT) %in% PCGENE[,4]),]

pbmc <- CreateSeuratObject(counts = MAT,  project = "hTSC", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features =all.genes)
pbmc <- RunPCA(pbmc, features = all.genes, npcs = 10)
DimPlot(pbmc, reduction = "pca", pt.size=3) #+ NoLegend()

DimPlot(pbmc, reduction = "pca", pt.size=3 ,label=TRUE)+ NoLegend()
#################################


DATA=as.matrix(pbmc@assays$RNA@data)
VAR=apply(DATA,1,var)
SORT_VAR=sort(VAR)
used_var=which(SORT_VAR >= abs(sort(-VAR)[1000]) )

plot(x=1:length(SORT_VAR),y=SORT_VAR, pch=16,col='grey70', xlab='Index', ylab='VAR')
points(x=(1:length(SORT_VAR))[used_var], y=SORT_VAR[used_var] , pch=16,col='black')



################################


S.DATA=as.matrix(pbmc@assays$RNA@scale.data)
mat=S.DATA[which(rownames(S.DATA) %in% names(used_var)),]
'''
"4.1" "5.1" "2.1"  "3.1"  "6.1"  "7.1"
"8.1"  "6.2"  "7.2"  "8.2"  "4.2" "5.2"
"2.2"  "3.2"  "1.1"    "1.2"

'''


TO=c("4.1","5.1", "2.1" , "3.1"  ,"6.1" , "7.1",
"8.1" , "6.2" , "7.2" , "8.2" , "4.2", "5.2",
"2.2" , "3.2" , "1.1"  ,  "1.2")

TO=as.numeric(TO)

mat=mat[,order(TO)]


library('ComplexHeatmap')
library('circlize')
library('seriation')



color_fun =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))

set.seed(12345)
K=6
C=kmeans(mat, K)$cluster
NEW.C=C
NEW.C[which(C==4)]=1
NEW.C[which(C==2)]=2
NEW.C[which(C==5)]=3
NEW.C[which(C==1)]=4
NEW.C[which(C==3)]=5
NEW.C[which(C==6)]=6

split <- paste0("C", NEW.C)

o.mat=mat[order(GENE_CLUSTER[,2]),]
o.split=split[order(split)]
Heatmap(o.mat,row_title='',name="Exp",
        split=o.split, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K)))
        )
        
        

pdf('./RESULT/Heatmap.pdf',width=5,height=4)

HM <- Heatmap(o.mat,row_title='',name="Exp",
        split=o.split, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K)))
        )
       
HM


dev.off()   


OUT=cbind(o.split,o.mat)
colnames(OUT)[1]='GENE\tCLUSTER'
write.table(OUT,file='./RESULT/OMAT.txt',sep='\t',row.names=T,col.names=T,quote=F)

write.table(OUT[,1],file='./RESULT/GENE_CLUSTER.txt',sep='\t',row.names=T,col.names=T,quote=F)



C1=apply(o.mat[which(o.split=='C1'),],2,mean)
C2=apply(o.mat[which(o.split=='C2'),],2,mean)
C3=apply(o.mat[which(o.split=='C3'),],2,mean)      
C4=apply(o.mat[which(o.split=='C4'),],2,mean)    
C5=apply(o.mat[which(o.split=='C5'),],2,mean)    
C6=apply(o.mat[which(o.split=='C6'),],2,mean)    
         
ha = HeatmapAnnotation(C1 = anno_lines(C1, smooth = TRUE, add_points = FALSE, gp = gpar(lwd=2)), gap = unit(1, "mm") ,
                       C2 = anno_lines(C2, smooth = TRUE, add_points = FALSE, gp = gpar(lwd=2)) ,
                       C3 = anno_lines(C3, smooth = TRUE, add_points = FALSE, gp = gpar(lwd=2)) ,
                       C4 = anno_lines(C4, smooth = TRUE, add_points = FALSE, gp = gpar(lwd=2)) ,
                       C5 = anno_lines(C5, smooth = TRUE, add_points = FALSE, gp = gpar(lwd=2)) ,
                       C6 = anno_lines(C6, smooth = TRUE, add_points = FALSE, gp = gpar(lwd=2))
                       )        
       
Heatmap(o.mat,row_title='',name="Exp",
        split=o.split, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K))),
        bottom_annotation=ha
        )       
       

pdf('./RESULT/Heatmap_trend.pdf',width=5,height=4)

ha = HeatmapAnnotation( gap = unit(1, "mm") ,
                       C1 = anno_lines(C1, smooth = TRUE, add_points = TRUE, gp = gpar(lwd=2), ylim = c(-2.5, 2.5)), 
                       C2 = anno_lines(C2, smooth = TRUE, add_points = TRUE, gp = gpar(lwd=2), ylim = c(-2.5, 2.5)) ,
                       C3 = anno_lines(C3, smooth = TRUE, add_points = TRUE, gp = gpar(lwd=2), ylim = c(-2.5, 2.5)) ,
                       C4 = anno_lines(C4, smooth = TRUE, add_points = TRUE, gp = gpar(lwd=2), ylim = c(-2.5, 2.5)) ,
                       C5 = anno_lines(C5, smooth = TRUE, add_points = TRUE, gp = gpar(lwd=2), ylim = c(-2.5, 2.5)) ,
                       C6 = anno_lines(C6, smooth = TRUE, add_points = TRUE, gp = gpar(lwd=2), ylim = c(-2.5, 2.5))
                       )       
                       
HM=Heatmap(o.mat,row_title='',name="Exp",
        #split=o.split, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        #left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K))),
        bottom_annotation=ha
        ) 
        
draw(HM, heatmap_legend_side = "left", annotation_legend_side = "bottom")

dev.off()          
       
        
        
       
       
       
       
       
        
        
