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

TO=c("4.1","5.1", "2.1" , "3.1"  ,"6.1" , "7.1",
"8.1" , "6.2" , "7.2" , "8.2" , "4.2", "5.2",
"2.2" , "3.2" , "1.1"  ,  "1.2")

TO=as.numeric(TO)

MAT=MAT[,order(TO)]
TAG=TAG[order(TO)]
TO=TO[order(TO)]

      
PCGENE=read.table('/home/database/annotation/hg19/Homo_sapiens.GRCh37.75.chr.pc.gene.SYM.bed',sep='\t',header=F)
MAT=MAT[which(rownames(MAT) %in% PCGENE[,4]),]
rownames(PCGENE)=PCGENE[,4]

RSUM=apply(MAT,1,sum)
MAT=MAT[which(RSUM>0),]
#########################################

COM=.simple_combine(MAT,PCGENE)
COUNT=COM$exp_sc_mat1
GINFO=COM$exp_sc_mat2

.normRPKM <- function(x, L){
    x_sum=sum(x)
    #print(x_sum)
    if(x_sum==0){
        y=x
    }else{
        y=x  / L * 1000 / x_sum * 1000000
        }
    return(y)
    }
    
  

RPKM=apply(COUNT, 2, .normRPKM, GINFO[,5])
RPKM.LOG=log(RPKM+1,10)
RPKM.LOG.SCALE=t(apply(RPKM.LOG,1,scale))
colnames(RPKM.LOG.SCALE)=colnames(COUNT)
rownames(RPKM.LOG.SCALE)=rownames(COUNT)
#####################################################



pca=prcomp(t(RPKM.LOG.SCALE),center = FALSE, scale. = FALSE)

max_pc1=max(pca$x[,1])
min_pc1=min(pca$x[,1])
max_pc2=max(pca$x[,2])
min_pc2=min(pca$x[,2])
add_l=20

plot(pca$x[,1],pca$x[,2],xlim=c(min_pc1-add_l,max_pc1+add_l),ylim=c(min_pc2-add_l,max_pc2+add_l),pch=16,
     xlab='PC1',ylab='PC2')
set.seed(123)
text(label=rownames(pca$x),x=pca$x[,1], y=pca$x[,2], pos=sample(c(1,2,3,4),nrow(pca$x),replace = TRUE))



####################################
STAR_GENE=c('POU5F1','NANOG','SOX2','TFAP2C','KRT8','KRT18','LRP2','FLT1','KDR','CDX2','KRT7','VGLL1','GATA2','GATA3','TFAP2A','HAND1','MYC','ZIC2','NR2F2','ZIC5')
STAR_GENE_INDEX=which(rownames(RPKM) %in% STAR_GENE)
######################



########################
MAX_RPKM=apply(RPKM,1,max)
MAX_RPKM[STAR_GENE_INDEX]
##################
IN_REP_DIFF= abs(RPKM.LOG.SCALE[,(1:8)*2]-RPKM.LOG.SCALE[,((1:8)*2-1)])
MAX_IN_REP_DIFF=apply(IN_REP_DIFF,1,max)
MAX_IN_REP_DIFF[STAR_GENE_INDEX]
###################
OUT_REP_RPKM=(RPKM[,(1:8)*2]  +  RPKM[,((1:8)*2-1)]   )/2

LOGFD=c()
i=1
#while(i<=ncol(OUT_REP_RPKM)){
while(i<=1){
   j=i+1
   while(j<=ncol(OUT_REP_RPKM)){
      this_log_fd=log((OUT_REP_RPKM[,i]+0.01)/(OUT_REP_RPKM[,j]+0.01),2)
      LOGFD=cbind(LOGFD, this_log_fd)
      j=j+1
      }
   i=i+1
   }

ABSLOGFD=abs(LOGFD)

MAX_ABSLOGFD=apply(ABSLOGFD,1,max)
MAX_FD=2**MAX_ABSLOGFD
MAX_FD[STAR_GENE_INDEX]


##############################################


used_gene_index=which(MAX_RPKM >1 & MAX_IN_REP_DIFF < 2 & MAX_FD > 2)

length(used_gene_index)

length(which(STAR_GENE %in% names(used_gene_index)))



library('ComplexHeatmap')
library('circlize')
library('seriation')

mat=RPKM.LOG.SCALE[which(rownames(RPKM.LOG.SCALE) %in% names(used_gene_index)),]

color_fun =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))

set.seed(12345)
K=6
C=kmeans(mat, K)$cluster
NEW.C=C

NEW.C[which(C==4)]=1
NEW.C[which(C==1)]=2
NEW.C[which(C==5)]=3
NEW.C[which(C==6)]=4
NEW.C[which(C==3)]=5
NEW.C[which(C==2)]=6


split <- paste0("C", NEW.C)

o.mat=mat[order(NEW.C),]
o.split=split[order(split)]

Heatmap(o.mat,row_title='',name="Exp",
        split=o.split, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K)))
        )

OUT=cbind(o.split,o.mat)
colnames(OUT)[1]='GENE\tCLUSTER'
write.table(OUT,file='./RESULT/Heatmap_zvalue_vsWT.txt',sep='\t',row.names=T,col.names=T,quote=F)


pdf('./RESULT/Heatmap_vsWT.pdf',width=5,height=4)

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

OUT=NEW.C[order(NEW.C)]
write.table(OUT,file='./RESULT/GENE_CLUSTER_vsWT.txt',sep='\t',row.names=T,col.names=F,quote=F)

####################

# DEG

DEG=as.character(read.table('DEG_HuiLi.txt',header=FALSE,row.names=NULL,sep='\t')[,1])
used_gene_index=which(rownames(RPKM) %in% DEG)
names(used_gene_index)=rownames(RPKM)[used_gene_index]
length(which(STAR_GENE %in% names(used_gene_index)))


library('ComplexHeatmap')
library('circlize')
library('seriation')

mat=RPKM.LOG.SCALE[which(rownames(RPKM.LOG.SCALE) %in% names(used_gene_index)),]

color_fun =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))

set.seed(12345)
K=6
C=kmeans(mat, K)$cluster
NEW.C=C

NEW.C[which(C==2)]=1
NEW.C[which(C==4)]=2
NEW.C[which(C==3)]=3
NEW.C[which(C==1)]=4
NEW.C[which(C==6)]=5
NEW.C[which(C==5)]=6


split <- paste0("C", NEW.C)

o.mat=mat[order(NEW.C),]
o.split=split[order(split)]

Heatmap(o.mat,row_title='',name="Exp",
        split=o.split, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K)))
        )

OUT=cbind(o.split,o.mat)
colnames(OUT)[1]='GENE\tCLUSTER'
write.table(OUT,file='./RESULT/Heatmap_zvalue_DEG.txt',sep='\t',row.names=T,col.names=T,quote=F)


pdf('./RESULT/Heatmap_DEG.pdf',width=5,height=4)

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

OUT=NEW.C[order(NEW.C)]
write.table(OUT,file='./RESULT/GENE_CLUSTER_DEG.txt',sep='\t',row.names=T,col.names=F,quote=F)



