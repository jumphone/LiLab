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



