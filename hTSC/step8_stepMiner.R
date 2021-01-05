setwd('/home/zhangfeng/project/hTSC/stepminer')


OMAT=read.table('STEP_MINER_EXP_ZVALUE.txt',sep='\t',header=T, row.names=1)

MAT=OMAT[,3:ncol(OMAT)]
MAT=apply(MAT,2,as.numeric)

DIFF=MAT[,16]-MAT[,1]

G1=c(1:1188)
G2=c(1189:2508)


MAT1=MAT[G1,]
MAT2=MAT[G2,]

BMAT1=MAT1
BMAT2=MAT2

BMAT1[which(MAT1>0)]=1
BMAT1[which(MAT1<0)]=0

BMAT2[which(MAT2>0)]=1
BMAT2[which(MAT2<0)]=0


PBMAT1=apply(BMAT1,2,sum)/nrow(BMAT1)*100
MPBMAT1=(PBMAT1[c(1:8)*2]+PBMAT1[c(1:8)*2-1]) /2
barplot(MPBMAT1,ylim=c(0,100))

PBMAT2=apply(BMAT2,2,sum)/nrow(BMAT2)*100
MPBMAT2=(PBMAT2[c(1:8)*2]+PBMAT2[c(1:8)*2-1]) /2
barplot(MPBMAT2,ylim=c(0,100))





library('ComplexHeatmap')
library('circlize')
library('seriation')

mat=MAT2
color_fun =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))

Heatmap(mat,row_title='',name="Exp",
        cluster_column_slices = FALSE, cluster_row_slices = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE
        
        )

