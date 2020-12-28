


setwd('/home/database/data/backup_20201022_new/H9DMSO_Len_lingjie-b67e54/process/DESeq2_fc2')

ORIG_MAT=read.table('DESeq2_result.txt',row.names=1,header=TRUE,sep='\t')
head(ORIG_MAT)

MAT_RPM=as.matrix(ORIG_MAT[,5:8])
head(MAT_RPM)


MAT_RPM_STD=t(apply(MAT_RPM, 1, scale))
MAT_RPM_STD[1:4,1:4]
colnames(MAT_RPM_STD)=colnames(MAT_RPM)
rownames(MAT_RPM_STD)=rownames(MAT_RPM)


library('ComplexHeatmap')
library('circlize')
library('seriation')
mat=MAT_RPM_STD
color_fun =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))

set.seed(12345)
K=2
C=kmeans(mat, K)$cluster
NEW.C=C
split <- paste0("C", NEW.C)

Heatmap(mat,row_title='',name="Exp",split=split,
         cluster_column_slices = FALSE, cluster_row_slices = FALSE,
         cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
       left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K)))
        )

pdf('~/demo/2heatmap.pdf',width=5,height=8)
Heatmap(mat,row_title='',name="Exp",split=split,
         cluster_column_slices = FALSE, cluster_row_slices = FALSE,
         cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE,
       left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:K)))
        )
dev.off()



############################
LOGRPM=log(apply(MAT_RPM,1,mean))
LOG2FC=ORIG_MAT$log2FoldChange


X=LOGRPM 
Y=LOG2FC
plot(X, Y, pch=16,xlab='log2RPM',ylab='log2FC')
select_point1=which(LOG2FC>1)
points(X[select_point1],Y[select_point1],pch=16,col='indianred3')
select_point2=which(LOG2FC< -1)
points(X[select_point2],Y[select_point2],pch=16,col='royalblue3')

pdf('~/demo/3MA.pdf',width=7,height=7)
X=LOGRPM 
Y=LOG2FC
plot(X, Y, pch=16,xlab='log2RPM',ylab='log2FC')
select_point1=which(LOG2FC>1)
points(X[select_point1],Y[select_point1],pch=16,col='indianred3')
select_point2=which(LOG2FC< -1)
points(X[select_point2],Y[select_point2],pch=16,col='royalblue3')
dev.off()





