#setwd('C:/Users/cchmc/Desktop/Lingjie')
#MAT=read.table('DATA.txt',header=T,row.names=1,sep='\t')

setwd('/home/zhangfeng/project/CD73')
MAT=read.table('OUTPUT_COMPARE.txt',header=T,row.names=1,sep='\t')



DIFF_TCMTEM=read.table('OUTPUT_DIFF_TCMTEM.txt',header=TRUE,sep='\t',row.names=NULL)
DIFF_CD73=read.table('INPUT_DIFF_CD73.txt',header=TRUE,sep='\t',row.names=NULL)


USED_GENE1=names(which(table(DIFF_TCMTEM[,1])==1))
USED_GENE2=names(which(table(DIFF_CD73[,1])==1))

DIFF_TCMTEM=DIFF_TCMTEM[which(DIFF_TCMTEM[,1] %in% USED_GENE1),]
DIFF_CD73=DIFF_CD73[which(DIFF_CD73[,1] %in% USED_GENE2),]

rownames(DIFF_TCMTEM)=DIFF_TCMTEM[,1]
rownames(DIFF_CD73)=DIFF_CD73[,1]

DIFF_TCMTEM=DIFF_TCMTEM[,2:ncol(DIFF_TCMTEM)]
DIFF_CD73=DIFF_CD73[,2:ncol(DIFF_CD73)]

colnames(DIFF_TCMTEM)=paste0('TCMTEM_',colnames(DIFF_TCMTEM))
colnames(DIFF_CD73)=paste0('CD73_',colnames(DIFF_CD73))

COM=.simple_combine(DIFF_TCMTEM,DIFF_CD73)$combine
COM=COM[22:nrow(COM),]


log2_TEM_div_TCM = -COM$TCMTEM_log2FoldChange 
log2_POS_div_NEG = COM$CD73_log2FoldChange


FCUT=1.5
PCUT=0.05


TEM=rownames(COM)[which( log2_TEM_div_TCM > log(FCUT,2)   & COM$TCMTEM_pvalue < PCUT)]
TCM=rownames(MAT)[which( log2_TEM_div_TCM < -log(FCUT,2)  & COM$TCMTEM_pvalue < PCUT)]

CD73POS=rownames(MAT)[which( log2_POS_div_NEG > log(FCUT,2)  & COM$CD73_pvalue < PCUT )]
CD73NEG=rownames(MAT)[which( log2_POS_div_NEG < -log(FCUT,2) & COM$CD73_pvalue < PCUT )]





#install.packages("VennDiagram")
library(VennDiagram)

x=list(TEM=TEM, TCM=TCM, CD73POS=CD73POS, CD73NEG=CD73NEG)
venn.diagram(x, filename='VENN_FOLD_1.5_PV_0.05.tiff', height = 3000, width = 3000, resolution =
    500, imagetype = "tiff", units = "px", compression =
    "lzw", na = "stop", main = NULL, sub = NULL, main.pos
    = c(0.5, 1.05), main.fontface = "plain",
    main.fontfamily = "serif", main.col = "black",
    main.cex = 1, main.just = c(0.5, 1), sub.pos = c(0.5,
    1.05), sub.fontface = "plain", sub.fontfamily =
    "serif", sub.col = "black", sub.cex = 1, sub.just =
    c(0.5, 1), category.names = names(x), force.unique =
    TRUE, print.mode = "raw", sigdigs = 3, direct.area =
    FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, 
    lower.tail = TRUE)





TEM_and_CD73POS=TEM[which(TEM %in% CD73POS)]
TEM_and_CD73NEG=TEM[which(TEM %in% CD73NEG)]

TCM_and_CD73POS=TCM[which(TCM %in% CD73POS)]
TCM_and_CD73NEG=TCM[which(TCM %in% CD73NEG)]


write.table(TEM_and_CD73POS,file='OVERLAP_TEM_and_CD73POS.txt',col.names=F,row.names=F,sep='\t',quote=F)
write.table(TEM_and_CD73NEG,file='OVERLAP_TEM_and_CD73NEG.txt',col.names=F,row.names=F,sep='\t',quote=F)
write.table(TCM_and_CD73POS,file='OVERLAP_TCM_and_CD73POS.txt',col.names=F,row.names=F,sep='\t',quote=F)
write.table(TCM_and_CD73NEG,file='OVERLAP_TCM_and_CD73NEG.txt',col.names=F,row.names=F,sep='\t',quote=F)


X=log2_TEM_div_TCM
Y=log2_POS_div_NEG
X[is.na(X)]=0
Y[is.na(Y)]=0


VCUT=8
X[which(X>VCUT)]=VCUT
X[which(X< -VCUT)]= -VCUT
Y[which(Y>VCUT)]=VCUT
Y[which(Y< -VCUT)]= -VCUT

plot(X,Y,pch=16,col='grey50',xlab='log2( TEM / TCM )',ylab='log2( CD73+ / CD73- )')

pdf('TEM_TCM_CD73pos_neg_plot.pdf',width=7,height=7)
plot(X,Y,pch=16,col='grey50',xlab='log2( TEM / TCM )',ylab='log2( CD73+ / CD73- )',main=paste0('COR=',cor(X,Y)))
dev.off()










































TEM_index=which( log2_TEM_div_TCM > log(FCUT,2)   & COM$TCMTEM_pvalue < PCUT)
TCM_index=which( log2_TEM_div_TCM < -log(FCUT,2)  & COM$TCMTEM_pvalue < PCUT)

CD73POS_index=which( log2_POS_div_NEG > log(FCUT,2)  & COM$CD73_pvalue < PCUT )
CD73NEG_index=which( log2_POS_div_NEG < -log(FCUT,2) & COM$CD73_pvalue < PCUT )


SIG_index=which(COM$TCMTEM_pvalue < PCUT | COM$CD73_pvalue < PCUT)


plot(X[SIG_index],Y[SIG_index] ,pch=16,col='grey90')
points(X[TEM_index],Y[TEM_index],pch=16,col='red')
points(X[TCM_index],Y[TCM_index],pch=16,col='red')
points(X[CD73POS_index],Y[CD73POS_index],pch=16,col='red')
points(X[CD73NEG_index],Y[CD73NEG_index],pch=16,col='red')













































FCUT=1.5
PCUT=0.05

#TEM=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PADJ< PCUT )]
#TCM=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PADJ< PCUT )]

#CD73POS=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PADJ < PCUT )]
#CD73NEG=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PADJ < PCUT )]

TEM=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PVALUE< PCUT )]
TCM=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PVALUE< PCUT )]

CD73POS=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PVALUE < PCUT )]
CD73NEG=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PVALUE < PCUT )]


#install.packages("VennDiagram")
library(VennDiagram)

x=list(TEM=TEM, TCM=TCM, CD73POS=CD73POS, CD73NEG=CD73NEG)
venn.diagram(x, filename='VENN_FOLD_1.5_PV_0.05.tiff', height = 3000, width = 3000, resolution =
    500, imagetype = "tiff", units = "px", compression =
    "lzw", na = "stop", main = NULL, sub = NULL, main.pos
    = c(0.5, 1.05), main.fontface = "plain",
    main.fontfamily = "serif", main.col = "black",
    main.cex = 1, main.just = c(0.5, 1), sub.pos = c(0.5,
    1.05), sub.fontface = "plain", sub.fontfamily =
    "serif", sub.col = "black", sub.cex = 1, sub.just =
    c(0.5, 1), category.names = names(x), force.unique =
    TRUE, print.mode = "raw", sigdigs = 3, direct.area =
    FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, 
    lower.tail = TRUE)



X=MAT$TEMdivTCM
Y=MAT$POSdivNEG

plot(log(X,FCUT),log(Y,FCUT))


plot(X,Y)
























































FCUT=2
PCUT=0.05

TEM=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PADJ< PCUT )]
TCM=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PADJ< PCUT )]

CD73POS=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PADJ < PCUT )]
CD73NEG=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PADJ < PCUT )]


#install.packages("VennDiagram")
library(VennDiagram)

x=list(TEM=TEM, TCM=TCM, CD73POS=CD73POS, CD73NEG=CD73NEG)
venn.diagram(x, filename='VENN_FOLD_2.tiff', height = 3000, width = 3000, resolution =
    500, imagetype = "tiff", units = "px", compression =
    "lzw", na = "stop", main = NULL, sub = NULL, main.pos
    = c(0.5, 1.05), main.fontface = "plain",
    main.fontfamily = "serif", main.col = "black",
    main.cex = 1, main.just = c(0.5, 1), sub.pos = c(0.5,
    1.05), sub.fontface = "plain", sub.fontfamily =
    "serif", sub.col = "black", sub.cex = 1, sub.just =
    c(0.5, 1), category.names = names(x), force.unique =
    TRUE, print.mode = "raw", sigdigs = 3, direct.area =
    FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, 
    lower.tail = TRUE)








FCUT=1.2
PCUT=0.05

TEM=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PADJ< PCUT )]
TCM=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PADJ< PCUT )]

CD73POS=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PADJ < PCUT )]
CD73NEG=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PADJ < PCUT )]


#install.packages("VennDiagram")
library(VennDiagram)

x=list(TEM=TEM, TCM=TCM, CD73POS=CD73POS, CD73NEG=CD73NEG)
venn.diagram(x, filename='VENN_FOLD_1.2.tiff', height = 3000, width = 3000, resolution =
    500, imagetype = "tiff", units = "px", compression =
    "lzw", na = "stop", main = NULL, sub = NULL, main.pos
    = c(0.5, 1.05), main.fontface = "plain",
    main.fontfamily = "serif", main.col = "black",
    main.cex = 1, main.just = c(0.5, 1), sub.pos = c(0.5,
    1.05), sub.fontface = "plain", sub.fontfamily =
    "serif", sub.col = "black", sub.cex = 1, sub.just =
    c(0.5, 1), category.names = names(x), force.unique =
    TRUE, print.mode = "raw", sigdigs = 3, direct.area =
    FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, 
    lower.tail = TRUE)






FCUT=1
PCUT=0.05

TEM=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PADJ< PCUT )]
TCM=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PADJ< PCUT )]

CD73POS=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PADJ < PCUT )]
CD73NEG=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PADJ < PCUT )]


#install.packages("VennDiagram")
library(VennDiagram)

x=list(TEM=TEM, TCM=TCM, CD73POS=CD73POS, CD73NEG=CD73NEG)
venn.diagram(x, filename='VENN_FOLD_1.tiff', height = 3000, width = 3000, resolution =
    500, imagetype = "tiff", units = "px", compression =
    "lzw", na = "stop", main = NULL, sub = NULL, main.pos
    = c(0.5, 1.05), main.fontface = "plain",
    main.fontfamily = "serif", main.col = "black",
    main.cex = 1, main.just = c(0.5, 1), sub.pos = c(0.5,
    1.05), sub.fontface = "plain", sub.fontfamily =
    "serif", sub.col = "black", sub.cex = 1, sub.just =
    c(0.5, 1), category.names = names(x), force.unique =
    TRUE, print.mode = "raw", sigdigs = 3, direct.area =
    FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, 
    lower.tail = TRUE)
