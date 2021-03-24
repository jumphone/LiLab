#setwd('C:/Users/cchmc/Desktop/Lingjie')

#MAT=read.table('DATA.txt',header=T,row.names=1,sep='\t')

setwd('/home/zhangfeng/project/CD73')
MAT=read.table('OUTPUT_COMPARE.txt',header=T,row.names=1,sep='\t')

FCUT=1.5
PCUT=0.05

#TEM=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PADJ< PCUT )]
#TCM=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PADJ< PCUT )]

#CD73POS=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PADJ < PCUT )]
#CD73NEG=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PADJ < PCUT )]

TEM=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PADJ< PCUT )]
TCM=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PADJ< PCUT )]

CD73POS=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PADJ < PCUT )]
CD73NEG=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PADJ < PCUT )]


#install.packages("VennDiagram")
library(VennDiagram)

x=list(TEM=TEM, TCM=TCM, CD73POS=CD73POS, CD73NEG=CD73NEG)
venn.diagram(x, filename='VENN_FOLD_1.5.tiff', height = 3000, width = 3000, resolution =
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
