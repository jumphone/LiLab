setwd('C:/Users/cchmc/Desktop/Lingjie')

MAT=read.table('DATA.txt',header=T,row.names=1,sep='\t')


FCUT=1.5
PCUT=0.05



TEM_GENE=rownames(MAT)[which( MAT$TEMdivTCM > FCUT & MAT$TCMTEM_PADJ< PCUT )]
TCM_GENE=rownames(MAT)[which( MAT$TCMdivTEM > FCUT & MAT$TCMTEM_PADJ< PCUT )]

POS_GENE=rownames(MAT)[which( MAT$POSdivNEG > FCUT & MAT$CD73_PADJ < PCUT )]
NEG_GENE=rownames(MAT)[which( MAT$NEGdivPOS > FCUT & MAT$CD73_PADJ < PCUT )]


