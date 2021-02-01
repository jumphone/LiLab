
setwd('E:/Lilab/Project/Fengqin')

COUNT_MAT_TCMTEM=read.table('INPUT_MAT_TCM_TEM.txt',header=TRUE,sep='\t',row.names=NULL)

USED_GENE=names(which(table(COUNT_MAT_TCMTEM[,1])==1))
COUNT_MAT_TCMTEM=COUNT_MAT_TCMTEM[which(COUNT_MAT_TCMTEM[,1] %in% USED_GENE),]
rownames(COUNT_MAT_TCMTEM)=COUNT_MAT_TCMTEM[,1]
COUNT_MAT_TCMTEM=COUNT_MAT_TCMTEM[,2:ncol(COUNT_MAT_TCMTEM)]

head(COUNT_MAT_TCMTEM)


library(DESeq2)
mydata=COUNT_MAT_TCMTEM
type<-factor(c(rep("TCM",3),rep("TEM",3)),levels = c("TCM","TEM"))
database<-round(as.matrix(mydata))
coldata<-data.frame(row.names = colnames(database),type)
dds<-DESeqDataSetFromMatrix(database,coldata,design = ~type)
dds<-DESeq(dds)
res = results(dds, contrast=c("type", "TCM", "TEM"))
res=as.matrix(res)

###############
.normRPM <- function(x){
    x.sum=sum(x)
    if(x.sum==0){
        y=x
    }else{
        y=x / x.sum *1000000
        }
    return(y)
    }
#####################

RPM_MAT_TCMTEM=apply(COUNT_MAT_TCMTEM,2,.normRPM)

mydata_rpm=RPM_MAT_TCMTEM
colnames(mydata_rpm)=paste0('RPM_',colnames(mydata_rpm))

###########################
all_res=cbind(mydata, mydata_rpm, res)

write.table(all_res,'OUTPUT_DIFF_TCMTEM.txt',row.names=T,col.names=T,quote=F,sep='\t')
########################################################################################





DIFF_TCMTEM=read.table('OUTPUT_DIFF_TCMTEM.txt',header=TRUE,sep='\t',row.names=NULL)
DIFF_CD73=read.table('INPUT_DIFF_CD73.txt',header=TRUE,sep='\t',row.names=NULL)

colnames(DIFF_TCMTEM)
colnames(DIFF_CD73)

DIFF_TCMTEM_GENE=DIFF_TCMTEM[,1]
DIFF_CD73_GENE=DIFF_CD73[,1]

DIFF_TCMTEM_PVALUE=DIFF_TCMTEM[,18]
DIFF_CD73_PVALUE=DIFF_CD73[,12]

DIFF_TCMTEM_PADJ=DIFF_TCMTEM[,19]
DIFF_CD73_PADJ=DIFF_CD73[,13]


DIFF_TCMTEM_RPM_TCM=apply(DIFF_TCMTEM[,2:4],1,mean)
DIFF_TCMTEM_RPM_TEM=apply(DIFF_TCMTEM[,5:7],1,mean)

DIFF_TCMTEM_TCMdivTEM = DIFF_TCMTEM_RPM_TCM / DIFF_TCMTEM_RPM_TEM
DIFF_TCMTEM_TEMdivTCM = DIFF_TCMTEM_RPM_TEM / DIFF_TCMTEM_RPM_TCM 


DIFF_CD73_FPKM_POS=apply(DIFF_CD73[,2:4],1,mean)
DIFF_CD73_FPKM_NEG=apply(DIFF_CD73[,5:7],1,mean)

DIFF_CD73_POSdivNEG= DIFF_CD73_FPKM_POS / DIFF_CD73_FPKM_NEG
DIFF_CD73_NEGdivPOS= DIFF_CD73_FPKM_NEG / DIFF_CD73_FPKM_POS


MAT1=cbind(DIFF_TCMTEM_RPM_TCM, DIFF_TCMTEM_RPM_TEM, DIFF_TCMTEM_TCMdivTEM, DIFF_TCMTEM_TEMdivTCM, DIFF_TCMTEM_PVALUE, DIFF_TCMTEM_PADJ)
rownames(MAT1)=DIFF_TCMTEM_GENE

MAT2=cbind(DIFF_CD73_FPKM_POS, DIFF_CD73_FPKM_NEG, DIFF_CD73_POSdivNEG, DIFF_CD73_NEGdivPOS, DIFF_CD73_PVALUE, DIFF_CD73_PADJ)
MAT2_USED_INDEX=which(DIFF_CD73_GENE %in% names(which(table(DIFF_CD73_GENE)==1)))
MAT2=MAT2[MAT2_USED_INDEX,]
rownames(MAT2)=DIFF_CD73_GENE[MAT2_USED_INDEX]

################

.simple_combine <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){    
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){ 
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(0,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(0,ncol=ncol(exp_ref_mat),nrow=length(gene21))
        rownames(exp_ref_mat_add)=gene21
        colnames(exp_ref_mat_add)=colnames(exp_ref_mat)
        exp_sc_mat=rbind(exp_sc_mat, exp_sc_mat_add)
        exp_ref_mat=rbind(exp_ref_mat, exp_ref_mat_add)
    }
    ############################################ 
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }
#################

COM=.simple_combine(MAT1,MAT2)$combine

write.table(COM,'OUTPUT_COMPARE.txt',row.names=T,col.names=T,quote=F,sep='\t')

