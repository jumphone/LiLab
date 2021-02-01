
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

DIFF_TCMTEM_PADJ=DIFF_TCMTEM[,19]
DIFF_CD73_PADJ=DIFF_CD73[,13]

DIFF_TCMTEM_TCMdivTEM=apply(DIFF_TCMTEM[,2:4],1,mean) / apply(DIFF_TCMTEM[,5:7],1,mean)
DIFF_TCMTEM_TEMdivTCM=apply(DIFF_TCMTEM[,5:7],1,mean) / apply(DIFF_TCMTEM[,2:4],1,mean)

DIFF_CD73_POSdivNEG=


