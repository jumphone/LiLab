SAM='./collect_pair_tag.interaction.sam'
###############
DATA=data.table::fread(paste0(SAM,'.bam.bedpe'),sep='\t')
CHR=cbind(DATA[,1],DATA[,4])
MAT=cbind(DATA[,2],DATA[,3],DATA[,5],DATA[,6])

getL=function(X){
    A=X[1]
    B=X[2]
    C=X[3]
    D=X[4]
    if(A<C){
        L=C-A-(B-A)
        }else{
        L=A-C-(D-C)
        }
    return(L)
    }

LLL=apply(MAT,1,getL)

OUT=DATA[which(LLL>1 & CHR[,1]==CHR[,2]),]

write.table(OUT,file='./collect_pair_tag.interaction.sam.filter.txt',sep='\t',quote=F,row.names=F,col.names=F)










