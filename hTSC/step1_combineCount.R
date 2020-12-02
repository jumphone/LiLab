#source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
source('/home/zhangfeng/project/EVOL/source/MultiTools.R')

setwd('/home/zhangfeng/project/hTSC/data/hTSC_paper_RNA_seq/htseq_count')

FILE_NAME=as.character(read.table('COUNT_FILE_NAME.txt',sep='\t',header=FALSE)[,1])

i=1
this_file_name=FILE_NAME[i]
this_count=read.table(this_file_name, row.names=1,header=TRUE,sep='\t') 
this_count=cbind(this_count, this_count)
this_count=this_count[1:(nrow(this_count)-5),]
COM=this_count

i=2
while(i<=length(FILE_NAME)){
    this_file_name=FILE_NAME[i]
    this_count=read.table(this_file_name, row.names=1,header=TRUE,sep='\t') 
    this_count=cbind(this_count, this_count)
    this_count=this_count[1:(nrow(this_count)-5),]
    this_com=.simple_combine(COM, this_count, FILL=TRUE)
    COM=this_com$combine
    print(this_file_name)
    print(dim(COM))
    print(i)
    i=i+1}


MAT=COM[,(1:length(FILE_NAME)) *2]
COL_NAME=str_replace(FILE_NAME,pattern='.unique.sort.dedupped.gene_count.txt','')
COL_NAME=str_replace(COL_NAME,pattern='_RNAseq','')
#COL_NAME=str_replace(COL_NAME,pattern='RNA_','')
#COL_NAME=str_replace(COL_NAME,pattern='-RNA','')

colnames(MAT)=COL_NAME

write.table(MAT, file='MAT_ORI.NAME.txt',sep='\t',quote=F,row.names=TRUE, col.names=TRUE)
#################

colnames(MAT)

NEW.COL_NAME=c(
  'CT_second',  'EVT_TSCT27', 'NA_CT27-P11',
  'ST_TSCT27-2D', 'NA_CT27-ST-3D', 'EVT_TSCT29',
  'NA_CT29-P13', 'NA_CT29-ST-3D', 'Stroma_endometrium_3',
  'StromaEndometrium_4', 'NA_ES8P74-WTD3','CT_term-1-1',
   'CT_term-1-2','CT_term-2-1','CT_term-2-2',
  'CT_first-1-1','CT_first-1-2','CT_first-2-1',
  'CT_first-2-2','CT_second-1','CT_second-2',
  'CT_first-1','CT_first-2','Epi_secretory-1',
  'Epi_secretory-2','Epi_follicula-2','Stroma_endometrium-1',
   'Stroma_endometrium-2','CT_first-3-1','CT_first-3-2',
  'EVT_first-1','EVT_first-2','NA_JKU048',
  'NA_JKU049','EVT_TSblast-1','ST_TSblast-1-3D',
  'TS_blast-1','EVT_TSblast-2','ST_TSblast-2-3D',
  'TS_blast-2','Stroma_first-1','Stroma_first-2',
  'NA_V1012h','NA_V1018h','NA_V103h',
  'NA_V106h','NA_V10-D1-6','NA_V10-D2-6',
  'NA_V10-D3-6','NA_V1-D1-5','NA_V1-D2-5',
  'NA_V1-D3-5','NA_V712h','NA_V718h',
  'NA_V73h','NA_V76h', 'NA_WT3h',
  'NA_WT6h','NA_WT-DMSOD3'  
 )

FULL.COL_NAME=paste0(NEW.COL_NAME, '__',COL_NAME )

FULL.MAT=MAT
colnames(FULL.MAT)=FULL.COL_NAME
write.table(FULL.MAT, file='MAT_FULL.NAME.txt',sep='\t',quote=F,row.names=TRUE, col.names=TRUE)
#####################







