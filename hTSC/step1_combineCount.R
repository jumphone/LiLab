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
  'CT_second',  'EVT_TSCT27', 'CT_CT27-P11',
  'ST_TSCT27-2D', 'ST_CT27-ST-3D', 'EVT_TSCT29',
  'CT_CT29-P13', 'ST_CT29-ST-3D', 'Stroma_endometrium_3',
  'Stroma_endometrium_4', 'NA_WT-ES8P74-D3','CT_term-1-1',
   'CT_term-1-2','CT_term-2-1','CT_term-2-2',
  'CT_first-1-1','CT_first-1-2','CT_first-2-1',
  'CT_first-2-2','CT_second-1','CT_second-2',
  'CT_first-1','CT_first-2','Epi_secretory-1',
  'Epi_secretory-2','Epi_follicula-2','Stroma_endometrium-1',
   'Stroma_endometrium-2','CT_first-3-1','CT_first-3-2',
  'EVT_first-1','EVT_first-2',
  'EVT_TSblast-1','ST_TSblast-1-3D',
  'TS_blast-1','EVT_TSblast-2','ST_TSblast-2-3D',
  'TS_blast-2','Stroma_first-1','Stroma_first-2',
  'NA_V10-12h','NA_V10-18h','NA_V10-3h',
  'NA_V10-6h','NA_V10-D1','NA_V10-D2',
  'NA_V10-D3','NA_V1-D1','NA_V1-D2',
  'NA_V1-D3','NA_V7-12h','NA_V7-18h',
  'NA_V7-3h','NA_V7-6h', 'NA_WT-3h',
  'NA_WT-6h','NA_WT-DMSO-D3'  
 )


TYPE=c(
  'CT',  'EVT', 'CT',
  'ST', 'ST', 'EVT',
  'CT', 'ST', 'Stroma',
  'Stroma', 'NA','CT',
   'CT','CT','CT',
  'CT','CT','CT',
  'CT','CT','CT',
  'CT','CT','Epi',
  'Epi','Epi','Stroma',
   'Stroma','CT','CT',
  'EVT','EVT',
  'EVT','ST',
  'TS','EVT','ST',
  'TS','Stroma','Stroma',
  'NA','NA','NA',
  'NA','NA','NA',
  'NA','NA','NA',
  'NA','NA','NA',
  'NA','NA', 'NA',
  'NA','NA'  
 )
saveRDS(TYPE, 'TYPE.rds')
FULL.COL_NAME=NEW.COL_NAME #paste0(NEW.COL_NAME, '__',COL_NAME )

FULL.MAT=MAT
colnames(FULL.MAT)=FULL.COL_NAME
write.table(FULL.MAT, file='MAT_FULL.NAME.txt',sep='\t',quote=F,row.names=TRUE, col.names=TRUE)
#####################







