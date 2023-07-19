


library(ArchR)
set.seed(1)

addArchRThreads(threads = 20)
addArchRGenome("hg19")
proj <- loadArchRProject(path = "ArchRSubset")

TMP=getMatrixFromProject(proj,useMatrix ='PeakMatrix')
CHR=TMP@rowRanges@seqnames
START=TMP@rowRanges@ranges@start
END=TMP@rowRanges@ranges@start+TMP@rowRanges@ranges@width-1
PEAK_MAT=TMP@assays@data$PeakMatrix
colnames(PEAK_MAT)=stringr::str_replace_all(rownames(TMP@colData),'#','_')
rownames(PEAK_MAT)=paste0(CHR,':',START,'-',END)

counts=PEAK_MAT


