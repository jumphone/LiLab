BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")




/home/toolkit/tools/R4.0.3/bin/R

setwd("/home/database/data/D7_SC_RNA_ATAC")

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)


library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)






