
/home/toolkit/tools/R4.0.3/bin/R


#######################################
setwd("/home/zhangfeng/project/CombineRnaAtac/data/10x_multiomic")

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
set.seed(1234)
