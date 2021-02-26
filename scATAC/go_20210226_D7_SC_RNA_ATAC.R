BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


/home/toolkit/tools/R4.0.3/bin/R


#######################################
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
############################

r_counts = Read10X_h5(filename = "RNA_filtered_feature_bc_matrix.h5")



##############################################
# ATAC
###############################################



a_counts = Read10X_h5(filename = "ATAC_filtered_peak_bc_matrix.h5")


a_metadata <- read.csv(
  file = "ATAC_singlecell.csv",
  header = TRUE,
  row.names = 1
)

a_chrom_assay <- CreateChromatinAssay(
  counts = a_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'ATAC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

a_pbmc <- CreateSeuratObject(
  counts = a_chrom_assay,
  assay = "peaks",
  meta.data = a_metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(a_pbmc) <- annotations

# compute nucleosome signal score per cell
a_pbmc <- NucleosomeSignal(object = a_pbmc)

# compute TSS enrichment score per cell
a_pbmc <- TSSEnrichment(object = a_pbmc, fast = FALSE)


# add blacklist ratio and fraction of reads in peaks
a_pbmc$pct_reads_in_peaks <- a_pbmc$peak_region_fragments / a_pbmc$passed_filters * 100
a_pbmc$blacklist_ratio <- a_pbmc$blacklist_region_fragments / a_pbmc$peak_region_fragments

##############################

saveRDS(a_pbmc, file='a_pbmc_beforeQC.rds')

#######################################
#a_pbmc=readRDS('a_pbmc_beforeQC.rds')


a_pbmc$high.tss <- ifelse(a_pbmc$TSS.enrichment > 2, 'High', 'Low')
a_pbmc$nucleosome_group <- ifelse(a_pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

TSSPlot(a_pbmc, group.by = 'high.tss') + NoLegend()
FragmentHistogram(object = a_pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = a_pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


a_pbmc <- subset(
  x = a_pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
a_pbmc

saveRDS(a_pbmc, file='a_pbmc_afterQC.rds')

#######################################

a_pbmc <- RunTFIDF(a_pbmc)
a_pbmc <- FindTopFeatures(a_pbmc, min.cutoff = 'q0')
a_pbmc <- RunSVD(a_pbmc)

DepthCor(a_pbmc)

a_pbmc <- RunUMAP(object = a_pbmc, reduction = 'lsi', dims = 2:30)
DimPlot(object = a_pbmc, label = TRUE) + NoLegend()

saveRDS(a_pbmc, file='a_pbmc_afterUMAP.rds')

###################################



pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
a_pbmc <- AddMotifs(
  object = a_pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)


a_chromvar=RunChromVAR(object = a_pbmc[["peaks"]], genome = BSgenome.Hsapiens.UCSC.hg38)
















a_pbmc <- FindNeighbors(object = a_pbmc, reduction = 'lsi', dims = 2:30)
a_pbmc <- FindClusters(object = a_pbmc, verbose = FALSE, algorithm = 3)















