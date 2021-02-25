

 ‘GenomeInfoDb’, ‘GenomicRanges’, ‘IRanges’, ‘Rsamtools’, ‘S4Vectors’, ‘BiocGenerics’, ‘Biostrings’, ‘ggbio’, ‘biovizBase’, ‘AnnotationFilter’

BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")

BiocManager::install("Rsamtools")
BiocManager::install("S4Vectors")
BiocManager::install("BiocGenerics")

BiocManager::install("Biostrings")
BiocManager::install("ggbio")
BiocManager::install("biovizBase")
BiocManager::install("AnnotationFilter")



BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("patchwork")



/home/toolkit/tools/R4.0.3/bin/R

setwd('/home/database/data/TMP_scATAC')

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=T)
py_config()

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)



counts <- Read10X_h5(filename = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)


chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

# add the gene information to the object
Annotation(pbmc) <- annotations




# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')


pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc




library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pfm
)


pbmc_chromvar=RunChromVAR(object = pbmc[["peaks"]], genome = BSgenome.Hsapiens.UCSC.hg19)


