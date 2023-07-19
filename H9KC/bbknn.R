
reticulate::use_python('/home/toolkit/local/bin/python3',required=TRUE)

.bbknn<-function(PCA, BATCH,NB=3, NT=10, DM=2){
    NB=NB
    NT=NT
    DM=DM
    pca.use=PCA
    batch=BATCH
    #mybeer=mybeer
    library(reticulate)
    #use_python("C:\Users\cchmc\Anaconda3\python")
    anndata = reticulate::import("anndata",convert=FALSE) #anndata==0.7
    bbknn = reticulate::import("bbknn", convert=FALSE) #1.5.1
    sc = reticulate::import("scanpy",convert=FALSE) #scanpy
    adata = anndata$AnnData(cbind(pca.use,pca.use), obs=batch)
    PCNUM=ncol(pca.use)
    sc$tl$pca(adata, n_comps=as.integer(PCNUM))
    adata$obsm$X_pca = pca.use
    bbknn$bbknn(adata,batch_key=0,neighbors_within_batch=as.integer(NB),n_pcs=as.integer(PCNUM), annoy_n_trees =as.integer(NT))
    sc$tl$umap(adata, n_components=as.integer(DM))
    umap = py_to_r(adata$obsm['X_umap'])
    rownames(umap)=rownames(PCA)
    colnames(umap)=paste0('UMAP_',c(1:ncol(umap)))#colnames(pbmc@reductions$umap@cell.embeddings)
    return(umap)
    }

