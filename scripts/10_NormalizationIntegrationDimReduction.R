#20221129 - Jayne Wiarda
# Integration, dimensionality reduciton, normalization

library(Seurat)
library(future)
library(SeuratDisk)


set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20221129_JEW_FilteredSeurat.h5seurat')
seu

seu.list <- SplitObject(seu, split.by = "sample_ID") # split by sample IDs

## Perform SCTransform normalization of data from each sample:
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

## Integrate the data from different samples:
seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top PCs
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT")

## Run multidimensional analyses on data:
seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 100 dimensions of the data
                         npcs = 100, 
                         verbose = TRUE) 
ElbowPlot(seu.integrated,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering

## Perform multidimensional visualization of data:
seu.integrated <- RunUMAP(seu.integrated, 
                          dims = 1:30, # use our calculated number of PCs
                          reduction = "pca", 
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          min.dist = 0.5,
                          spread = 0.2,
                          assay = "SCT") # create UMAP
seu.integrated <- RunTSNE(seu.integrated, 
                          dims = 1:30, 
                          reduction = "pca", 
                          assay = "SCT") # create tSNE plot 

## Also define nearest neighbors:
seu.integrated <- FindNeighbors(seu.integrated, 
                                dims = 1:30, 
                                verbose = TRUE) 

## Add normalized/scaled data to RNA assay:
#dim(seu.integrated[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA")
#dim(seu.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now

DefaultAssay(seu.integrated) <- "RNA"
SaveH5Seurat(seu.integrated, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20221129_JEW_IntegratedSeurat.h5seurat', overwrite = TRUE)

