library(Seurat)
library(tidyverse)
library(future)
library(SeuratDisk)


set.seed(4)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

# reading regarding how to split samples when SCTransforming
# https://github.com/satijalab/sctransform/issues/55
# https://github.com/satijalab/sctransform/issues/32


# Hi Alex,
# Thank you for your interest in the package.
# Generally, I'd combine all samples across batches by simply merging.
# If common cell types separate due to batch effect,
# merge samples per batch and try an integration approach as outline here.

## split by tissue

seu_filt <- SeuratDisk::LoadH5Seurat('outputs/seurat_QC_done.h5seurat')

# https://satijalab.org/seurat/articles/integration_introduction.html

# recommended in the 'integration' tutorial
# but multiple paths, could split by tissue for example

# function to run integration on a seurat object, need to specify the metadata col
# to "split" on as a string.
RUN_integration <- function(SPLIT_BY, SEURAT){
  All.list <- SplitObject(SEURAT, split.by = SPLIT_BY)

  for (i in 1:length(All.list)) { # normalize data using SCTransform method
    All.list[[i]] <- SCTransform(All.list[[i]],
                                 assay='originalexp',
                                 return.only.var.genes = TRUE,
                                 variable.features.n = 7500,
                                 # variable.features.n = NULL,   # Null to use rv.th
                                 # variable.features.rv.th = 1.3, # 1.3 = default
                                 verbose = TRUE,
                                 n_genes=NULL , # use all genes for sctransform::vst
                                 n_cells=NULL # use all cells for sctransform::vst
    )
  }

  All.features <- SelectIntegrationFeatures(All.list,
                                            verbose = TRUE,
                                            nfeatures=5000) # select the genes to use for integration
  All.list <- PrepSCTIntegration(All.list,
                                 anchor.features = All.features,
                                 verbose = TRUE)

  # 50 looked like it worked well....
  All.anchors <- FindIntegrationAnchors(All.list,
                                        normalization.method = "SCT",
                                        anchor.features = All.features,
                                        dims = 1:50) # identify anchors for integration from top 50 data dimensions
  All.integrated <- IntegrateData(All.anchors,
                                  normalization.method = "SCT",
                                  dims = 1:50) # integrate the data
  DefaultAssay(All.integrated) <- "integrated"

  return(All.integrated)

}


## run integration splitting by sample_ID
sample_ID_integration <- RUN_integration(SPLIT_BY = 'sample_ID', SEURAT = seu_filt)
# remove SCT assay because we dont need it after integration
sample_ID_integration[['SCT']] <- NULL
DefaultAssay(sample_ID_integration) <- "integrated"
SeuratDisk::SaveH5Seurat(sample_ID_integration, 'outputs/split_by_sample_ID', overwrite = T)


# clean up environment
rm(sample_ID_integration)
gc()

## run integration splitting by tissue
tissue_integration <- RUN_integration(SPLIT_BY = 'tissue', SEURAT = seu_filt)
# remove SCT assay because we dont need it after integration
tissue_integration[['SCT']] <- NULL
DefaultAssay(tissue_integration) <- "integrated"
SeuratDisk::SaveH5Seurat(tissue_integration, 'outputs/split_by_tissue', overwrite = T)




