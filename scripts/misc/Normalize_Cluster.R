library(Seurat)
library(tidyverse)

# install glmGamPoi
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("glmGamPoi")
# install sctransform from Github
# remotes::install_github("satijalab/sctransform", ref = "develop")



# reading regarding how to split samples when SCTransforming
# https://github.com/satijalab/sctransform/issues/55
# https://github.com/satijalab/sctransform/issues/32


# Hi Alex,
# Thank you for your interest in the package.
# Generally, I'd combine all samples across batches by simply merging.
# If common cell types separate due to batch effect,
# merge samples per batch and try an integration approach as outline here.

## split by tissue


##
seu_filt <- read_rds('outputs/seurat_QC_done.rds')



# https://satijalab.org/seurat/articles/integration_introduction.html
# WHY DOES JAYNE SPLIT THESE ?????

# recommended in the 'integration' tutorial
# but multiple paths, could split by tissue for example
seu_filt@meta.data$sample_ID %>% unique()

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
  # trying 100
  All.anchors <- FindIntegrationAnchors(All.list, 
                                        normalization.method = "SCT", 
                                        anchor.features = All.features, 
                                        dims = 1:50) # identify anchors for integration from top 30 data dimensions
  All.integrated <- IntegrateData(All.anchors, 
                                  normalization.method = "SCT", 
                                  dims = 1:50) # integrate the data
  DefaultAssay(All.integrated) <- "integrated"
  
  
  
  
  return(All.integrated)
  
}
  ### PASTE
  DefaultAssay(All.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  All.integrated <- ScaleData(All.integrated, verbose = TRUE)
  All.integrated <- RunPCA(All.integrated, npcs = 100, verbose = TRUE)
  All.integrated <- JackStraw(All.integrated, num.replicate = 100)
  All.integrated <- ScoreJackStraw(All.integrated, dims = 1:100)
  
  jack_straw_plot <- JackStrawPlot(All.integrated, dims = 1:15)
  # 
  # pct <- All.integrated[["pca"]]@stdev / sum(All.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
  # cumu <- cumsum(pct) # find cumulative percentages for PCs
  # co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
  # # co1 # list PC
  # co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
  # # co2 # list PC
  # pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
  # print(pcs)
  # PCdims <- 1:pcs 
  All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:50)
  All.integrated <- FindNeighbors(All.integrated, reduction = "pca", dims = 1:50)
  # All.integrated <- FindClusters(All.integrated, resolution = 0.5)
  All.integrated <- FindClusters(All.integrated,
                                 resolution = c(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                                 verbose = FALSE) 
  
return(All.integrated)


sample_ID_integration <- RUN_integration(SPLIT_BY = 'sample_ID', SEURAT = seu_filt)
write_rds(sample_ID_integration, 'outputs/split_by_sample_ID.rds')

rm(sample_ID_integration)
gc()

tissue_integration <- RUN_integration(SPLIT_BY = 'tissue', SEURAT = seu_filt)
write_rds(tissue_integration, 'outputs/split_by_tissue.rds')


# END FUNCTION

# 
# All.list <- SplitObject(seu_filt, split.by = "sample_ID") 
# 
# for (i in 1:length(All.list)) { # normalize data using SCTransform method
#   All.list[[i]] <- SCTransform(All.list[[i]], 
#                                assay='originalexp',
#                                return.only.var.genes = TRUE, 
#                                variable.features.n = 7500,
#                                # variable.features.n = NULL,   # Null to use rv.th
#                                # variable.features.rv.th = 1.3, # 1.3 = default
#                                verbose = TRUE, 
#                                n_genes=NULL , # use all genes for sctransform::vst
#                                n_cells=NULL # use all cells for sctransform::vst
#   ) 
# }
# 
# All.features <- SelectIntegrationFeatures(All.list, 
#                                           verbose = TRUE, 
#                                           nfeatures=7500) # select the genes to use for integration
# All.list <- PrepSCTIntegration(All.list, 
#                                anchor.features = All.features,
#                                verbose = TRUE)
# 
# # 50 looked like it worked well....
# # trying 100
# All.anchors <- FindIntegrationAnchors(All.list, 
#                                       normalization.method = "SCT", 
#                                       anchor.features = All.features, 
#                                       dims = 1:100) # identify anchors for integration from top 30 data dimensions
# All.integrated <- IntegrateData(All.anchors, 
#                                 normalization.method = "SCT", 
#                                 dims = 1:100) # integrate the data





### PASTE
# DefaultAssay(sample_ID_integration) <- "integrated"
# 
# # Run the standard workflow for visualization and clustering
# All.integrated <- ScaleData(All.integrated, verbose = TRUE)
# All.integrated <- RunPCA(All.integrated, npcs = 100, verbose = TRUE)
# All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:100)
# All.integrated <- FindNeighbors(All.integrated, reduction = "pca", dims = 1:100)
# All.integrated <- FindClusters(All.integrated, resolution = 0.5)
# 
# write_rds(All.integrated, 'outputs/split_by_sample_ID.rds')
# 
# 




p1 <- DimPlot(tissue_integration, reduction = "umap", group.by = "tissue")
p2 <- DimPlot(sample_ID_integration, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(All.integrated, reduction = "umap", split.by = "tissue")


DefaultAssay(All.integrated) <- "originalexp"
nk.markers <- FindConservedMarkers(All.integrated, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)


# USE MARKER GENES FROM TABLE
FeaturePlot(All.integrated, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                          "CCL2", "PPBP"), min.cutoff = "q9")



### END PASTE
### SPLIT INVESTIGATION FROM THE SCTRANSFORM VIGNET
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
# InstallData("ifnb")
# load dataset
# ifnb <- LoadData("ifnb")

# ifnb@meta.data$orig.ident %>% table()


# split the dataset into a list of two seurat objects (stim and CTRL)
# ifnb.list <- SplitObject(ifnb, split.by = "stim")

# ctrl <- ifnb.list[["CTRL"]]
# stim <- ifnb.list[["STIM"]]




### traditional clustering route

seu_filt <- NormalizeData(seu_filt, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)

# ID highly variable features
seu_filt <- FindVariableFeatures(seu_filt, selection.method = "vst", nfeatures = 7500)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu_filt), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu_filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# scale data
all.genes <- rownames(seu_filt)
seu_filt <- ScaleData(seu_filt, features = all.genes)


seu_filt <- RunPCA(seu_filt, features = VariableFeatures(object = seu_filt))

# lots of ribosomal proteins for PC1....
VizDimLoadings(seu_filt, dims = 1:2, reduction = "pca")


DimPlot(seu_filt, reduction = "pca", group.by = 'individual', split.by = 'tissue')

# seu_filt@reductions$pca
DimHeatmap(seu_filt, dims = 1:6, cells = 500, balanced = TRUE)


seu_filt <- JackStraw(seu_filt, dims=50,num.replicate = 100)
seu_filt <- ScoreJackStraw(seu_filt, dims = 1:50)


JackStrawPlot(seu_filt, dims = 20:50)
ElbowPlot(seu_filt)

###

seu_filt <- FindNeighbors(seu_filt, dims = 1:50)



seu_filt <- FindClusters(seu_filt, resolution = 0.95)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.95=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.95



seu_filt <- FindClusters(seu_filt, resolution = 0.85)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.85=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.85





seu_filt <- FindClusters(seu_filt, resolution = 0.75)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.75=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.75



seu_filt <- FindClusters(seu_filt, resolution = 0.5)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.5=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.5

seu_filt <- FindClusters(seu_filt, resolution = 0.25)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.25=seurat_clusters)

seu_filt <- FindClusters(seu_filt, resolution = 0.125)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.125=seurat_clusters)

seu_filt <- FindClusters(seu_filt, resolution = 0.0625)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.0625=seurat_clusters)

seu_filt <- RunUMAP(seu_filt, dims = 1:50)
DimPlot(seu_filt, reduction = "umap", group.by = 'seurat_clusters_0.95', label = T)
DimPlot(seu_filt, reduction = "umap", group.by = 'seurat_clusters_0.95', label = T, split.by = 'tissue')
#####
# pretty much all clusters are exclusive to either milk or blood
seu_filt@meta.data %>% 
  group_by(seurat_clusters_0.95, tissue) %>% 
  tally() %>% 
  ggplot(aes(y=n, x=tissue, fill=tissue)) +
  geom_col(position = position_dodge(preserve = 'single'), color='black') + 
  facet_wrap(~seurat_clusters_0.95, nrow=1) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
