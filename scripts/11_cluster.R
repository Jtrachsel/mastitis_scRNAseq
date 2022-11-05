library(tidyverse)
library(Seurat)
library(future)
library(glue)


if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

options(future.globals.maxSize = 200000 * 1024^2)

## functions 


# this was taken from one of Jayne's script, 2 methods to determine
# an appropriate number of dimensions to use in clustering and visualization
determine_dimensionality <- 
  function(seurat_object){
    pct <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100 # find standard deviation for each PC
    cumu <- cumsum(pct) # find cumulative percentages for PCs
    co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
    co1 # list PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
    co2 # list PC
    pcs <- c(co1, co2)
    return(pcs) 
  }

# runs clustering and visualization with the specified number of dimensions
run_cluster_viz <- function(seurat_obj, PCdims){
  seurat_obj <- RunTSNE(seurat_obj, reduction='pca', dims=PCdims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = PCdims)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = PCdims)
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = c(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                             verbose = FALSE) 
  return(seurat_obj)
}

# for a given resolution, returns the number of clusters present
# used to investigate different clustering results
num_clusters_by_resolution <- 
  function(seurat_object){
    seurat_object@meta.data %>% 
      select(starts_with('integrated')) %>% 
      rownames_to_column(var='cell_id') %>% 
      as_tibble() %>% 
      pivot_longer(-cell_id,
                   names_to = 'resolution', 
                   values_to = 'cluster', 
                   names_prefix = 'integrated_snn_res.') %>% 
      group_by(resolution) %>% 
      summarise(num_clusters=length(unique(cluster)))  
  }


library(SeuratDisk)
#

## two different seurat objects, one integrated by sample, one integrated by tissue
# I think integrated by sample_ID is probably the way to go
integrated_by_sample_ID <- LoadH5Seurat('outputs/split_by_sample_ID.h5seurat')
#read_rds('outputs/split_by_sample_ID.rds')
integrated_by_tissue <- LoadH5Seurat('outputs/split_by_tissue.h5seurat')

#
DefaultAssay(integrated_by_sample_ID) <- "integrated"
DefaultAssay(integrated_by_tissue) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated_by_sample_ID <- ScaleData(integrated_by_sample_ID, verbose = TRUE)
integrated_by_sample_ID <- RunPCA(integrated_by_sample_ID, npcs = 100, verbose = TRUE)

# Run the standard workflow for visualization and clustering
integrated_by_tissue <- ScaleData(integrated_by_tissue, verbose = TRUE)
integrated_by_tissue <- RunPCA(integrated_by_tissue, npcs = 100, verbose = TRUE)

# Method for determining dimensionality from Jayne's scripts
pcs_ID <- determine_dimensionality(integrated_by_sample_ID)
pcs_tissue <- determine_dimensionality(integrated_by_tissue)

# look at various cutoffs on the elbow plots

ElbowPlot(integrated_by_sample_ID, ndims = 100) +
  geom_vline(xintercept = c(pcs_ID[1],30,pcs_ID[2])) + 
  ylim(0,15)+
  annotate(geom = 'label', x=pcs_ID[2],y=10, label=pcs_ID[2])+
  annotate(geom = 'label', x=30,y=10, label=30)+
  annotate(geom = 'label', x=pcs_ID[1],y=10, label=pcs_ID[1])+
  ggtitle('dimensionality choices for integrated by sample')

ggsave(filename = 'outputs/figures/int_by_samp_elbow.jpeg', width = 7, height = 5, units = 'in')

ElbowPlot(integrated_by_tissue, ndims = 100) +
  geom_vline(xintercept = c(pcs_tissue[1],30,pcs_tissue[2])) +
  annotate(geom = 'label', x=pcs_tissue[2],y=10, label=pcs_tissue[2])+
  annotate(geom = 'label', x=30,y=10, label=30)+
  annotate(geom = 'label', x=pcs_tissue[1],y=10, label=pcs_tissue[1])+
  ylim(0,15)+
  ggtitle('dimensionality choices for integrated by tissue')

ggsave(filename = 'outputs/figures/int_by_tissue_elbow.jpeg', width = 7, height = 5, units = 'in')

# run clustering and vis on the various dimensionality choices
RESULTS <- 
  tribble(~input_data, ~seurat_objects, 
        glue('integrated_by_sample_ID_{pcs_ID[1]}'), run_cluster_viz(integrated_by_sample_ID, PCdims = 1:pcs_ID[1]),
        glue('integrated_by_sample_ID_{pcs_ID[2]}'), run_cluster_viz(integrated_by_sample_ID, PCdims = 1:pcs_ID[2]),
        glue('integrated_by_tissue_{pcs_tissue[1]}'), run_cluster_viz(integrated_by_sample_ID, PCdims = 1:pcs_tissue[1]),
        glue('integrated_by_tissue_{pcs_tissue[2]}'), run_cluster_viz(integrated_by_sample_ID, PCdims = 1:pcs_tissue[2]),
        glue('integrated_by_sample_ID_30'), run_cluster_viz(integrated_by_sample_ID, PCdims = 1:30),
        glue('integrated_by_tissue_30'), run_cluster_viz(integrated_by_tissue, PCdims = 1:30),
        )

# make some plots describing the number of clusters with different methods
#HERE
RESULTS <- 
  RESULTS %>% 
  mutate(num_clusts_dat=map(seurat_objects, ~num_clusters_by_resolution(.x)))

# plot the number of clusters formed with the different integration and
# dimensionality choices

RESULTS %>% 
  select(-seurat_objects) %>% 
  unnest(num_clusts_dat) %>% 
  ggplot(aes(x=resolution, y=num_clusters, color=input_data)) +
  geom_point() +
  geom_line(aes(group=input_data))# +
ggsave('outputs/figures/num_clusters_by_resolution.jpeg', width = 7, height = 5, units = 'in')
### USE THIS FIG ###

# generate DimPlots for choices of integration and dimensionality
RESULTS <- 
  RESULTS %>%
  mutate(dimplots_clusters=map(.x=seurat_objects, 
                      .f=~(DimPlot(.x, reduction='umap',group.by = 'integrated_snn_res.0.5', split.by = 'tissue'))), 
         dimplots_individual=map(.x=seurat_objects, 
                      .f=~(DimPlot(.x, reduction='umap',group.by = 'individual', split.by = 'tissue'))))

usethis::use_directory('outputs/figures')

# save figures
RESULTS %>% 
  select('input_data', 'dimplots_clusters', 'dimplots_individual') %>% 
  pivot_longer(-'input_data',names_to = 'type', values_to = 'plot' ) %>% 
  mutate(path=glue('outputs/figures/{input_data}_{type}.jpeg'), 
         save_res=map2(.x=plot,
                       .y=path,
                       .f=~ggsave(filename = .y, plot = .x, width = 9, height = 5, units = 'in')))


# using more dimensions gives more structure in the dimensionality reduction plots
# going with integrated by sampleid with 86 PCs used for clustering and viz
All.integrated <- RESULTS$seurat_objects[[1]]
DefaultAssay(All.integrated) <- "originalexp"

# integrated has the dim reductions, so cant remove?
# All.integrated[['integrated']] <- NULL
All.integrated[['SCT']] <- NULL


All.integrated <- NormalizeData(All.integrated,  
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "originalexp")

All.integrated <- ScaleData(All.integrated, 
                            assay = "originalexp")

# any(All.integrated@assays$originalexp@counts < 0)

SaveH5Seurat(All.integrated, filename = 'outputs/All.integrated_86', overwrite = TRUE)

### NOW FOR 30 DIMS

All.integrated <- RESULTS$seurat_objects[[5]]

# integrated has the dim reductions so cant remove?
# All.integrated[['integrated']] <- NULL
All.integrated[['SCT']] <- NULL

DefaultAssay(All.integrated) <- "originalexp"

All.integrated <- NormalizeData(All.integrated,  
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "originalexp")

All.integrated <- ScaleData(All.integrated, 
                            assay = "originalexp")

# any(All.integrated@assays$originalexp@counts < 0)

SaveH5Seurat(All.integrated, filename = 'outputs/All.integrated_30', overwrite=TRUE)

All.integrated

