library(Seurat)
library(SeuratDisk)
library(tictoc)
library(tidyverse)

tic()
seu <- LoadH5Seurat('outputs/classified_clusters_30.h5seurat')
project_time <- toc()



seu@meta.data$tissue_clust <- 
  paste(seu@meta.data$dot_plot_group,
        seu@meta.data$integrated_snn_res.0.5, 
        seu@meta.data$tissue,
        sep = '_')

seu@meta.data$class_clust <- 
  paste(seu@meta.data$dot_plot_group,
        seu@meta.data$integrated_snn_res.0.5,
        sep = '_')


# check current cell identiies
Idents(seu)

# set cell identities for seurat object
Idents(seu) <- seu@meta.data$class_clust




# this function is a wrapper for Seurat FindMarkers that
# renames the results and adds a column indicating which tissue
# the marker is enriched in.
get_tissue_markers <-
  function(seurat_object, SUBSET_IDENT){
    
    FindMarkers(seu, densify = TRUE, ident.1 = 'milk', subset.ident = SUBSET_IDENT, group.by = 'tissue')%>%
    rownames_to_column('gene') %>% 
    transmute(gene, p_val_adj, avg_log2FC,
              pct_milk=pct.1,
              pct_blood=pct.2,
              enriched_in = ifelse(avg_log2FC > 0 , 'milk', 'blood')) %>% 
    arrange(desc(avg_log2FC)) %>% 
      filter(p_val_adj < 0.05)
    
    
  }



# find tissue markers for each cluster (lowest resolution)
tissue_markers <- 
  seu@meta.data %>% 
  group_by(class_clust) %>% 
  tally() %>% 
  mutate(tissue_markers = map(.x=class_clust,
                              .f=~get_tissue_markers(seurat_object = seu,
                                                     SUBSET_IDENT = .x)))

# write out a table for each cluster
system('mkdir outputs/tissue_marker_genes')

tissue_markers %>% 
  mutate(PATH=paste0('outputs/tissue_marker_genes/', class_clust, '.tsv')) %>% 
  mutate(map2(.x=PATH,
              .y=tissue_markers,
              .f=~write_tsv(file =.x, x = .y )))



seu@meta.data
