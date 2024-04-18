library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(glue)
library(future)


if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

options(future.globals.maxSize = 64000 * 1024^2)



# read in 2 seurat objects under consideration
jayne <- LoadH5Seurat('jayne_dat/20230130_JEW_IntegratedSeurat_AnnotatedClusters_PhyloOrder.h5seurat')




### subset to just neutrophil clusters



jayne@meta.data$celltype %>% table() %>% enframe() %>% 
  ggplot(aes(x=name, y=value)) + geom_col()
jayne@meta.data$phyloorder
jayne@meta.data %>% colnames()


granulocytes <- subset(jayne, subset = celltype == 'granulocyte')

#
#######################


# remove bad cells
# seu_filt <- subset(seu, subset = REMOVE == 'KEEP')

Assays(granulocytes)
DefaultAssay(object = granulocytes) <- "RNA"
# granulocytes@assays$SCT <- NULL
# granulocytes@assays$integrated <- NULL
# remove genes with no expression in remaining cells
# non_zero_Features <- names(which(!rowSums(granulocytes) == 0))

# granulocytes <- subset(granulocytes, features=non_zero_Features)

# remove genes expressed in fewer than 10 cells

# counts <- GetAssayData(object = granulocytes, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
# nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
# keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
# filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
# filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = granulocytes@meta.data)


# write out seurat object
seu_filt <- SeuratDisk::SaveH5Seurat(granulocytes, 'outputs/granulocytes', overwrite = TRUE)


