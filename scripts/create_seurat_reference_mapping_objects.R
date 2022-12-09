# install.packages("anndata")
library(anndata)
library(tidyverse)

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(Matrix)


# some growing pains in the communtiy regarding anndata to seurat formats
# https://github.com/mojaveazure/seurat-disk/issues/109

tabula_blood_anndata <- anndata::read_h5ad('reference_mapping_data/tabula_10X_blood.h5ad')

tabula_blood <- CreateSeuratObject(counts = Matrix::t(tabula_blood_anndata$X), meta.data = tabula_blood_anndata$obs)

non_zero_genes <- rowSums(tabula_blood) > 0
tabula_blood_filt <- tabula_blood[non_zero_genes,]

SaveH5Seurat(tabula_blood_filt, 'reference_mapping_data/tabula_sapiens_blood')


### MIT human breast milk study 

# uses weird genenames
# some are gene symbols and some are embl ids?

Nyquist_counts <- read_tsv('reference_mapping_data/MIT_Milk_Study_Raw_counts.txt') %>%
  column_to_rownames(var = 'GENE') %>% as.matrix() %>% Matrix(sparse=TRUE)



Nyquist_meta <- read_csv('reference_mapping_data/MIT_milk_study_metadata.csv')
# the first row is information about the type of each column, remove it here
Nyquist_meta <- Nyquist_meta[-1,] %>% column_to_rownames(var = 'NAME')

Nyquist_milk <- Seurat::CreateSeuratObject(counts=Nyquist_counts, meta.data = Nyquist_meta)


SaveH5Seurat(Nyquist_milk, 'reference_mapping_data/Nyquist_milk')
seu <- LoadH5Seurat('reference_mapping_data/Nyquist_milk.h5seurat')

detected_embl <- sub('(.*)\\.[0-9]+','\\1',rownames(seu))



look <- embl_mapping %>% filter(embl %in% detected_embl) 


detected_embl[!(detected_embl %in% look$embl)]




