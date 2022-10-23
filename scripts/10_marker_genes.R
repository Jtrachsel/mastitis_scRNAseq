library(tidyverse)
library(Seurat)
library(SeuratDisk)


# right now this is integrated by sampleID with 86 dimensions used for clustering
# All.integrated <- LoadH5Seurat('outputs/All.integrated.h5seurat')

All.integrated <- LoadH5Seurat('outputs/All.integrated_30.h5seurat')


# USE MARKER GENES FROM TABLE
DOT_PLOT <- read_csv('raw_data/dot_plot_markers.csv') 
MARKERS <- read_tsv('outputs/marker_genes_long.tsv')

MARKERS <- 
  MARKERS %>%
  mutate(amatches=map(.x=gene_name, .f=~agrep(.x, rownames(All.integrated), value = T)), 
         matches=map(.x=paste0('^',gene_name, '$'), .f=~grep(.x, rownames(All.integrated), value = T)))

MARKERS_dot <- 
  DOT_PLOT %>%
  mutate(amatches=map(.x=gene_name, .f=~agrep(.x, rownames(All.integrated), value = T)), 
         matches=map(.x=paste0('^',gene_name, '$'), .f=~grep(.x, rownames(All.integrated), value = T)))


UN_matched_markers <- 
  MARKERS %>% 
  filter(map_lgl(.x=matches, .f=~identical(.x, character(0)) ))

# grep('BOLA', rownames(All.integrated), value=T)

# PAX5 == ENSBTAG00000012498
# cd16 synonym
# itgam = cd18
# CD31 = PECAM1
# CD64 = FCGR1A
# glut1 = SLC2A1
# SELL, CD62L, LAM1, LECAM1, LEU8, LNHR, LSEL, LYAM1, PLNHR, TQ1, selectin L

matched_markers <- 
  MARKERS %>% 
  filter(!map_lgl(.x=matches, .f=~identical(.x, character(0)) )) %>% 
  select(-amatches, -matches)

#
feature_plots <- 
  matched_markers %>%
  group_by(type) %>% 
  summarise(marker_genes=list(gene_name)) %>% 
  mutate(feature_plot=map(.x=marker_genes, .f=~FeaturePlot(All.integrated, features = .x)))

# Bcell markers
feature_plots$feature_plot[[1]]
# XBP1 is a bad bcell marker
# PRDM1 is a bad bcell marker


#monocyte markers
feature_plots$feature_plot[[2]]
# remove CCR2

#myeloid markers
feature_plots$feature_plot[[3]]
# CD93
# FLT3


#neutrophil markers
feature_plots$feature_plot[[4]]

# MPO
# FCGR1A

#Other markers
feature_plots$feature_plot[[5]]
# both AHSP and HBM are bad


# T cell markers
feature_plots$feature_plot[[6]]
# all good
# matched_markers$type %>% unique()

matched_markers_filt <- 
  matched_markers %>%
  filter(!(gene_name %in% c('XBP1', 'PRDM1') & type == 'B/ASC')) %>% 
  filter(type != 'Others') %>% 
  filter(!(type == 'neutrophil' & gene_name %in% c('MPO', 'FCGR1A'))) %>% 
  filter(!(type == 'myeloid' & gene_name %in% c('CD93', 'FLT3'))) %>% 
  filter(!(type == 'monocyte' & gene_name %in% c('CCR2')))


matched_markers_filt %>% write_tsv('outputs/matched_markers_filt.tsv')

feature_plots_filt <- 
  matched_markers_filt %>%
  group_by(type) %>% 
  summarise(marker_genes=list(gene_name)) %>% 
  mutate(feature_plot=map(.x=marker_genes, .f=~FeaturePlot(All.integrated, features = .x)))

# B
feature_plots_filt$feature_plot[[1]]

#monocyte
feature_plots_filt$feature_plot[[2]]

# myeloid
feature_plots_filt$feature_plot[[3]]

#neutrophil
feature_plots_filt$feature_plot[[4]]

# Tcell
feature_plots_filt$feature_plot[[5]]



##
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5742132/

number_cells_expressing_markers <- 
  tibble(gene_name=names(rowSums(All.integrated[matched_markers_filt$gene_name,]@assays$originalexp@counts > 0)),
         num_positive=rowSums(All.integrated[matched_markers_filt$gene_name,]@assays$originalexp@counts > 0))

## make a kable table from this
number_cells_expressing_markers %>% arrange((num_positive))

# THIS IS TCRD
# grep('ENSBTAG00000055197', All_genes$ensembl_gene_id)


matched_markers


### find all positive markers for each cluster
Idents(All.integrated) <- 'integrated_snn_res.0.5'

FOUND_MARKERS <- FindAllMarkers(All.integrated, assay = 'originalexp',only.pos = TRUE )

FOUND_MARKERS %>% write_tsv('outputs/FOUND_MARKERS.tsv')

