library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(glue)


# This script takes pre-identified marker genes and matches them to the
# genes detected in the Seurat object

# also makes feature plots of these marker genes and removes
# marker genes that are not very informative

# finally runs FindALlMarkers() and outputs a tsv with each cluster's positive marker genes.


# makes 'feature' plots using specified markers
make_feature_plots <- function(SEURAT, matched_markers,label ){
  # browser()
  feature_plots <-
    matched_markers %>%
    group_by(type) %>%
    summarise(marker_genes=list(gene_name)) %>%
    mutate(feature_plot=map(.x=marker_genes, .f=~FeaturePlot(SEURAT, features = .x)),
           plot_path=glue('outputs/figures/{label}_{type}_feature_plot.jpeg'),
           ggsave_res=map2(.x=feature_plot,
                           .y=plot_path,
                           .f=~ggsave(plot = .x, filename = .y, width=9, height=7, units='in', bg='white')))

  return(feature_plots)

}

# runs FindAllMarkers and outputs a tsv
get_cluster_markers <- function(SEURAT, label){

  output_name <- glue('outputs/{label}_markers.tsv')

  Idents(SEURAT) <- 'integrated_snn_res.0.5'
  FOUND_MARKERS <- FindAllMarkers(SEURAT, assay = 'originalexp',only.pos = TRUE )
  FOUND_MARKERS %>% write_tsv(output_name)
  return(FOUND_MARKERS)

}


# read in 2 seurat objects under consideration
All.integrated_30 <- LoadH5Seurat('outputs/All.integrated_30.h5seurat')
All.integrated_86 <- LoadH5Seurat('outputs/All.integrated_86.h5seurat')


# USE MARKER GENES FROM TABLE
DOT_PLOT <- read_csv('raw_data/dot_plot_markers.csv')

# read in markers and also removes '/' and ' ' characters from the type column
MARKERS <- read_tsv('outputs/marker_genes_long.tsv') %>%
  mutate(type=sub('[/ ]','_',type))

# check the seurat object to see which marker genes are present
# amatch uses aproximate grep to see if typos were responsible for not matching
MARKERS <-
  MARKERS %>%
  mutate(amatches=map(.x=gene_name, .f=~agrep(.x, rownames(All.integrated_30), value = T)),
         matches=map(.x=paste0('^',gene_name, '$'), .f=~grep(.x, rownames(All.integrated_30), value = T)))

# matches genes specified in the 'dot plot' table
MARKERS_dot <-
  DOT_PLOT %>%
  mutate(amatches=map(.x=gene_name, .f=~agrep(.x, rownames(All.integrated_30), value = T)),
         matches=map(.x=paste0('^',gene_name, '$'), .f=~grep(.x, rownames(All.integrated_30), value = T)))

# these requested markers were not matched
UN_matched_markers <-
  MARKERS %>%
  filter(map_lgl(.x=matches, .f=~identical(.x, character(0)) ))


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
  mutate(feature_plot=map(.x=marker_genes, .f=~FeaturePlot(All.integrated_30, features = .x)))


# I used these feature plots to check which marker genes were useful
# I remove those that are not

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

# removes non-informative markers
matched_markers_filt <-
  matched_markers %>%
  filter(!(gene_name %in% c('XBP1', 'PRDM1') & type == 'B_ASC')) %>%
  filter(type != 'Others') %>%
  filter(!(type == 'neutrophil' & gene_name %in% c('MPO', 'FCGR1A'))) %>%
  filter(!(type == 'myeloid' & gene_name %in% c('CD93', 'FLT3'))) %>%
  filter(!(type == 'monocyte' & gene_name %in% c('CCR2')))


matched_markers_filt %>% write_tsv('outputs/matched_markers_filt.tsv')

# pub for some marker genes
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5742132/

number_cells_expressing_markers <-
  tibble(gene_name=names(rowSums(All.integrated_30[matched_markers_filt$gene_name,]@assays$originalexp@counts > 0)),
         num_positive=rowSums(All.integrated_30[matched_markers_filt$gene_name,]@assays$originalexp@counts > 0))

## make a kable table from this
number_cells_expressing_markers %>%
  arrange((num_positive)) %>%
  write_tsv('outputs/num_cells_expressing_markers.tsv')

# THIS IS TCRD
# grep('ENSBTAG00000055197', All_genes$ensembl_gene_id)


### find all positive markers for each cluster at coarsest resolution
feature_plots_30 <- make_feature_plots(SEURAT = All.integrated_30, matched_markers = matched_markers_filt, label = '30')
feature_plots_86 <- make_feature_plots(SEURAT = All.integrated_86, matched_markers = matched_markers_filt, label = '86')

markers_30 <- get_cluster_markers(SEURAT = All.integrated_30, label='30')
markers_86 <- get_cluster_markers(SEURAT = All.integrated_86, label='86')
