library(tidyverse)
library(Seurat)
library(SeuratDisk)

All.integrated <- LoadH5Seurat('outputs/All.integrated_30.h5seurat')

# need to make sure this matches the seurat object read in 
FOUND_MARKERS <- read_tsv('outputs/FOUND_MARKERS.tsv')

matched_markers_filt <- read_tsv('outputs/matched_markers_filt.tsv')

FOUND_MARKERS %>%
  rownames_to_column(var='gene_name') %>% 
  group_by(cluster) %>% 
  arrange(desc(pct.1))

FOUND_MARKERS$pct.1 %>% hist()
FOUND_MARKERS$avg_log2FC %>% hist(breaks=100)

FOUND_MARKERS %>%
  filter(gene %in% matched_markers_filt$gene_name) %>% 
  left_join(matched_markers_filt, by=c('gene' = 'gene_name')) %>% 
  group_by(cluster, type) %>% tally() %>%
  ggplot(aes(x=factor(cluster, levels = c(0:20)),y=n, fill=type)) + 
  geom_col(color='black')


top5_markers <- 
  FOUND_MARKERS %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n=5) %>% ungroup()


cluster_identities <- 
  FOUND_MARKERS %>%
  filter(gene %in% matched_markers_filt$gene_name) %>%
  left_join(matched_markers_filt, by=c('gene' = 'gene_name')) %>% 
  filter(type != 'monocyte') %>% 
  group_by(cluster, type) %>% tally() %>% 
  arrange(desc(n)) %>% 
  summarise(max_ident=type[which.max(n)], 
            all_ident=paste0(type, collapse = '_')) %>% 
  mutate(integrated_snn_res.0.5=factor(cluster)) %>% 
  select(-cluster)

# All.integrated@meta.data$integrated_snn_res.0.5
# cluster_identities$cluster


# SeuratObject::AddMetaData()
# this works
DimPlot(All.integrated, group.by = 'integrated_snn_res.0.5')


feature_plots_filt <- 
  matched_markers_filt %>%
  group_by(type) %>% 
  summarise(marker_genes=list(gene_name)) %>% 
  mutate(feature_plot=map(.x=marker_genes, .f=~FeaturePlot(All.integrated, features = .x)), 
         split_features=map(.x=marker_genes, .f=~FeaturePlot(All.integrated, features = .x, split.by = 'tissue')))

feature_plots_filt$feature_plot[[1]]
feature_plots_filt$feature_plot[[5]]

feature_plots_filt$split_features[[1]]



DOTS1 <- matched_markers_filt$gene_name %>% unique()

DOTS <- read_csv('raw_data/dot_plot_markers.csv') %>%
  filter(!(gene_name %in% c('BOLA-DRB1','FTL3'))) %>% 
  pull(gene_name) %>% unique()

DOTS2 <- c(DOTS, DOTS1) %>% unique()


# DotPlot(All.integrated,
#         features = DOTS1,
#         group.by = 'integrated_snn_res.0.5', 
#         cluster.idents = T) + 
#   theme(axis.text.x = element_text(angle=-45, size=7.5, hjust = 0))

DotPlot(All.integrated,
        features = DOTS2,
        group.by = 'integrated_snn_res.0.5', 
        cluster.idents = T) + 
  theme(axis.text.x = element_text(angle=-45, size=7.5, hjust = 0))

ggsave('outputs/figures/initial_dot_plot.jpeg', height=5, width = 9, units = 'in')

cluster_identities <- 
  cluster_identities %>% 
  mutate(manual_ID = case_when(
    integrated_snn_res.0.5  == 2 ~ paste0('B_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 20 ~ paste0('GD_T_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 15 ~ paste0('GD_T_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 1 ~ paste0('CD4_T_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 11 ~ paste0('CD8_T_',integrated_snn_res.0.5), 
    TRUE ~ paste0(max_ident,'_',integrated_snn_res.0.5), 
  )) %>% 
  mutate(manual_ID_coarse = case_when(
    integrated_snn_res.0.5  == 2 ~ paste0('B_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 20 ~ paste0('GD_T_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 15 ~ paste0('GD_T_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 1 ~ paste0('CD4_T_',integrated_snn_res.0.5), 
    integrated_snn_res.0.5  == 11 ~ paste0('CD8_T_',integrated_snn_res.0.5), 
    TRUE ~ max_ident, 
  ))

# 2 == B cells

# 20 == gamma delta T cells
# 15 == gamma delta T cells?
# 1 == CD4 T cells
# 11 == CD8 T cells

# rest are myeloid?


All.integrated@meta.data <- 
  All.integrated@meta.data %>% 
  rownames_to_column(var = 'CELL') %>% 
  left_join(cluster_identities) %>% 
  column_to_rownames(var = 'CELL')

DimPlot(All.integrated, group.by = 'integrated_snn_res.0.5')
DimPlot(All.integrated, group.by = 'max_ident', split.by = 'tissue')
DimPlot(All.integrated, group.by = 'all_ident', split.by = 'tissue')


# THIS ONE
DimPlot(All.integrated, group.by = 'manual_ID', split.by = 'tissue',label = T)
ggsave('outputs/figures/dim_plot_manual_ID.jpeg', width = 7, height = 5, units = 'in')

# THIS ONE
DimPlot(All.integrated, group.by = 'manual_ID_coarse', split.by = 'tissue',label = T)
ggsave('outputs/figures/dim_plot_manual_ID_coarse.jpeg', width = 7, height = 5, units = 'in')

#PBMC single cell dairy cattle
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08562-0


UMAP_COORDS <- 
  All.integrated@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'CELL')

All.integrated@meta.data %>% 
  rownames_to_column(var = 'CELL') %>% 
  left_join(UMAP_COORDS) %>% 
  filter(UMAP_1 >4 & UMAP_2 < -5) %>% 
  group_by(all_ident, manual_ID) %>% 
  tally()

