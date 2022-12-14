library(tidyverse)
library(Seurat)
library(SeuratDisk)


All.integrated <- LoadH5Seurat('outputs/All.integrated_30.h5seurat')

# need to make sure this matches the seurat object read in
FOUND_MARKERS <- read_tsv('outputs/30_markers.tsv')

matched_markers_filt <- read_tsv('outputs/matched_markers_filt.tsv')

FOUND_MARKERS %>%
  # rownames_to_column(var='gene_name') %>%
  group_by(cluster) %>%
  arrange(desc(pct.1))

# FOUND_MARKERS$pct.1 %>% hist()
# FOUND_MARKERS$avg_log2FC %>% hist(breaks=100)

FOUND_MARKERS %>%
  filter(gene %in% matched_markers_filt$gene_name) %>%
  left_join(matched_markers_filt, by=c('gene' = 'gene_name')) %>%
  group_by(cluster, type) %>% tally() %>%
  ggplot(aes(x=factor(cluster, levels = c(0:20)),y=n, fill=type)) +
  geom_col(color='black')

# make dot plots from this?
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
DimPlot(All.integrated, group.by = 'integrated_snn_res.0.5', label = T)

ggsave('outputs/figures/dimplot_30_cluster_labels.jpeg', width = 9, height = 5, units = 'in', bg='white')
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


# NAMED DOTS
# need a named list of feature vectors
NDOTS1 <-
  matched_markers_filt %>%
  filter(!(type == 'myeloid' & gene_name == 'TCF4')) %>%
  filter(type != 'monocyte') %>% #group_by(gene_name) %>% tally() %>% arrange(desc(n))
  group_by(type) %>%
  summarise(data=list(gene_name)) %>%
  ungroup()

NDOT_LIST <- setNames(as.list(NDOTS1$data), NDOTS1$type)


# I also found a spreadsheet titled 'dot plots' 
# may want to look at these genes too.  

# DOTS <- read_csv('raw_data/dot_plot_markers.csv') %>%
#   filter(!(gene_name %in% c('BOLA-DRB1','FTL3'))) %>%
#   pull(gene_name) %>% unique()
#
# DOTS2 <- c(DOTS, DOTS1) %>% unique()


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

ggsave('outputs/figures/initial_dot_plot.jpeg', height=5, width = 9, units = 'in', bg='white')


pdot <-
  DotPlot(All.integrated,
        features = NDOT_LIST,
        group.by = 'integrated_snn_res.0.5',
        cluster.idents = T) +
  theme(axis.text.x = element_text(angle=-45, size=7.5, hjust = 0),
        panel.border = element_rect(fill=NA, color='black'),
        panel.grid.major.y = element_line(color='grey90')) +
  ggtitle('Pre-identified marker genes for cell types')
pdot
ggsave('outputs/figures/grouped_dot_plot.jpeg', height=6, width = 11, units = 'in', bg='white')

# this function takes:
# 1) a dot plot object that was generated with a named list of feature vectors
#   - the names of the list define the group names,
# 2)  the column name that defines the clustering
# returns the cluster IDs with the max scoring group (from the list names)
# scores are (Percent Expressed) * (Average Expression) summed for each gene in the group
cluster_groups_from_dotplot <- function(dot_plot, cluster_column_name){
  result <-
    dot_plot$data %>%
    group_by(feature.groups, id) %>%
    summarise(INDEX=pct.exp * avg.exp, .groups = 'drop') %>%
    group_by(id) %>%
    summarise(dot_plot_group=feature.groups[which.max(INDEX)],
              .groups = 'drop')
  result[[cluster_column_name]] <- result$id
  result <- result[,- which(colnames(result) == 'id')]
  return(result)
}


dot_plot_assignments <-
  cluster_groups_from_dotplot(pdot, cluster_column_name = 'integrated_snn_res.0.5')

# best to let the immunologists do the manual IDs...
cluster_identities <-
  cluster_identities %>%
  # mutate(manual_ID = case_when(
  #   integrated_snn_res.0.5  == 2 ~ paste0('B_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 14 ~ paste0('GD_T_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 14 ~ paste0('GD_T_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 4 ~ paste0('CD4_T_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 11 ~ paste0('CD8_T_',integrated_snn_res.0.5),
  #   TRUE ~ paste0(max_ident,'_',integrated_snn_res.0.5),
  # )) %>%
  # mutate(manual_ID_coarse = case_when(
  #   integrated_snn_res.0.5  == 4 ~ paste0('B_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 18 ~ paste0('GD_T_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 14 ~ paste0('GD_T_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 2 ~ paste0('CD4_T_',integrated_snn_res.0.5),
  #   integrated_snn_res.0.5  == 11 ~ paste0('CD8_T_',integrated_snn_res.0.5),
  #   TRUE ~ max_ident,
  # )) %>%
  left_join(dot_plot_assignments)

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

# DimPlot(All.integrated, group.by = 'integrated_snn_res.0.5')
# DimPlot(All.integrated, group.by = 'max_ident', split.by = 'tissue')
# DimPlot(All.integrated, group.by = 'all_ident', split.by = 'tissue')


# THIS ONE
# DimPlot(All.integrated, group.by = 'manual_ID', split.by = 'tissue',label = T, label.size = 3)
# ggsave('outputs/figures/dim_plot_manual_ID.jpeg', width = 7, height = 5, units = 'in', bg='white')

# THIS ONE
# DimPlot(All.integrated, group.by = 'manual_ID_coarse', split.by = 'tissue',label = T, label.size = 3) +
  # theme(panel.border = element_rect(fill=NA, color='black'))
# ggsave('outputs/figures/dim_plot_manual_ID_coarse.jpeg', width = 7, height = 5, units = 'in', bg='white')

DimPlot(All.integrated, group.by = 'dot_plot_group', split.by = 'tissue',label = TRUE) +
  theme(panel.border = element_rect(fill=NA, color='black')) + 
  ggtitle('Cell types from marker genes')

ggsave('outputs/figures/dim_plot_dot_plot_group.jpeg', width = 7, height = 5, units = 'in', bg='white')


#PBMC single cell dairy cattle
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08562-0


SaveH5Seurat(All.integrated, 'outputs/classified_clusters_30', overwrite=TRUE)

### cell types proportions across samples:

cell_type_counts <- 
  All.integrated@meta.data %>%
  group_by(sample_ID,individual,tissue, dot_plot_group) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(sample_ID,individual,tissue) %>% 
  mutate(percent_abund = (n / sum(n)) * 100) 

cell_type_counts %>%
  ggplot(aes(x=individual, y=percent_abund, fill=dot_plot_group)) + 
  geom_col() +
  facet_wrap(~tissue) + 
  ggtitle('Percent abundance of cell types')

ggsave('outputs/figures/cell_type_relative_abundance.jpeg',
       width=7, height=5, units='in',bg='white')


cell_type_counts %>%
  ggplot(aes(x=tissue, y=percent_abund, fill=dot_plot_group)) + 
  geom_boxplot() 

cluster_cell_type_counts <- 
  All.integrated@meta.data %>%
  group_by(integrated_snn_res.0.5,sample_ID,individual,tissue, dot_plot_group) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(sample_ID,individual,tissue) %>% 
  mutate(percent_abund = (n / sum(n)) * 100) 


cluster_cell_type_counts %>%
  mutate(FACET=paste0(integrated_snn_res.0.5, '_',dot_plot_group)) %>% 
  ggplot(aes(x=tissue, y=percent_abund, fill=individual)) +
  geom_col(position = position_dodge(), color='white') +
  facet_wrap(~FACET, scales = 'free', nrow = 3) + 
  theme_bw() + ylab('percent abundance of sample') + 
  ggtitle('Cluster abundances')

ggsave('outputs/figures/cluster_relative_abundances.jpeg', height = 5, width = 9, units = 'in', bg='white')
# 1634 has relatively high levels of some T cells and B cells in the blood
  # clusters 11, 12, 14, 15, 16 

cluster_cell_type_counts %>% ggplot(aes(x=individual, y=percent_abund, fill=integrated_snn_res.0.5)) + 
  geom_col() + 
  facet_wrap(~tissue)


# cell_type_counts %>%
#   ggplot(aes(x=dot_plot_group, y=percent_abund, fill=tissue)) + 
#   geom_boxplot()+
#   facet_wrap(~dot_plot_group)


#rmarkdown::render("mastitis_scRNAseq_report.Rmd")
