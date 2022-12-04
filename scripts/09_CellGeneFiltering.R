#20221129 - Jayne Wiarda
# Cell & gene filtering

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
#library(ggridges)

seu <- LoadH5Seurat('outputs/pre_QC.h5seurat')
seu[['RNA']] <- seu[['originalexp']] 
DefaultAssay(seu) <- 'RNA'
seu[['originalexp']] <- NULL

seu@meta.data %>%
  ggplot(aes(x=nCount_RNA)) +
  geom_histogram(bins=100)+ 
  scale_x_log10()+
  facet_wrap(~tissue+individual) + #xlim(0,5000) +
  geom_vline(xintercept = c(500)) +
  ggtitle('UMIs per cell (log scale)', 'cutoff = 500')

seu@meta.data %>%
  ggplot(aes(x=nFeature_RNA)) +
  geom_histogram(bins = 100)+ 
  scale_x_log10()+
  facet_wrap(~tissue+individual) + #xlim(0,2000) +
  geom_vline(xintercept = c(250)) +
  ggtitle('Genes per cell (log scale)', 'cutoff = 250')

seu@meta.data %>%
  ggplot(aes(x=subsets_mitochondria_percent)) + 
  geom_histogram(bins=100)+
  facet_wrap(~tissue+individual) + xlim(0,25) +
  geom_vline(xintercept = c(12.5))+
  ggtitle('Percent mitochondrial reads per cell', 'cutoff = 12.5%')

seu@meta.data %>% 
  ggplot(aes(y=nFeature_RNA,
             x=nCount_RNA,
             color=subsets_mitochondria_percent)) + 
  geom_point() +
  scale_colour_gradient(low = "grey", high = "red4") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)+
  facet_wrap(~tissue+individual) 

seu@meta.data %>% 
  ggplot(aes(y=nFeature_RNA,
             x=nCount_RNA,
             color=subsets_mitochondria_percent)) + 
  geom_point(color = dplyr::case_when(seu@meta.data$scDblFinder.class == 'singlet' ~ "grey70",
                                      seu@meta.data$scDblFinder.class == 'doublet' ~ "red4")) +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)+
  facet_wrap(~tissue+individual) 

table(seu$scDblFinder.class, seu$sample_ID)

seu@meta.data %>% 
  ggplot(aes(y=nFeature_RNA,
             x=nCount_RNA,
             color=subsets_mitochondria_percent)) + 
  geom_point(color = dplyr::case_when(seu@meta.data$subsets_mitochondria_percent < 12.5 & 
                                        seu@meta.data$nCount_RNA > 500 & 
                                        seu@meta.data$nFeature_RNA > 250 & 
                                        seu@meta.data$scDblFinder.class == 'singlet' ~ "grey70",
                                      TRUE ~ "red4")) +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)+
  facet_wrap(~tissue+individual) 

seu@meta.data %>%
  ggplot(aes(x=individual, y=subsets_rRNA_percent, fill=tissue)) + geom_violin() +
  ggtitle('rRNA depletion appears to have failed for one sample')

# num features
seu@meta.data %>%
  ggplot(aes(x=individual, y=nFeature_originalexp, fill=tissue)) +
  geom_violin()

# num features by doublet status
seu@meta.data %>%
  ggplot(aes(x=scDblFinder.class, y=nFeature_originalexp, fill=tissue)) +
  geom_violin()

# this one for singlets vs doublets
seu@meta.data %>%
  group_by(sample_ID, tissue, individual, scDblFinder.class) %>%
  tally() %>%
  ggplot(aes(x=scDblFinder.class, y=n, fill=tissue)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~individual)

Idents(seu) <- seu$sample_ID
RidgePlot(seu, 'log10GenesPerUMI', same.y.lims = TRUE)


### REMOVE rRNA genes here
non_rRNA_genes <- GENE_IDS %>% filter(!grepl('_rRNA', name)) %>% pull(seurat_IDs)
seu <- subset(seu, features=non_rRNA_genes)

seu@meta.data <-
  seu@meta.data %>%
  mutate(REMOVE=case_when(
    scDblFinder.class == 'doublet'   ~ 'DOUBLET',
    subsets_mitochondria_percent > 12.5 ~ 'PCT_MIT_ABOVE_12.5',
    nFeature_originalexp < 250       ~  'GENE_COUNT_BELOW_250',
    nCount_originalexp < 500        ~  'UMIs_BELOW_500',
    TRUE                             ~ 'KEEP'),
    REMOVE=factor(REMOVE, levels = c('KEEP', 'DOUBLET', 'GENE_COUNT_BELOW_250',
                                     'PCT_MIT_ABOVE_12.5', 'UMIs_BELOW_500')))

seu@meta.data %>%
  group_by(tissue, individual) %>%
  count(REMOVE) %>%
  ggplot(aes(y=REMOVE, x=n, fill=REMOVE)) +
  geom_col(color='black') +
  geom_text(aes(label=n),hjust=0 )+
  facet_wrap(~individual+tissue, ncol=2) +
  xlab('number of cells') +
  theme(legend.position = 'none') +
  xlim(0,10000) +
  ggtitle('Cells removed by different QC metrics')

# remove bad cells
seu_filt <- subset(seu, subset = REMOVE == 'KEEP')

# remove genes with no expression in remaining cells
non_zero_Features <- names(which(!rowSums(seu_filt) == 0))

seu_filt <- subset(seu_filt, features=non_zero_Features)

SaveH5Seurat(seu_filt, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20221129_JEW_FilteredSeurat.h5seurat', overwrite = TRUE)

