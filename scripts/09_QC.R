library(SingleCellExperiment)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(scater)

# read in the single cell experiment object output by scDblFinder
scDbl_out <- read_rds('outputs/scDblFinder_out.rds')
# colnames(scDbl_out)

# all detected genes and some data about them
GENE_IDS <- read_tsv('outputs/gene_ID_mapping.tsv')

# many genes with functions associated with the mitochondira are not on the
# mitochondrial genome, many migrate to the host chromosome.

# detecting 'mitochondrial genes' using mitochondrial key words in protein descriptions
# will over estimate the presence of transcripts from genes encoded on the mitochondrial genome
# mito_genes <- GENE_IDS %>% filter(grepl('mitochon', description))

# I appended 'MT' to mitochondrial genes based on genomic location from ensebml
# when I first created a single-cell experiment object in the script: 08_scDblFinder.R
MITO_GENES <- grep('^MT-', rownames(scDbl_out), value = T)

# rRNA that should probably be pretty low in most cells because of rRNA depletion
rRNAs <- GENE_IDS %>% filter(grepl('_rRNA', name)) %>% pull(ROWNAMES)

ribosome_assoc <- GENE_IDS %>%
  filter(grepl('ribosom', description))%>%
  filter(!grepl('_rRNA', name)) %>%
  pull(ROWNAMES)

# QC_SUBSETS <- list(mitochondria=MITO_GENES, ribosome=ribosome_assoc)
QC_SUBSETS <- list(mitochondria=MITO_GENES, rRNA=rRNAs, ribosome=ribosome_assoc)

# use the scater package to add per cell QC metrics and %s for the subsets above
scDbl_out <- scater::addPerCellQC(scDbl_out, subsets=QC_SUBSETS)

# convert to seurat object
seu <- scDbl_out %>% as.Seurat(data=NULL)

GENE_IDS$seurat_IDs <- rownames(seu)
write_tsv(GENE_IDS, 'outputs/gene_ID_mapping.tsv')


seu@meta.data %>%
  mutate(CpF=log(nCount_originalexp/nFeature_originalexp)) %>%
  pull(CpF) %>%
  hist(xlab='log(nCount/nFeature)',
       main='log(total counts / number of genes) per cell' )

# from https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
# to help id low complexit cells
seu@meta.data$log10GenesPerUMI <- log10(seu@meta.data$nFeature_originalexp) / log10(seu@meta.data$nCount_originalexp)

seu@meta.data %>% 
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_histogram() + 
  geom_vline(xintercept = .8)+
  facet_wrap(~individual + tissue) + 
  ggtitle("Log10 Genes / Log10 UMI", 'proxy for cell complexity, cutoff = 0.8')
ggsave('outputs/figures/genes_per_umi_cutoff.jpeg', width = 6, height = 4, units = 'in', bg = 'white')


# UMIs per cell
seu@meta.data %>%
  ggplot(aes(x=nCount_originalexp)) +
  geom_histogram(bins=50)+ 
  scale_x_log10()+
  facet_wrap(~tissue+individual) + #xlim(0,5000) +
  geom_vline(xintercept = c(1200)) +
  ggtitle('UMIs per cell (log scale)', 'cutoff = 1200')

ggsave('outputs/figures/UMIs_per_cell_cutoff.jpeg', width = 6, height = 4, units = 'in', bg = 'white')

#
# Genes per cell
seu@meta.data %>%
  ggplot(aes(x=nFeature_originalexp)) +
  geom_histogram()+ scale_x_log10()+
  facet_wrap(~tissue+individual) + #xlim(0,2000) +
  geom_vline(xintercept = c(500)) +
  ggtitle('Genes per cell (log scale)', 'cutoff = 500')

ggsave('outputs/figures/genes_per_cell_cutoff.jpeg', width = 6, height = 4, units = 'in', bg = 'white')

#
seu@meta.data %>% 
  ggplot(aes(y=nFeature_originalexp,
             x=nCount_originalexp,
             color=subsets_mitochondria_percent)) + 
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1200) +
  geom_hline(yintercept = 500)+
  facet_wrap(~tissue)  
  
  


# percent mt
seu@meta.data %>%
  ggplot(aes(x=individual, y=subsets_mitochondria_percent, fill=tissue)) + geom_violin()

# cutoff for mitochondira percent
seu@meta.data %>%
  ggplot(aes(x=subsets_mitochondria_percent)) + geom_histogram(bins=50)+
  facet_wrap(~tissue+individual) + xlim(0,25) +
  geom_vline(xintercept = c(10))
ggsave('outputs/figures/mitochondria_percent.jpeg', width = 7, height = 5, units = 'in', bg='white')

## SHOW THIS ONE
seu@meta.data %>%
  ggplot(aes(x=individual, y=subsets_rRNA_percent, fill=tissue)) + geom_violin() +
  ggtitle('rRNA depletion appears to have failed for one sample')

ggsave('outputs/figures/rRNA_percent.jpeg', width = 7, height = 5, units = 'in', bg='white')
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

ggsave('outputs/figures/num_doublets.jpeg', width = 7, height = 5, units = 'in', bg='white')


seu@meta.data %>%
  ggplot(aes(x=individual, y=subsets_ribosome_percent, fill=tissue)) + geom_violin()

# THIS ONE TO SHOW rRNA DEPLETION DIDNT WORK IN ONE MILK SAMPLE

### REMOVE rRNA genes here
non_rRNA_genes <- GENE_IDS %>% filter(!grepl('_rRNA', name)) %>% pull(seurat_IDs)
seu <- subset(seu, features=non_rRNA_genes)

# 10-8-2022 CHANGED PCT MITO FROM 5 TO 10!!!!
# calculate cells to remove
seu@meta.data <-
  seu@meta.data %>%
  mutate(REMOVE=case_when(
    scDblFinder.class == 'doublet'   ~ 'DOUBLET',
    subsets_mitochondria_percent > 10 ~ 'PCT_MIT_ABOVE_10',
    nFeature_originalexp < 500       ~  'GENE_COUNT_BELOW_500',
    nCount_originalexp < 1200        ~  'UMIs_BELOW_1200',
    log10GenesPerUMI < 0.8           ~  'COMPLEXITY_BELOW_0.8',
    # nFeature_originalexp > 3000      ~ 'GENE_COUNT_ABOVE_3000',
    TRUE                             ~ 'KEEP'),
    REMOVE=factor(REMOVE, levels = c('KEEP', 'DOUBLET', 'GENE_COUNT_BELOW_500',
                                     'PCT_MIT_ABOVE_10', 'UMIs_BELOW_1200','COMPLEXITY_BELOW_0.8' )))


seu@meta.data %>%
  group_by(tissue, individual) %>%
  count(REMOVE) %>%
  ggplot(aes(y=REMOVE, x=n, fill=REMOVE)) +
  geom_col(color='black') +
  geom_text(aes(label=n),hjust=0 )+
  facet_wrap(~individual+tissue, ncol=2) +
  xlab('number of cells') +
  theme(legend.position = 'none') +
  xlim(0,8500) +
  ggtitle('Cells removed by different QC metrics')

ggsave('outputs/figures/removed_cells.jpeg', width = 7, height = 5, units = 'in', bg='white')

# remove bad cells
seu_filt <- subset(seu, subset = REMOVE == 'KEEP')

# remove genes with no expression in remaining cells
non_zero_Features <- names(which(!rowSums(seu_filt) == 0))

seu_filt <- subset(seu_filt, features=non_zero_Features)

# remove genes expressed in fewer than 10 cells

counts <- GetAssayData(object = seu_filt, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = seu_filt@meta.data)


# write out seurat object
seu_filt <- SeuratDisk::SaveH5Seurat(filtered_seurat, 'outputs/seurat_QC_done', overwrite = TRUE)






