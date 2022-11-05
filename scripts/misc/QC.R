# doublet removal and qc
# generally following the seurat pbmc tutorial
library(SingleCellExperiment)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(scater)


scDbl_out <- read_rds('outputs/scDblFinder_out.rds') 


GENE_IDS <- read_tsv('outputs/gene_ID_mapping.tsv')

# mito_genes <- GENE_IDS %>% filter(grepl('mitochon', description))

MITO_GENES <- grep('^MT-', rownames(scDbl_out), value = T)
rRNAs <- GENE_IDS %>% filter(grepl('_rRNA', name)) %>% pull(ROWNAMES)
ribosome_assoc <- GENE_IDS %>%
  filter(grepl('ribosom', description))%>% 
  filter(!grepl('_rRNA', name)) %>% 
  pull(ROWNAMES)

# QC_SUBSETS <- list(mitochondria=MITO_GENES, ribosome=ribosome_assoc)
QC_SUBSETS <- list(mitochondria=MITO_GENES, rRNA=rRNAs, ribosome=ribosome_assoc)

# remove rRNA genes

###.

scDbl_out <- scater::addPerCellQC(scDbl_out, subsets=QC_SUBSETS)

# REMOVE rRNAs,
# scDbl_out <- scDbl_out[!rownames(scDbl_out) %in% rRNAs,]
##
# convert to seurat object
seu <- scDbl_out %>% as.Seurat(data=NULL)

GENE_IDS$seurat_IDs <- rownames(seu)
write_tsv(GENE_IDS, 'outputs/gene_ID_mapping.tsv')


#
# rownames(seu)
#

# rm(scDbl_out)
# calculate percent mitochondrial transcripts
# seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# SHOULD DO NUM FEATURES DIVIDED BY TOTAL READS IN THE CELL

seu@meta.data %>% 
  mutate(CpF=log(nCount_originalexp/nFeature_originalexp)) %>%
  pull(CpF) %>% 
  hist(xlab='log(nCount/nFeature)', 
       main='log(total counts / number of genes) per cell' )

# percent mt
seu@meta.data %>% 
  ggplot(aes(x=individual, y=subsets_mitochondria_percent, fill=tissue)) + geom_violin()

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
#
# CHECK SOUPX #

seu@meta.data %>% 
  ggplot(aes(x=individual, y=subsets_ribosome_percent, fill=tissue)) + geom_violin()

# THIS ONE TO SHOW rRNA DEPLETION DIDNT WORK IN ONE MILK SAMPLE
seu@meta.data %>% 
  ggplot(aes(x=individual, y=subsets_rRNA_percent, fill=tissue)) + geom_violin()

### REMOVE rRNA genes here
non_rRNA_genes <- GENE_IDS %>% filter(!grepl('_rRNA', name)) %>% pull(seurat_IDs)
seu <- subset(seu, features=non_rRNA_genes)


# calculate cells to remove
seu@meta.data <- 
  seu@meta.data %>%
  mutate(REMOVE=case_when(
    scDblFinder.class == 'doublet'   ~ 'DOUBLET',     
    subsets_mitochondria_percent > 5 ~ 'PCT_MIT_ABOVE_5',
    nFeature_originalexp < 150       ~  'GENE_COUNT_BELOW_150',
    nFeature_originalexp > 3000      ~ 'GENE_COUNT_ABOVE_3000', 
    TRUE                             ~ 'KEEP'), 
    REMOVE=factor(REMOVE, levels = c('KEEP', 'DOUBLET', 'GENE_COUNT_BELOW_150',
                                     'GENE_COUNT_ABOVE_3000', 'PCT_MIT_ABOVE_5')))


seu@meta.data %>% 
  group_by(tissue, individual) %>% 
  count(REMOVE) %>%
  ggplot(aes(y=REMOVE, x=n, fill=REMOVE)) + 
  geom_col(color='black') + 
  geom_text(aes(label=n),hjust=0 )+
  facet_wrap(~individual+tissue, ncol=2) +
  xlab('number of cells') + 
  theme(legend.position = 'none') + 
  xlim(0,6500) + 
  ggtitle('Cells removed by differnt QC metrics')

# remove bad cells
seu_filt <- subset(seu, subset = REMOVE == 'KEEP')

# remove genes with no expression in remaining cells
non_zero_Features <- names(which(!rowSums(seu_filt) == 0))

seu_filt <- subset(seu_filt, features=non_zero_Features)


seu_filt <- write_rds(seu_filt, 'outputs/seurat_QC_done.rds')
### SPLIT SCRIPT HERE






