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
granulocytes@assays$SCT <- NULL
granulocytes@assays$integrated <- NULL
# remove genes with no expression in remaining cells
non_zero_Features <- names(which(!rowSums(granulocytes) == 0))

granulocytes <- subset(granulocytes, features=non_zero_Features)

# remove genes expressed in fewer than 10 cells

counts <- GetAssayData(object = granulocytes, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = granulocytes@meta.data)


# write out seurat object
seu_filt <- SeuratDisk::SaveH5Seurat(filtered_seurat, 'outputs/granulocytes_filtered', overwrite = TRUE)



#######################



ellie_genes <- read_tsv('raw_data/NeutrophilGenesofInterest_EP_24May2023.tsv')



library(biomaRt)
library(tidyverse)


# scdblout <- read_rds('outputs/scDblFinder_out.rds')
# 
All_genes <- read_tsv('outputs/gene_ID_mapping.tsv')

All_genes <- read_tsv('outputs/05_gene_ID_mapping.tsv')
# All_genes <- cbind(All_genes,ROWNAMES=rownames(scdblout))

'outputs/gene_synonyms.tsv'
'outputs/human_homologs.tsv'
'outputs/mouse_homologs.tsv'
ellie_genes <- 
  ellie_genes %>% 
  mutate(clean_gene=str_to_upper(`Gene Name`))

ellie_found <- ellie_genes %>% filter(clean_gene %in% All_genes$name)
ellie_not_found <- ellie_genes %>% filter(!(clean_gene %in% All_genes$name))


# I manually looked up the genes that werent matched to the detected genes

# MALAT1 not found
# NEAT1 not found

# ENSBTAG00000054461 ITGAX FOUND
# ENSBTAG00000040134 CTSG FOUND
# ENSBTAG00000039046 CD24 FOUND
# ENSBTAG00000054230 IFR7 FOUND
# ENSBTAG00000011511 IFI16 FOUND
# ENSBTAG00000019017 IFITM2 FOUND


MOUSE <- read_tsv('outputs/mouse_homologs.tsv')
HUMAN <- read_tsv('outputs/human_homologs.tsv')

'Hepcarcin' %in% MOUSE$mmusculus_homolog_associated_gene_name
'ENSG00000245532' %in% HUMAN$hsapiens_homolog_ensembl_gene

'ENSBTAG00000019017' %in% All_genes$ensembl_gene_id


All_genes$name
All_genes

synonyms <- read_tsv('outputs/gene_synonyms.tsv')

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
