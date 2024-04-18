library(tidyverse)
# library(Seurat)
# library(SeuratDisk)
library(glue)
library(future)


if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

options(future.globals.maxSize = 64000 * 1024^2)

#######################


#seu_filt <- SeuratDisk::LoadH5Seurat('outputs/granulocytes_filtered')

ellie_genes <- read_tsv('raw_data/NeutrophilGenesofInterest_EP_24May2023.tsv')



library(biomaRt)
library(tidyverse)


# scdblout <- read_rds('outputs/scDblFinder_out.rds')
# 
# All_genes <- read_tsv('outputs/gene_ID_mapping.tsv')

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

alternative_names <- 
  c(ITGAX = 'ENSBTAG00000054461',
  CTSG =  'ENSBTAG00000040134',
  CD24A =  'ENSBTAG00000039046',
  IFI16 = 'ENSBTAG00000011511',
  IFITM2 = 'ENSBTAG00000019017')



ellie_not_found <- 
  ellie_not_found %>%
  mutate(gene_name = alternative_names[clean_gene]) %>% 
  filter(!is.na(gene_name))

ellie_found <- ellie_found %>% mutate(gene_name=clean_gene)


ellie_final <- bind_rows(ellie_found, ellie_not_found) %>% 
  transmute(Name, `Gene Name`, `Source of interest`, Priority, `Notes:`, gene_name)
ellie_final %>% write_tsv('outputs/ellie_genes.tsv')

# 
# 
# MOUSE <- read_tsv('outputs/mouse_homologs.tsv')
# HUMAN <- read_tsv('outputs/human_homologs.tsv')
# 
# 'Hepcarcin' %in% MOUSE$mmusculus_homolog_associated_gene_name
# 'ENSG00000245532' %in% HUMAN$hsapiens_homolog_ensembl_gene
# 
# 'ENSBTAG00000019017' %in% All_genes$ensembl_gene_id
# 
# 
# All_genes$name
# All_genes
# 
# synonyms <- read_tsv('outputs/gene_synonyms.tsv')
# 
# # USE MARKER GENES FROM TABLE
# DOT_PLOT <- read_csv('raw_data/dot_plot_markers.csv')
# 
# # read in markers and also removes '/' and ' ' characters from the type column
# MARKERS <- read_tsv('outputs/marker_genes_long.tsv') %>%
#   mutate(type=sub('[/ ]','_',type))
# 
# # check the seurat object to see which marker genes are present
# # amatch uses aproximate grep to see if typos were responsible for not matching
# MARKERS <-
#   MARKERS %>%
#   mutate(amatches=map(.x=gene_name, .f=~agrep(.x, rownames(All.integrated_30), value = T)),
#          matches=map(.x=paste0('^',gene_name, '$'), .f=~grep(.x, rownames(All.integrated_30), value = T)))
# 
# # matches genes specified in the 'dot plot' table
# MARKERS_dot <-
#   DOT_PLOT %>%
#   mutate(amatches=map(.x=gene_name, .f=~agrep(.x, rownames(All.integrated_30), value = T)),
#          matches=map(.x=paste0('^',gene_name, '$'), .f=~grep(.x, rownames(All.integrated_30), value = T)))
# 
# # these requested markers were not matched
# UN_matched_markers <-
#   MARKERS %>%
#   filter(map_lgl(.x=matches, .f=~identical(.x, character(0)) ))
# 
# 
# # PAX5 == ENSBTAG00000012498
# # cd16 synonym
# # itgam = cd18
# # CD31 = PECAM1
# # CD64 = FCGR1A
# # glut1 = SLC2A1
# # SELL, CD62L, LAM1, LECAM1, LEU8, LNHR, LSEL, LYAM1, PLNHR, TQ1, selectin L
# 
# matched_markers <-
#   MARKERS %>%
#   filter(!map_lgl(.x=matches, .f=~identical(.x, character(0)) )) %>%
#   select(-amatches, -matches)
