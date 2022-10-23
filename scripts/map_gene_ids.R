# map ensembl geneids to gene names, GO terms, etc.

library(biomaRt)
library(tidyverse)

### collect all detected genes from cellranger output


cellranger_dirs <- 
  list.dirs(recursive = FALSE, path = 'cellranger_out/') #%>%
  # grep('-[1-4]$', ., value = TRUE)


All_genes <- 
  tibble(
    DIR=cellranger_dirs, 
    genes_path=paste0(DIR,'/outs/filtered_feature_bc_matrix/features.tsv.gz'),
    GENES=map(.x=genes_path, .f=~read_tsv(gzfile(.x), col_names = c('ID', 'name', 'type')))
  ) %>% 
  unnest(GENES) %>% 
  select(ID, name) %>%
  unique()

# 26006 total transcripts detected across all samples

####

ensembl <- useEnsembl(biomart = 'genes', dataset = 'btaurus_gene_ensembl')


DATASETS <- listDatasets(mart = useMart('ensembl'))
DATASETS %>% filter(grepl('btaurus',dataset)) %>% pull(dataset)

mart <- useDataset("btaurus_gene_ensembl", useMart("ensembl"))
FILTERS <- listFilters(mart)
ATTRIBUTES <- listAttributes(mart)
grep('name',ATTRIBUTES$description, value = TRUE)


# entrezgene_accession  NCBI gene (formerly Entrezgene) accession(s) [e.g. 20ALPHA-HSD]

# 214 Human gene name
# hgnc_symbol
# HGNC symbol
# 50 external_gene_name Gene Name(s) [e.g. 5S_rRNA]
# 61 hgnc_symbol HGNC symbol(s) [e.g. ACAP2]
G_list <- getBM(filters= c("ensembl_gene_id"),
                attributes= c("ensembl_gene_id",'chromosome_name', "external_gene_name",'description'),
                values=All_genes$ID,mart= mart)

G_list %>% pull(chromosome_name) %>% unique()
G_list %>% filter(chromosome_name == 'MT')

# G_list %>% write_tsv('outputs/ensembl_ID_mapping.tsv')


GO_terms <- getBM(filters= c("ensembl_gene_id"),
                attributes= c("ensembl_gene_id",'go_id'),
                values=All_genes$ID,mart= mart)

GO_terms <- 
  GO_terms %>% group_by(ensembl_gene_id) %>% 
  summarise(GO_ID=paste(go_id, collapse = ';'))

reactome <- getBM(filters= c("ensembl_gene_id"),
                  attributes= c("ensembl_gene_id",'reactome_gene'),
                  values=All_genes$ID,mart= mart)


All_genes <- All_genes %>%transmute(ensembl_gene_id=ID, name) %>%  left_join(G_list)


write_tsv(All_genes, 'outputs/gene_ID_mapping.tsv')
write_tsv(GO_terms, 'outputs/gene_ID_GO_terms.tsv')
write_tsv(reactome, 'outputs/gene_ID_reactome.tsv')



#### marker gene checks
# G_list <- G_list %>% as_tibble()

All_genes <- read_tsv('outputs/gene_ID_mapping.tsv')

grep('SIRPA', All_genes$external_gene_name, value = T)

MARKER_GENES <- read_csv('raw_data/marker_genes_wide.csv') %>% 
  pivot_longer(cols = everything(), names_to = 'type', values_to = 'gene_name') %>% 
  arrange(type) %>% filter(!is.na(gene_name))


# all external gene names I got from ensmbl were identical to the ones cellranger
# pulled from the gtf file.  No new info gained there.  But mapped GO terms so thats nice
MARKER_GENES %>% filter(gene_name %in% All_genes$external_gene_name) %>% unique()

MARKER_GENES %>% filter(!grepl('T|B', type))

MARKER_GENES %>% unique() %>% write_tsv('outputs/marker_genes_long.tsv')

All_genes %>% filter(ensembl_gene_id == 'ENSBTAG00000055197')

NOT_FOUND_MARKER_GENES <- 
  MARKER_GENES %>%
  filter(!(gene_name %in% All_genes$name)) %>% 
  unique()

NOT_FOUND_MARKER_GENES

All_genes %>% filter(grepl('BOLA', external_gene_name))
All_genes %>% filter(grepl('FLT', external_gene_name))

All_genes %>% filter(grepl('FTL', external_gene_name))




grep('CD11', All_genes$external_gene_name, value = T)

# CD11a  == ITGAL
grep('ITGAL', All_genes$external_gene_name, value = T)

grep('ITGAM', All_genes$external_gene_name, value = T)

grep('ITGAX', All_genes$external_gene_name, value = T)
grep('SLEB6', All_genes$external_gene_name, value = T)

# ITGAX, CD11C, SLEB6, integrin subunit alpha X
# CD64
# FCGR1B
# FCGR1C
grep('FCGR1C', All_genes$external_gene_name, value = T)

grep('TCRD', All_genes$external_gene_name, value = T)
# THIS IS TCRD or TRDC
grep('ENSBTAG00000055197', All_genes$ensembl_gene_id)


# XBP1 and XBP2 may be the same thing?
# The XBP2 protein is a reported synonym for the human gene XBP1, 
# NOT_FOUND_MARKER_GENES %>% 
#   mutate(SYNON)







All_genes %>% filter(external_gene_name != name) %>% filter(external_gene_name != '')
