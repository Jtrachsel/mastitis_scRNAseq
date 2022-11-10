#probably should do this after QC...

# map ensembl geneids to gene names, GO terms, homologs etc.

library(biomaRt)
library(tidyverse)

### collect all detected genes from cellranger output


cellranger_dirs <- 
  list.dirs(recursive = FALSE, path = 'cellranger_out/') #%>%
  # grep('-[1-4]$', ., value = TRUE)

# get a dataframe of all genes detected
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
# sometimes ensembl servers are busy and non-responsive.
# should maybe have a way of re-trying if this causes script to fail
ensembl <- useEnsembl(biomart = 'genes', dataset = 'btaurus_gene_ensembl')


DATASETS <- listDatasets(mart = useMart('ensembl'))
DATASETS %>% filter(grepl('btaurus',dataset)) %>% pull(dataset)

mart <- useDataset("btaurus_gene_ensembl", useMart("ensembl"))
FILTERS <- listFilters(mart)

# names for orthologs
ortholog_filters <- 
  bind_rows(
  FILTERS %>% filter(grepl('Orthologous Human', description)),
  FILTERS %>% filter(grepl('Orthologous Pig - Duroc Genes', description)),
  FILTERS %>% filter(grepl('Orthologous Mouse Genes', description))
  ) %>% 
  pull(name)

ATTRIBUTES <- listAttributes(mart)

human_homologs <- 
  ATTRIBUTES %>% 
  filter(page == 'homologs') %>% 
  filter(grepl('Human', description)) %>% 
  pull(name)

mouse_homologs <- 
  ATTRIBUTES %>%
  filter(page == 'homologs') %>% 
  filter(grepl('mmusculus', name)) %>% 
  pull(name)

pig_homologs <- 
  ATTRIBUTES %>%
  filter(page == 'homologs') %>% 
  filter(grepl('Pig - Duroc', description)) %>% 
  pull(name)


# ATTRIBUTES %>% filter(grepl('Synon', description))

G_list <- getBM(filters= c("ensembl_gene_id"),
                attributes= c("ensembl_gene_id",
                              'chromosome_name',
                              "external_gene_name",
                              'description'),
                values=All_genes$ID,
                mart= mart)
# G_list %>% pull(chromosome_name) %>% unique()
# G_list %>% filter(chromosome_name == 'MT')

synonyms <- getBM(filters= c("ensembl_gene_id"),
                 attributes= c("ensembl_gene_id",
                               "external_gene_name",
                               'external_synonym',
                               'description'),
                 values=All_genes$ID,
                 mart= mart)

# synonyms

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



### homologs 

human_homolog_df <- 
  getBM(filters= c("ensembl_gene_id"),
                  attributes= c("ensembl_gene_id",
                              human_homologs),
                  values=All_genes$ID,
                  mart= mart) %>% 
  write_tsv('outputs/human_homologs.tsv')

mouse_homolog_df <-
  getBM(filters= c("ensembl_gene_id"),
                  attributes= c("ensembl_gene_id",
                                mouse_homologs),
                  values=All_genes$ID,
                  mart= mart)%>% 
  write_tsv('outputs/mouse_homologs.tsv')

pig_homologs_df <- 
  getBM(filters= c("ensembl_gene_id"),
                  attributes= c("ensembl_gene_id",
                                pig_homologs),
                  values=All_genes$ID,
                  mart= mart) %>% 
  write_tsv('outputs/pig_homologs.tsv')

synonyms %>% write_tsv('outputs/gene_synonyms.tsv')


All_genes <- 
  All_genes %>%
  transmute(ensembl_gene_id=ID, name) %>%
  left_join(G_list)


write_tsv(All_genes, 'outputs/gene_ID_mapping.tsv')
write_tsv(GO_terms, 'outputs/gene_ID_GO_terms.tsv')
write_tsv(reactome, 'outputs/gene_ID_reactome.tsv')

#### marker gene checks

All_genes <- read_tsv('outputs/gene_ID_mapping.tsv')

# All_genes %>% count(ensembl_gene_id) %>% filter(n != 1)

# grep('SIRPA', All_genes$external_gene_name, value = T)

# dot plot genes

dot_plot <- 
  read_csv('raw_data/dot_plot_markers.csv')

not_found_dot <- dot_plot %>% filter(!gene_name %in% All_genes$name) %>% pull(gene_name)

MARKER_GENES <- 
  read_tsv('raw_data/marker_genes_updated.tsv') %>% 
  pivot_longer(cols = everything(), names_to = 'type', values_to = 'gene_name') %>% 
  arrange(type) %>% filter(!is.na(gene_name))


# all external gene names I got from ensmbl were identical to the ones cellranger
# pulled from the gtf file.  No new info gained there.  But mapped GO terms so thats nice
# MARKER_GENES %>% filter(gene_name %in% All_genes$external_gene_name) %>% unique()

# MARKER_GENES %>% filter(!grepl('T|B', type))

MARKER_GENES %>% unique() %>% write_tsv('outputs/marker_genes_long.tsv')

# All_genes %>% filter(ensembl_gene_id == 'ENSBTAG00000055197')

NOT_FOUND_MARKER_GENES <- 
  MARKER_GENES %>%
  filter(!(gene_name %in% All_genes$name)) %>% 
  unique()


