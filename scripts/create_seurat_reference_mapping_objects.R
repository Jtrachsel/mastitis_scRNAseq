# install.packages("anndata")
library(anndata)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(biomaRt)


gtf_parse <- function(path){
  #browser()
  gtf <- read_delim(path,
                    delim = '\t',
                    col_names = c("seqid", "source", "type", "start", "end", "score", "strand","phase","attributes"),
                    comment = '#', progress = FALSE, col_types = c('cccddcccc')) %>%
    filter(type == 'gene') %>%
    tidyr::extract(attributes,
                   into = c('gene_id', 'gene_type', 'gene_name'),
                   regex ='gene_id "(.*)"; gene_type "(.*)"; gene_name "(.*)"; level',
                   remove = FALSE) 
  
  
  
  
  return(gtf)
}

# function to help map IDs to ensembl
# iteratively goes through a vector of attributes and maps only those that didnt 
# map in the previous iteration

map_IDs_to_ensembl <- 
  function(ID_vector, attributes_vector, mart){
    
    # filter the attributes to only those that can be used as a filter
    FILTERS <- listFilters(mart) %>% as_tibble()
    attributes_vector <- attributes_vector[attributes_vector %in% FILTERS$name]
    print('these attributes not used')
    print(attributes_vector[!(attributes_vector %in% FILTERS$name)])
    
    # browser()
    ID_tibble <- tibble(ID=ID_vector)
    
    mapped_tib <- tibble()
    
    for (attribute in attributes_vector){
      # exception for the 1st iteration with no mapping results yet
      if (nrow(mapped_tib) == 0){
        unmapped_tibble <- ID_tibble
        
      } else {
        
        unmapped_tibble <- ID_tibble %>% filter(!(ID %in% mapped_tib$ID))
      }
      print(paste0('mapping ',nrow(unmapped_tibble),' IDs to ', attribute))
      ensembl_mapping <- 
        getBM(values=unmapped_tibble$ID,
              filters = attribute,
              attributes= c("ensembl_gene_id",
                            attribute),
              mart= mart, uniqueRows = TRUE) %>%
        as_tibble()
      
      if (nrow(ensembl_mapping) > 0){
        ensembl_mapping[['ID']] <- ensembl_mapping[[attribute]]
        ensembl_mapping <- ensembl_mapping[,-which(colnames(ensembl_mapping) == attribute)]
        ensembl_mapping[['mapping_type']] <- attribute
        
        mapped_tib <- bind_rows(ensembl_mapping, mapped_tib)
        
        
      } 
      
      
      # wrap for loop iteration here
      
      
    }
    # finalize return here
    
    
    unmapped_tibble <- ID_tibble %>% filter(!(ID %in% mapped_tib$ID))
    result <- bind_rows(mapped_tib, unmapped_tibble)
    
    percent_mapped <- 
      
      return(result)
    
    
  }

# public AWS S3 bucket (https://registry.opendata.aws/tabula-sapiens/)
# tabula used ENCODE GRCh38
# nyquist used GRCh37
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz

#### for reference mapping ####

ENCODE30_GTF <- gtf_parse('gencode.v30.annotation.gtf')

# maps our gene_IDs to ensembl ID
our_genes <- 
  read_tsv('outputs/gene_ID_mapping.tsv') %>% 
  transmute(BOVINE_ID=ensembl_gene_id, our_id=seurat_IDs)


# human homologs to only our detected ensembl_gene_ids
human_homologs <- read_tsv('outputs/human_homologs.tsv') %>% 
  filter(hsapiens_homolog_orthology_type == 'ortholog_one2one')

BOVINIZER <- human_homologs %>% transmute(BOVINE_ID=ensembl_gene_id, ensembl_gene_id=hsapiens_homolog_ensembl_gene)


#######


# Tabula sapiens whole blood data

# some growing pains in the communtiy regarding anndata to seurat formats
# https://github.com/mojaveazure/seurat-disk/issues/109

tabula_blood_anndata <- anndata::read_h5ad('reference_mapping_data/tabula_10X_blood.h5ad')

tabula_blood_counts <- Matrix::t(tabula_blood_anndata$X)

# remove 0 count genes
tabula_blood_counts <- tabula_blood_counts[rowSums(tabula_blood_counts) != 0,]
## RENAME HERE
# listMarts()
# listEnsembl()
# useEnsembl(GRCh = 37, 'genes', dataset=)

# mart=useEnsembl('regulation')
# listDatasets(mart)

# JOIN TO DATA FROM GTF
ENCODE30_mapping <-
  ENCODE30_GTF %>%
  dplyr::transmute(gene_name, ensembl_gene_id=sub('([A-Z0-9]+)\\.[0-9]+','\\1',gene_id))

tabula_blood_mapping <- 
  tibble(gene_name=rownames(tabula_blood_counts)) %>% 
  left_join(ENCODE30_mapping)


bovinized_tabula_blood <- 
  tabula_blood_mapping %>% 
  left_join(BOVINIZER) %>%
  filter(!is.na(BOVINE_ID)) %>% 
  left_join(our_genes)



tabula_swapper <- bovinized_tabula_blood$our_id
names(tabula_swapper) <- bovinized_tabula_blood$gene_name

# subset to only genes with bovine homologs
tabula_blood_counts <- tabula_blood_counts[rownames(tabula_blood_counts) %in% bovinized_tabula_blood$gene_name,]


# swap the gene names to the cattle gene names
rownames(tabula_blood_counts) <- tabula_swapper[rownames(tabula_blood_counts)]

# remove the names from the rownames
names(rownames(tabula_blood_counts)) <- NULL
rownames(tabula_blood_counts)

# bovinized tabula whole blood
# no need for biomaRt because we used the gtf from encode
# tabula_blood_counts


tabula_blood <- CreateSeuratObject(counts = tabula_blood_counts, meta.data = tabula_blood_anndata$obs)

non_zero_genes <- rowSums(tabula_blood) > 0
tabula_blood_filt <- tabula_blood[non_zero_genes,]

SaveH5Seurat(tabula_blood_filt, 'reference_mapping_data/tabula_sapiens_blood', overwrite=TRUE)


### MIT human breast milk study 

# uses weird genenames
# some are gene symbols and some are embl ids?
# info on embl accession in this article
# https://academic.oup.com/nar/article/28/1/19/2384388
Nyquist_counts <- 
  data.table::fread('reference_mapping_data/MIT_Milk_Study_Raw_counts.txt') %>% 
  Matrix::as.matrix(rownames='GENE')


Nyquist_meta <- read_csv('reference_mapping_data/MIT_milk_study_metadata.csv')
# the first row is information about the type of each column, remove it here
Nyquist_meta <- Nyquist_meta[-1,] %>% column_to_rownames(var = 'NAME')







## PREVIOUS VERSION OF HUMAN GENOME BECAUSE THIS IS WHAT MIT MILK DATA USES
# IS THIS NECESSARY?

# The nyquist MIT milk data uses an old version of the human genome
# clone based assembly.
# sometimes genes are specified via a clone-based system
# AC000003 is the clone AC000003.1 is one feature on the clone AC000003.2 is another
# AC000003.2 (Clone-based (Ensembl) gene)


grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                  path="/biomart/martservice",dataset="hsapiens_gene_ensembl")

# 
# 
ATTRIBUTES <- listAttributes(grch37)
# ALL_IDs <- ATTRIBUTES %>% 
#   filter(page =='feature_page') %>% 
#   dplyr::select(name, description) %>% 
#   filter(grepl('ID', description)) %>% 
#   # filter(grepl('ensembl'))
#   pull(name)

ACCESSIONS <- 
  ATTRIBUTES %>% 
  filter(page =='feature_page') %>% 
  dplyr::select(name, description) %>%
  # filter(!(name %in% NOT_THESE)) %>%
  # filter(grepl('name|symbol', name)) %>% 
  # filter(grepl('ID', description)) %>% 
  filter(grepl('accession', name)) %>% 
  pull(name)

NAMES <- 
  ATTRIBUTES %>% 
  filter(page =='feature_page') %>% 
  dplyr::select(name, description) %>%
  # filter(!(name %in% NOT_THESE)) %>%
  # filter(grepl('name|symbol', name)) %>% 
  # filter(grepl('ID', description)) %>% 
  filter(grepl('name', name)) %>% 
  pull(name)

CLONES <- 
  ATTRIBUTES %>% 
  filter(page =='feature_page') %>% 
  dplyr::select(name, description) %>%
  # filter(!(name %in% NOT_THESE)) %>%
  # filter(grepl('name|symbol', name)) %>% 
  # filter(grepl('ID', description)) %>% 
  filter(grepl('clone_based', name)) %>% 
  pull(name)




ATTRIBUTES <- 
  c('hgnc_symbol',
    'entrezgene_accession',
    'external_gene_name',
    'rfam',
    'embl',
    'mirbase_id',
    CLONES,
    ACCESSIONS,
    NAMES
  )







nyquist_mapping_results <-
  map_IDs_to_ensembl(ID_vector = rownames(Nyquist_counts), 
                     attributes_vector = ATTRIBUTES, mart = grch37)


NUM_UNMAPPED <- nyquist_mapping_results %>% filter(is.na(mapping_type)) %>% nrow()
NUM_MAPPED <- nrow(nyquist_mapping_results)
PERCENT_MAPPED <- (NUM_MAPPED / (NUM_MAPPED + NUM_UNMAPPED)) * 100

# 
# grch37_mapping <- 
#   getBM(values=rownames(Nyquist_counts),
#         attributes= c("ensembl_gene_id",
#                       'external_gene_name',
#                       'clone_based_ensembl_gene',
#                       'clone_based_vega_gene',
#                       'description'),
#         mart= grch37) %>% as_tibble()
# 
# 
# nyquist_mapping <- tibble(nyquist_ID=rownames(Nyquist_counts), 
#                           external_gene_name=ifelse(nyquist_ID %in% grch37_mapping$external_gene_name,nyquist_ID, NA), 
#                           clone_based_ensembl_gene=ifelse(nyquist_ID %in% grch37_mapping$clone_based_ensembl_gene,nyquist_ID, NA), 
#                           clone_based_vega_gene=ifelse(nyquist_ID %in% grch37_mapping$clone_based_vega_gene,nyquist_ID, NA)) 
# 
# ensembl_external_map <- 
#   grch37_mapping %>% 
#   dplyr::select(ensembl_gene_id, external_gene_name) %>% 
#   unique() 
# 
# 
# nyquist_mapping <- nyquist_mapping %>% left_join(ensembl_external_map) 


# all but 183 gene_IDs from the nyquist study are mapped to ensembl's external_gene_name vector
# I couldn't find the remaining 183 genes in clone based mappings...
# oh well, cant be perfect... moving on.
nyquist_bovine_mapping <- 
  nyquist_mapping_results %>%
  transmute(nyquist_ID=ID, ensembl_gene_id) %>% 
  # dplyr::select(nyquist_ID, ensembl_gene_id) %>% 
  left_join(BOVINIZER) %>% 
  filter(!is.na(BOVINE_ID)) %>% 
  dplyr::select(nyquist_ID, BOVINE_ID) %>%
  left_join(our_genes)

multi_mapped <- 
  nyquist_bovine_mapping %>%
  group_by(our_id) %>%
  tally() %>%
  arrange(desc(n)) %>% 
  filter(n>1) %>%
  pull(our_id)

nyquist_bovine_mapping <- 
  nyquist_bovine_mapping %>% 
  filter(!(our_id %in% multi_mapped))
#when provided an id from the nyquiust object the nyquist swapper will return the corresponding ID 
# from our study
nyquist_swapper <- nyquist_bovine_mapping$our_id
names(nyquist_swapper) <- nyquist_bovine_mapping$nyquist_ID
# nyquist_swapper['ASPM']

# subset the nyquist counts to only contain genes that have mappings to bovine orthologs
Nyquist_counts <- Nyquist_counts[rownames(Nyquist_counts) %in% nyquist_bovine_mapping$nyquist_ID,]
# colSums(Nyquist_counts) %>% log() %>% hist()
rownames(Nyquist_counts) <- nyquist_swapper[rownames(Nyquist_counts)]
names(rownames(Nyquist_counts)) <- NULL
rownames(Nyquist_counts)



# Nyquist MIT milk data set, features subset to only bovine homologs and renamed to our bovine geneIDs
Nyquist_milk <- Seurat::CreateSeuratObject(counts=Nyquist_counts, meta.data = Nyquist_meta)
SaveH5Seurat(Nyquist_milk, 'reference_mapping_data/Nyquist_milk', overwrite = TRUE)


