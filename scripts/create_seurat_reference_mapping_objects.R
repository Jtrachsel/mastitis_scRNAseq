# install.packages("anndata")
library(anndata)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(biomaRt)

# function to parse the gencode.v30.annotation.gtf
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
# tabula used GENCODE GRCh38

#### for reference mapping ####


# This gtf annotation file was used by the tabula sapiens atlas
# it contains information to convert between the gene IDs used and emsembl_IDs
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz

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

# some problems converting anndata to seurat formats
# https://github.com/mojaveazure/seurat-disk/issues/109

tabula_blood_anndata <- anndata::read_h5ad('reference_mapping_data/tabula_10X_blood.h5ad')

tabula_blood_counts <- Matrix::t(tabula_blood_anndata$X)

# remove 0 count genes
tabula_blood_counts <- tabula_blood_counts[rowSums(tabula_blood_counts) != 0,]


# strips the version decimal off the ensembl gene id
ENCODE30_mapping <-
  ENCODE30_GTF %>%
  dplyr::transmute(gene_name, ensembl_gene_id=sub('([A-Z0-9]+)\\.[0-9]+','\\1',gene_id))

# join the gtf based mappings to the IDs present in the data
tabula_blood_mapping <- 
  tibble(gene_name=rownames(tabula_blood_counts)) %>% 
  left_join(ENCODE30_mapping)

# join in the data linking one to one homologs with our data
bovinized_tabula_blood <- 
  tabula_blood_mapping %>% 
  left_join(BOVINIZER) %>%
  filter(!is.na(BOVINE_ID)) %>% 
  left_join(our_genes)


# make a little "dictionary" that converts 'tabula sapiens IDs' to 'our IDs'
tabula_swapper <- bovinized_tabula_blood$our_id
names(tabula_swapper) <- bovinized_tabula_blood$gene_name

# subset the data to only genes with bovine homologs
tabula_blood_counts <- tabula_blood_counts[rownames(tabula_blood_counts) %in% bovinized_tabula_blood$gene_name,]


# swap the gene names to the cattle gene names
rownames(tabula_blood_counts) <- tabula_swapper[rownames(tabula_blood_counts)]

# remove the names from the rownames
names(rownames(tabula_blood_counts)) <- NULL
rownames(tabula_blood_counts)

# make a seurat object from the renamed counts and the metadata from the anndata object
tabula_blood <- CreateSeuratObject(counts = tabula_blood_counts, meta.data = tabula_blood_anndata$obs)

# double check that we dont have any genes with all zeros
non_zero_genes <- rowSums(tabula_blood) > 0
tabula_blood_filt <- tabula_blood[non_zero_genes,]

# write out object
SaveH5Seurat(tabula_blood_filt, 'reference_mapping_data/tabula_sapiens_blood', overwrite=TRUE)



###

### MIT human breast milk study 
# nyquist used GRCh37
# past version of human genome

###

Nyquist_counts <- 
  data.table::fread('reference_mapping_data/MIT_Milk_Study_Raw_counts.txt') %>% 
  Matrix::as.matrix(rownames='GENE')


Nyquist_meta <- read_csv('reference_mapping_data/MIT_milk_study_metadata.csv')
# the first row is information about the type of each column, remove it here
Nyquist_meta <- Nyquist_meta[-1,] %>% column_to_rownames(var = 'NAME')




## PREVIOUS VERSION OF HUMAN GENOME BECAUSE THIS IS WHAT MIT MILK DATA USES

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


ACCESSIONS <- 
  ATTRIBUTES %>% 
  filter(page =='feature_page') %>% 
  dplyr::select(name, description) %>%
  filter(grepl('accession', name)) %>% 
  pull(name)

NAMES <- 
  ATTRIBUTES %>% 
  filter(page =='feature_page') %>% 
  dplyr::select(name, description) %>%
  filter(grepl('name', name)) %>% 
  pull(name)

CLONES <- 
  ATTRIBUTES %>% 
  filter(page =='feature_page') %>% 
  dplyr::select(name, description) %>%
  filter(grepl('clone_based', name)) %>% 
  pull(name)



# list of attributes to try checking these IDS against
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






# run the mapping function
nyquist_mapping_results <-
  map_IDs_to_ensembl(ID_vector = rownames(Nyquist_counts), 
                     attributes_vector = ATTRIBUTES, mart = grch37)

# calculate % IDs mapped
NUM_UNMAPPED <- nyquist_mapping_results %>% filter(is.na(mapping_type)) %>% nrow()
NUM_MAPPED <- nrow(nyquist_mapping_results)
PERCENT_MAPPED <- (NUM_MAPPED / (NUM_MAPPED + NUM_UNMAPPED)) * 100

# 99%+

nyquist_bovine_mapping <- 
  nyquist_mapping_results %>%
  transmute(nyquist_ID=ID, ensembl_gene_id) %>% 
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


