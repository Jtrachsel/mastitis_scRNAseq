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
  filter(!is.na(BOVINE_ID)) %>%  # removes genes from the mapping file without a 1:1 bovine homolog
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
# SaveH5Seurat(tabula_blood_filt, 'reference_mapping_data/tabula_sapiens_blood', overwrite=TRUE)

# ref <- LoadH5Seurat("reference_mapping_data/tabula_sapiens_blood.h5seurat") # load in .h5seurat file 

# reference data is already subset to just whole blood cells and already bovinized
# meaning it's gene IDs have been converted to the same ones our data uses
tabula_blood_filt@meta.data %>% group_by(donor) %>% tally()
# 1974 cells is the fewest cells per donor




# Normalize & integrate data:
ref.list <- SplitObject(tabula_blood_filt, split.by = "donor") # split into the original samples that were processed for scRNA-seq
for (i in 1:length(ref.list)) { # for each sample individually, let's normalize the data and find the 2000 most highly variable features
  ref.list[[i]] <- NormalizeData(ref.list[[i]], 
                                 verbose = TRUE, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000, 
                                 assay = "RNA")
  ref.list[[i]] <- FindVariableFeatures(ref.list[[i]], 
                                        selection.method = "vst", 
                                        nfeatures = 2500, 
                                        verbose = TRUE)
}
ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                      dims = 1:40) # find integration anchors between samples based on variable features for each sample
tabula_blood_filt <- IntegrateData(anchorset = ref.anchors, 
                     dims = 1:40) # integrate the data together based on integration anchors found with default parameters

tabula_blood_filt <- ScaleData(tabula_blood_filt, 
                 verbose = TRUE, 
                 assay = 'integrated') # scale the genes in the integrated assay

# Calculate principle components:
tabula_blood_filt <- RunPCA(tabula_blood_filt, # calculate first 100 PCs
              npcs = 100, 
              verbose = TRUE)
ElbowPlot(tabula_blood_filt,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering
# pct <- ref[["pca"]]@stdev / sum(ref[["pca"]]@stdev) * 100 # find standard deviation for each PC
# cumu <- cumsum(pct) # find cumulative percentages for PCs
# co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
# pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
#plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
#                      cumu = cumu, 
#                      rank = 1:length(pct))
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
# PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
#length(PCdims) # how many significant PCs are there?
PCdims <- 1:30 # JT

# Run UMAP dimensionality reduction:
tabula_blood_filt <- RunUMAP(tabula_blood_filt, dims = PCdims, reduction = "pca") # create UMAP
DimPlot(tabula_blood_filt, 
        group.by = "donor") # plot by sample ID




# Save reference as an .h5seurat object:
SaveH5Seurat(tabula_blood_filt, 'reference_mapping_data/tabula_sapiens_blood', overwrite=TRUE)
###



##### MILK DATA #####


### MIT human breast milk study 
# nyquist used GRCh37
# past version of human genome
rm(ref)
rm(ref.anchors)
rm(ref.list)



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


# SaveH5Seurat(Nyquist_milk, 'reference_mapping_data/Nyquist_milk', overwrite = TRUE)


# ref <- LoadH5Seurat("reference_mapping_data/tabula_sapiens_blood.h5seurat") # load in .h5seurat file 

# reference data is already subset to just whole blood cells and already bovinized
# meaning it's gene IDs have been converted to the same ones our data uses
Nyquist_milk@meta.data %>% group_by(biosample_id) %>% tally() %>% arrange(desc(n))
# 172 cells is the fewest cells per donor
# Nyquist_milk@meta.data$biosample_id %>% unique()
Nyquist_milk@meta.data %>% colnames()
Nyquist_milk@meta.data %>% group_by(biosample_id)
Nyquist_milk@meta.data %>% group_by(donor_id) %>% tally()


Nyquist_milk@meta.data <- 
  Nyquist_milk@meta.data %>%
  mutate(DONOR_ID=ifelse(is.na(donor_id), 'BM00', donor_id))
# Normalize & integrate data:
ref.list <- SplitObject(Nyquist_milk, split.by = "DONOR_ID") # split into the original samples that were processed for scRNA-seq
for (i in 1:length(ref.list)) { # for each sample individually, let's normalize the data and find the 2000 most highly variable features
  ref.list[[i]] <- NormalizeData(ref.list[[i]], 
                                 verbose = TRUE, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000, 
                                 assay = "RNA")
  ref.list[[i]] <- FindVariableFeatures(ref.list[[i]], 
                                        selection.method = "vst", 
                                        nfeatures = 2000, 
                                        verbose = TRUE)
}
ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                      dims = 1:30) # find integration anchors between samples based on variable features for each sample
Nyquist_milk <- IntegrateData(anchorset = ref.anchors, 
                     dims = 1:30) # integrate the data together based on integration anchors found with default parameters

Nyquist_milk <- ScaleData(Nyquist_milk, 
                 verbose = TRUE, 
                 assay = 'integrated') # scale the genes in the integrated assay

# Calculate principle components:
Nyquist_milk <- RunPCA(Nyquist_milk, # calculate first 100 PCs
              npcs = 100, 
              verbose = TRUE)
ElbowPlot(Nyquist_milk,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering
# pct <- ref[["pca"]]@stdev / sum(ref[["pca"]]@stdev) * 100 # find standard deviation for each PC
# cumu <- cumsum(pct) # find cumulative percentages for PCs
# co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
# pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
#plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
#                      cumu = cumu, 
#                      rank = 1:length(pct))
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
# PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
#length(PCdims) # how many significant PCs are there?
PCdims <- 1:30 # JT

# Run UMAP dimensionality reduction:
Nyquist_milk <- RunUMAP(Nyquist_milk, dims = PCdims, reduction = "pca") # create UMAP
DimPlot(Nyquist_milk, group.by = 'DONOR_ID', shuffle = T)

Nyquist_milk


# Save reference as an .h5seurat object:
SaveH5Seurat(Nyquist_milk,filename = 'reference_mapping_data/nyquist_milk', overwrite=TRUE)


# create query objects from our data

QUERY <- LoadH5Seurat('outputs/Integrated_classified.h5seurat')


DefaultAssay(QUERY) <- 'RNA'
# remove unneeded assays to reduce the size of the object
QUERY[['SCT']] <- NULL
QUERY[['integrated']] <- NULL


# Filter query data to include only one-to-one gene orthologs 

bovine_human_homologs <- BOVINIZER %>% left_join(our_genes)

# this is our data subset to only genes that are one to one homologs with human genes
QUERY <- QUERY[rownames(QUERY) %in% bovine_human_homologs$our_id,]
QUERY$sample_ID
query <- SplitObject(QUERY, split.by = "sample_ID") # split into the original samples that were processed for scRNA-seq
for (i in 1:length(query)) { # for each sample individually, let's normalize the data and find the 2000 most highly variable features, scale the data, and find top 50 PCs
  query[[i]] <- NormalizeData(query[[i]], 
                              verbose = TRUE, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000, 
                              assay = "RNA")
  query[[i]] <- FindVariableFeatures(query[[i]], 
                                     selection.method = "vst", 
                                     nfeatures = 2000, 
                                     verbose = TRUE)
  query[[i]] <- ScaleData(query[[i]]) # scale the data
  query[[i]] <- RunPCA(query[[i]], # calculate 50 PCs
                       npcs = 50, 
                       verbose = TRUE)
  DefaultAssay(query[[i]]) <- 'RNA'
}


write_rds(query, 'reference_mapping_data/reference_mapping_query.rds')



# build integrated reference
# make same column name for DONOR_ID


tabula_blood_filt <- LoadH5Seurat('reference_mapping_data/tabula_sapiens_blood.h5seurat')
Nyquist_milk <- LoadH5Seurat('reference_mapping_data/nyquist_milk.h5seurat')

tabula_blood_filt$free_annotation %>% unique() # this one
Nyquist_milk$General_Celltype %>% unique() # this one


tabula_blood_filt@meta.data$DONOR_ID <- tabula_blood_filt@meta.data$donor
Nyquist_milk@meta.data$DONOR_ID

###
reference_combined <- merge(Nyquist_milk, y = tabula_blood_filt, add.cell.ids = c("MILK", "BLOOD"), project = "integrated_reference")
DefaultAssay(reference_combined) <- 'RNA'

reference_combined[['integrated']] <- NULL


### INTEGRATE BLOOD AND MILK REFERENCES INTO ONE REFERENCE

ref.list <- SplitObject(reference_combined, split.by = "DONOR_ID") # split into the original samples that were processed for scRNA-seq
for (i in 1:length(ref.list)) { # for each sample individually, let's normalize the data and find the 2000 most highly variable features
  ref.list[[i]] <- NormalizeData(ref.list[[i]], 
                                 verbose = TRUE, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000, 
                                 assay = "RNA")
  ref.list[[i]] <- FindVariableFeatures(ref.list[[i]], 
                                        selection.method = "vst", 
                                        nfeatures = 4000, 
                                        verbose = TRUE)
}
ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                      dims = 1:40) # find integration anchors between samples based on variable features for each sample
reference_combined <- IntegrateData(anchorset = ref.anchors, 
                                   dims = 1:40) # integrate the data together based on integration anchors found with default parameters

reference_combined <- ScaleData(reference_combined, 
                               verbose = TRUE, 
                               assay = 'integrated') # scale the genes in the integrated assay

# Calculate principle components:
reference_combined <- RunPCA(reference_combined, # calculate first 100 PCs
                            npcs = 100, 
                            verbose = TRUE)
ElbowPlot(reference_combined,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering
PCdims <- 1:30 # JT

# Run UMAP dimensionality reduction:
reference_combined <- RunUMAP(reference_combined, dims = PCdims, reduction = "pca") # create UMAP
DimPlot(reference_combined, 
        group.by = "donor") # plot by sample ID


reference_combined$General_Celltype %>% unique()
reference_combined$cell_type__ontology_label %>% unique()

reference_combined@meta.data <- 
  reference_combined@meta.data %>%
  mutate(general_cell_type=ifelse(!is.na(General_Celltype),General_Celltype, free_annotation))
reference_combined@meta.data %>% pull(general_cell_type) %>% unique()


SaveH5Seurat(reference_combined, filename = 'reference_mapping_data/integrated_reference')

