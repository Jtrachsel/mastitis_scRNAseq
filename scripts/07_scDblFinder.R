library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

sample_dirs <- list.dirs('SoupX', recursive = F)
names(sample_dirs) <- sub('SoupX/','',sample_dirs)

SCE <- Read10X(sample_dirs)

# add 'MT-' prefix to mitochondrial genes
# O should probably find a better place to do this....
# need to be careful that genes are in the correct order here...

# ID_map %>% filter(name != ROWNAMES)%>% select(name, ROWNAMES)

# write_tsv(ID_map, 'outputs/gene_ID_mapping.tsv')

# 
# mito_names <-
#   ID_map %>%
#   filter(chromosome_name == 'MT') %>%
#   mutate(new_name=paste0('MT-', ROWNAMES))
# 
# mito_swap <- mito_names$new_name
# names(mito_swap) <- mito_names$ROWNAMES
# 
# other_names <-
#   ID_map %>%
#   filter(chromosome_name != 'MT')
# 
# other_swap <- other_names$ROWNAMES
# names(other_swap) <- other_names$ROWNAMES
# 
# all_swap <- c(mito_swap, other_swap)
# 
# NEW_ROWNAMES <- all_swap[rownames(SCE)]
# 
# names(NEW_ROWNAMES) <- NULL
# rownames(SCE) <- NEW_ROWNAMES
# 
# ### end change gene names ###
# 
# # create SingleCellExperiment object
# SCE <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=SCE))

# extract metadata from colnames
COL_DAT <-
  tibble(cell_ID=colnames(SCE),
         sample_ID=sub('([0-9]+[a-z]+)_[ATCG]+-[0-9]+','\\1',cell_ID),
         individual=sub('([0-9]+)([a-z]+)','\\1',sample_ID),
         tissue=sub('([0-9]+)([a-z]+)','\\2',sample_ID)) %>%
  DataFrame()

# assign metadata
colData(SCE) <- COL_DAT

### Run scDblFinder
set.seed(12)
results <- scDblFinder(sce = SCE, samples = 'sample_ID', clusters=TRUE)

# write results object (SingleCellExperiment object)
write_rds(results,'outputs/scDblFinder_out.rds')


