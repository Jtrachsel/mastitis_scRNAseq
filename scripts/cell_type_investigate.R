library(tidyverse)
library(SeuratDisk)
library(Seurat)

library(readxl)

seu <- LoadH5Seurat('outputs/classified_clusters_30.h5seurat')


CLUSTER_MARKERS <- 
  read_tsv('outputs/30_markers.tsv')%>% 
  filter(p_val_adj < 0.05)

# cluster 0 has IL1RN as top 
system('mkdir CellMarker_DB')
download.file('http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Mouse.xlsx',
              destfile = 'CellMarker_DB/Cell_marker_Mouse.xlsx')


download.file('http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx',
              destfile = 'CellMarker_DB/Cell_marker_Human.xlsx')



download.file('http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Seq.xlsx',
              destfile = 'CellMarker_DB/Cell_marker_SINGLE_CELL.xlsx')


# filter to 1:1 homologues
# use cluster markers and map geneIDs to cell marker DB
# get consensus cluster types


human_homologs <- read_tsv('outputs/human_homologs.tsv') %>%
  filter(hsapiens_homolog_orthology_type == 'ortholog_one2one')



HUMAN_MARKERS <- read_xlsx('CellMarker_DB/Cell_marker_Human.xlsx')



MOUSE_MARKERS <- read_xlsx('CellMarker_DB/Cell_marker_Mouse.xlsx')



GENE_ID_MAPPING <- read_tsv('outputs/gene_ID_mapping.tsv') %>% mutate(gene=name)


HUMAN_CLUSTER_MARKERS <- 
  CLUSTER_MARKERS %>% 
  left_join(GENE_ID_MAPPING) %>%
  left_join(human_homologs) %>% 
  filter(!is.na(hsapiens_homolog_associated_gene_name))

FeaturePlot(seu, features = 'LCP1', split.by = 'tissue')



# DMXL2 drives epithelial to mesenchymal transition 
# mesenchymal stem cell surface markers CD44, CD29, SCA-1 and negative for CD33, CD34, CD45, CD73



# https://www.science.org/doi/10.1126/sciadv.abm6865
# Profiling of mature-stage human breast milk cells identifies six unique lactocyte subpopulations

# RESULTS
# Epithelial lactocytes are the most abundant cell population in breast milk

# Hassiotou et al. (7) reported that up to 94% of the total breast milk cell population consists of leukocytes in the case of infection,



# KRT18, the gene associated with CK18 and lactocytes