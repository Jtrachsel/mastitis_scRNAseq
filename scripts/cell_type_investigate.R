library(tidyverse)
library(SeuratDisk)
library(Seurat)

library(readxl)

seu <- LoadH5Seurat('outputs/classified_clusters_30.h5seurat')


CLUSTER_MARKERS <- read_tsv('outputs/30_markers.tsv')
CLUSTER_MARKERS

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


human_homologs <- read_tsv('outputs/human_homologs.tsv')



HUMAN_MARKERS <- read_xlsx('CellMarker_DB/Cell_marker_Human.xlsx')
GENE_ID_MAPPING <- read_tsv('outputs/gene_ID_mapping.tsv') %>% mutate(gene=name)


CLUSTER_MARKERS <- CLUSTER_MARKERS %>% left_join(GENE_ID_MAPPING)










