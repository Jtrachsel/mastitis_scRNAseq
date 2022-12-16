library(Seurat)
library(future)
library(SeuratDisk)
library(clustree)

set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

# seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20221129_JEW_IntegratedSeurat.h5seurat')
# seu <- LoadH5Seurat('outputs/20221129_JEW_IntegratedSeurat.h5seurat')
seu <- LoadH5Seurat('outputs/20221205_JEW_IntegratedSeurat.h5seurat')

seu


DefaultAssay(seu) <- 'integrated'
seu <- FindClusters(seu, 
                    resolution = c(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10), 
                    verbose = FALSE) 

# clustree(seu, # observe patterns of clustering with clustree
#          prefix = "integrated_snn_res.")

# ggsave()

DefaultAssay(seu) <- "RNA"
clustree(seu, 
         prefix = "integrated_snn_res.",
         node_colour = "CD3E",  #CD3E, CD79A, SIRPA, CSF3R
         node_colour_aggr = "mean")

# ggsave()
#DimPlot(seu, group.by = 'sample_ID')
#DimPlot(seu, group.by = 'tissue')
#FeaturePlot(seu,
#            features = c('nCount_RNA', 'nFeature_RNA'), max.cutoff = 20000)
#FeaturePlot(seu,
#            features = c('PTPRC', 
#                         'LALBA', 'CSN2', 
#                         'CD3E', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', 
#                         'CD79A', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 
#                         'BOLA-DRA', 'SIRPA', 'CD68', 'FCGR3A', 'CSF1R', 'CST3', 'NLRP3', 'CD163', 'TLR4', 
#                         'PCLAF', 'UBE2C'),
#            cols = c('tan', 'red4'),
#            ncol = 7) & NoLegend() & NoAxes()

DimPlot(seu, group.by = 'sample_ID', reduction = 'tsne', shuffle = TRUE)
# ggsave('outputs/figures/')

DimPlot(seu, group.by = 'tissue', reduction = 'tsne', shuffle = TRUE)
# ggsave()

FeaturePlot(seu,
            features = c('nCount_RNA', 'nFeature_RNA'), 
            max.cutoff = 10000, 
            reduction = 'tsne',
            cols = c('grey90', 'darkslateblue'))


ggsave('outputs/figures/005_UMIs_genes_FEATPLOT.jpeg', height=5, width = 8, units = 'in', bg='white')

# ADD  FCGR3A (aka CD16)
## jules add

##

MARKERS <- list(
  leukocyte = c('PTPRC'), 
  epithelial = c('LALBA', 'CSN2', 'LGALS3', 'SPP1', 'FABP5'), 
  T_ILC = c('CD3E', 'CD247', 'ZAP70', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A'), 
  B_ASC = c('CD79A', 'CD19', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 'XBP1'), 
  pan_myeloid = c('SIRPA', 'CD68', 'TYROBP'),
  macrophage = c('CSF1R', 'CST3', 'CD83', 'CD86', 'CD14', 'TLR4'), 
  DCs =        c('ENSBTAG00000038397', 'ENSBTAG00000009656', 'ENSBTAG00000013919', 'ENSBTAG00000037605', 'BOLA-DQA5', 'BOLA-DRA', 'BOLA-DQB','FCGR3A'),
  pDC = c('GZMB', 'IL3RA'), 
  neutrophil = c('CSF3R', 'CXCL8', 'SRGN'),
  RBCs   = c('AHSP', 'HBB', 'HBM'), 
  cell_cycle = c('PCLAF', 'UBE2C')
  )

FeaturePlot(seu,
            features = c('PTPRC', # leukocyte
                         'LALBA', 'CSN2', 'LGALS3', 'SPP1', 'FABP5', # epithelial...FABP5/SPP1 could also be unique macrophage markers???
                         'CD3E', 'CD247', 'ZAP70', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', # T/ILC
                         'CD79A', 'CD19', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 'XBP1', # B/ASC
                         'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                         'CSF1R', 'CST3', 'CD83', 'CD86', 'CD14', 'TLR4', # macrophage/monocyte
                         'ENSBTAG00000038397', 'ENSBTAG00000009656', 'ENSBTAG00000013919', 'ENSBTAG00000037605', 'BOLA-DQA5', 'BOLA-DRA', 'BOLA-DQB','FCGR3A', # DC
                         'GZMB', 'IL3RA', # pDC...also CD4+ CD14- 
                         'CSF3R', 'CXCL8', 'SRGN', #neutrophil
                         'AHSP', 'HBB', 'HBM', # RBC genes
                         'PCLAF', 'UBE2C'), # cell cycle
            cols = c('grey90', 'darkslateblue'),
            ncol = 5, 
            reduction = 'tsne') & NoLegend() & NoAxes()

ggsave('outputs/figures/006_ALL_FEATURES.jpeg', height=12, width = 7, units = 'in', bg='white')

### write out individual feature plots  

FEAT_TIB <- 
  enframe(MARKERS) %>% 
  unnest(cols = value) %>% 
  mutate(PATH=paste0('outputs/figures/FEATS_', value, '.jpeg'), 
         FEAT_PLOT=map(.x=value, .f=~FeaturePlot(seu,reduction = 'tsne', features = .x) & NoLegend() & NoAxes()), 
         GGSAVE=map2(.x=PATH, .y=FEAT_PLOT, .f=~ggsave(.x, .y, width = 3, height = 3, units = 'in', bg='white')))






names(MARKERS) != 'RBCs'
# grouped dotplot
DotPlot(seu,
        features = MARKERS[names(MARKERS) != 'RBCs'],
        cols = c('gold', 'red3'),
        group.by = 'integrated_snn_res.6') & RotatedAxis()

ggsave('outputs/figures/007_grouped_dotplot.jpeg', width = 15, height = 8, units = 'in',bg='white' )



# DotPlot(seu,
#             features = c('PTPRC', # leukocyte
#                          'LALBA', 'CSN2', 'LGALS3', 'SPP1', 'FABP5', # epithelial...FABP5/SPP1 could also be unique macrophage markers???
#                          'CD3E', 'CD247', 'ZAP70', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', # T/ILC
#                          'CD79A', 'CD19', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 'XBP1', # B/ASC
#                          'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
#                          'CSF1R', 'CST3', 'CD83', 'CD86', 'CD14', 'TLR4', # macrophage/monocyte
#                          'ENSBTAG00000038397', 'ENSBTAG00000009656', 'ENSBTAG00000013919', 'ENSBTAG00000037605', 'BOLA-DQA5', 'BOLA-DRA', 'BOLA-DQB','FCGR3A', # DC
#                          'GZMB', 'IL3RA', # pDC...also CD4+ CD14- 
#                          'CSF3R', 'CXCL8', 'SRGN', #neutrophil
#                          'AHSP', 'HBB', 'HBM', # RBC genes
#                          'PCLAF', 'UBE2C'), # cell cycle
#             cols = c('gold', 'red3'),
#         group.by = 'integrated_snn_res.6') & RotatedAxis()

DimPlot(seu, group.by = 'integrated_snn_res.6', reduction = 'tsne', label = TRUE)

ggsave('outputs/figures/008_Dimplot.jpeg', width = 8, height = 5, units='in', bg='white')

seu$celltype <- seu$integrated_snn_res.6
Idents(seu) <- seu$celltype
types <- rev(c('ep',
           'mac',
           'neu',
           'T',
           'unknown',
             'neu',
           'neu',
           'neu',
           'unknown',
             'neu',
           'T',
           'T',
           'ep',
           'mac',
           'neu',
           'neu',
           'neu',
           'B/mac mix proliferating',
           'T',
           'mac',
           'T',
           'cDC',
           'B',
           'pDC',
           'B',
           'T',
           'neu',
           'neu',
           'neu',
           'mac',
           'T',
           'T',
           'neu',
           'neu',
           'T',
           'mac',
           'neu',
           'B',
           'T',
           'B',
           'B',
           'T',
           'B',
           'neu',
           'neu',
           'B',
           'B',
           'B',
           'T',
           'neu',
           'neu',
           'T',
           'neu',
           'neu',
          'neu',
           'neu',
           'mac',
           'neu',
           'neu',
           'neu',
           'neu',
          'T',
          'T',
           'neu',
           'mac',
           'neu',
          'T',
           'neu',
           'neu',
           'B',
           'mac',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu',
           'neu')) # Rename clusters based on phenotype IDs
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$celltype <- Idents(seu)
Idents(seu) <- seu$celltype


DimPlot(seu, reduction = 'tsne', label = TRUE)

ggsave('outputs/figures/009_classified_TSNE.jpeg', width = 8, height=5, units = 'in', bg='white')

DimPlot(seu, reduction = 'umap', label = TRUE)

ggsave('outputs/figures/010_classified_UMAP.jpeg', width = 8, height=5, units = 'in', bg='white')


DimPlot(seu, reduction = 'tsne', label = TRUE, split.by = 'tissue')
ggsave('outputs/figures/011_classified_TSNE_tissue.jpeg', width = 8, height=5, units = 'in', bg='white')


Idents(seu) <- seu$celltype
tissueTotalCells <- prop.table(table(seu$tissue)) # What percent of total cells are from each tissue?
tissuePercents <- prop.table(table(Idents(seu),seu$tissue), 
                                margin = 1) # What percent of cells from each cluster belong to each region?
tissuePercents <- rbind(tissuePercents, tissueTotalCells) # add row of overall percentages to table
#rowSums(tissuePercents) # make sure all row sums are equal to 1
tissuePercents <- t(tissuePercents) # transpose the table
# barplot(tissuePercents, # create stacked bar plot
#         col = c('red', 'blue'), 
#         legend = rownames(tissuePercents),
#         ylab = "Frequency within cell type", 
#         las = 2,
#         border = NA,
#         space = 0.05,
#         legend.text = TRUE, 
#         args.legend = list(x = "topright", bty = "n"))


tissuePercents %>%
  as.data.frame() %>%
  rownames_to_column(var='tissue') %>% 
  pivot_longer(cols = -tissue, names_to = 'celltype', values_to = 'proportion') %>% 
  ggplot(aes(x=celltype, y=proportion, fill=tissue)) + geom_col(color='black') + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10, color='black')) + ggtitle('proportion cell types by tissue')

ggsave('outputs/figures/012_celltypes_tissue.jpeg', height=5, width =8, units = 'in', bg='white')

# Stacked bar of cell types within each original sample:
# Idents(seu) <- seu$sample_ID
# tissueTotalCells <- prop.table(table(seu$celltype)) # What percent of total cells are from each region?
# tissuePercents <- prop.table(table(Idents(seu),seu$celltype), 
#                                 margin = 1) # What percent of cells from each cluster belong to each region?
# tissuePercents <- rbind(tissuePercents, tissueTotalCells) # add row of overall percentages to table
# #rowSums(tissuePercents) # make sure all row sums are equal to 1
# tissuePercents <- t(tissuePercents) # transpose the table
# 


# barplot(tissuePercents, # create stacked bar plot
#         col = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta', 'black', 'grey', 'red'), 
#         legend = rownames(tissuePercents),
#         ylab = "Frequency within cell type", 
#         las = 2,
#         border = NA,
#         space = 0.05,
#         legend.text = TRUE, 
#         args.legend = list(x = "topright", bty = "n"))
# 
# 


library(tidyverse)


percent_abunds <- 
  seu@meta.data %>%
  group_by(sample_ID, tissue, individual, celltype) %>% 
  tally() %>% 
  mutate(percent_sample = (n/sum(n)) * 100)

percent_abunds %>%
  ggplot(aes(x=sample_ID, y=percent_sample, fill=celltype)) +
  geom_col(color='black')+
  facet_wrap(~ tissue, scales = 'free') +
  ggtitle('cell type percentages')


ggsave('outputs/figures/013_celltypes_sample.jpeg', height=5, width =8, units = 'in', bg='white')

seu@meta.data %>%
  filter(celltype != 'neu') %>% 
  group_by(sample_ID, tissue, individual, celltype) %>% 
  tally() %>% 
  mutate(percent_sample = (n/sum(n)) * 100) %>%
  ggplot(aes(x=sample_ID, y=percent_sample, fill=celltype)) +
  geom_col(color='black')+
  facet_wrap(~ tissue, scales = 'free') + 
  ggtitle('cell type percentages without neu')

ggsave('outputs/figures/014_celltypes_sample_NOnue.jpeg', height=5, width =8, units = 'in', bg='white')

# tissuePercents %>% as.data.frame() %>%
#   rownames_to_column(var='cell_type')
#   


Idents(seu) <- seu$integrated_snn_res.6
de61 <- FindMarkers(seu, ident.1 = 61)
de80 <- FindMarkers(seu, ident.1 = 80)














FeaturePlot(seu,
            features = c('NTRK2', 'ELF5', 'CHRDL2', 'OLAH', 'CA6', #LC2
                         'ATF3', 'MALAT1', 'JUN', 'MAFF', 'IER2', #LC1
                         'CTSL', 'CD68', 'CD163', 'LRP1', 'EMILIN2', 'LAPTM5', 'TYROBP', 'ITGB2', 'LRP1', #macrophages
                         'CSF3R', 'CXCL8', 'SRGN', 'TREM1', 'AQP9', #NEUTROPHIL
                         'BCL11A', 'DCK', 'ZBTB33', 'NOPCHAP1', #B CELL; NOPCHAP1 = C12orf45
                         'ETS1', 'SPOCK2', 'PTPRC', 'FYN', # T CELL
                         'NR4A3', 'CD86', 'GPR183', 'CD83', 'RAB31', 'CADM1',
'CLEC9A', 'IDO1', 'C1orf54', 'BATF3', 'SLAMF8', 'SNX22', 'CPNE3',
'GCSAM', 'THBD', 'CLNK', #DC
                         'HSPA6', 'PTN', 'THY1', 'SPARC', 'PTPRZ1', #FIBROBLASWT
                         'SORL1', 'CORO1A', 'SYNE1', 'CLC', 'EMR1'),
            reduction = 'tsne', ncol = 9) & NoLegend() & NoAxes()


ggsave('outputs/figures/015_celltypes_sample.jpeg', height=5, width =8, units = 'in', bg='white')


# SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20221129_JEW_IntegratedSeurat.h5seurat', overwrite = TRUE)
SaveH5Seurat(seu, 'outputs/20221129_JEW_IntegratedSeurat.h5seurat', overwrite = TRUE)


