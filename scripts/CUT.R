
### END PASTE
### SPLIT INVESTIGATION FROM THE SCTRANSFORM VIGNET
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
# InstallData("ifnb")
# load dataset
# ifnb <- LoadData("ifnb")

# ifnb@meta.data$orig.ident %>% table()


# split the dataset into a list of two seurat objects (stim and CTRL)
# ifnb.list <- SplitObject(ifnb, split.by = "stim")

# ctrl <- ifnb.list[["CTRL"]]
# stim <- ifnb.list[["STIM"]]




### traditional clustering route

seu_filt <- NormalizeData(seu_filt, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)

# ID highly variable features
seu_filt <- FindVariableFeatures(seu_filt, selection.method = "vst", nfeatures = 7500)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu_filt), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu_filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# scale data
all.genes <- rownames(seu_filt)
seu_filt <- ScaleData(seu_filt, features = all.genes)


seu_filt <- RunPCA(seu_filt, features = VariableFeatures(object = seu_filt))

# lots of ribosomal proteins for PC1....
VizDimLoadings(seu_filt, dims = 1:2, reduction = "pca")


DimPlot(seu_filt, reduction = "pca", group.by = 'individual', split.by = 'tissue')

# seu_filt@reductions$pca
DimHeatmap(seu_filt, dims = 1:6, cells = 500, balanced = TRUE)


seu_filt <- JackStraw(seu_filt, dims=50,num.replicate = 100)
seu_filt <- ScoreJackStraw(seu_filt, dims = 1:50)


JackStrawPlot(seu_filt, dims = 20:50)
ElbowPlot(seu_filt)

###

seu_filt <- FindNeighbors(seu_filt, dims = 1:50)



seu_filt <- FindClusters(seu_filt, resolution = 0.95)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.95=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.95



seu_filt <- FindClusters(seu_filt, resolution = 0.85)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.85=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.85





seu_filt <- FindClusters(seu_filt, resolution = 0.75)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.75=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.75



seu_filt <- FindClusters(seu_filt, resolution = 0.5)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.5=seurat_clusters)
seu_filt@meta.data$seurat_clusters_0.5

seu_filt <- FindClusters(seu_filt, resolution = 0.25)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.25=seurat_clusters)

seu_filt <- FindClusters(seu_filt, resolution = 0.125)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.125=seurat_clusters)

seu_filt <- FindClusters(seu_filt, resolution = 0.0625)
seu_filt@meta.data <- 
  seu_filt@meta.data %>%
  mutate(seurat_clusters_0.0625=seurat_clusters)

seu_filt <- RunUMAP(seu_filt, dims = 1:50)
DimPlot(seu_filt, reduction = "umap", group.by = 'seurat_clusters_0.95', label = T)
DimPlot(seu_filt, reduction = "umap", group.by = 'seurat_clusters_0.95', label = T, split.by = 'tissue')
#####
# pretty much all clusters are exclusive to either milk or blood
seu_filt@meta.data %>% 
  group_by(seurat_clusters_0.95, tissue) %>% 
  tally() %>% 
  ggplot(aes(y=n, x=tissue, fill=tissue)) +
  geom_col(position = position_dodge(preserve = 'single'), color='black') + 
  facet_wrap(~seurat_clusters_0.95, nrow=1) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
