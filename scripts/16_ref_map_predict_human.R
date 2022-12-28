# Title: Reference Mapping and Prediction to human milk and whole blood data
# Date: 12-19-2022
# Author: Jayne Wiarda -- modified by Julian Trachsel
library(tidyverse)
library(Seurat)
library(SeuratDisk)
# library(readxl)
library(scales)
library(future)

# This script performs reference mapping of our cells to 2 human datasets
# Whole blood : tabula sapiens
# Breast milk : Nyquist et. al. 

# 2 strategies.

  # 1. sequential mapping
    #- Map all of our cells (both milk and blood) to the human blood data
    #- Map all of our cells (both milk and blood) to the human milk data
  # 2. Integrated reference
    #- integrate the human milk and human blood data to create an integrated reference
    #- Map all of our cells to the single integrated reference



set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=16)
} else {
  future::plan(multisession, workers=16)
}

options(future.globals.maxSize = 1000000 * 1024^2)




# Load in reference data
ref <- LoadH5Seurat("reference_mapping_data/tabula_sapiens_blood.h5seurat") # load in .h5seurat file 
DefaultAssay(ref) <- 'integrated'

# reference data is already subset to just whole blood cells and already bovinized
# meaning it's gene IDs have been converted to the same ones our data uses
ref@meta.data %>% group_by(donor) %>% tally()
# 1974 cells is the fewest cells per donor


# Load in the query data:
query <- read_rds('reference_mapping_data/reference_mapping_query.rds')


# Filter query data to include only one-to-one gene orthologs & convert to human gene symbols:

## Perform label transfer and mapping:
MappingScores <- list()
CellTypePredictions <- list()
for(i in 1:length(query)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query[[i]],
    reduction = "cca", # opted to use cca since the method is recommended for cross-species mapping 
    dims = 1:30, 
    normalization.method = "LogNormalize",
    ) 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = ref$cell_ontology_class), # predict query dataset IDs at level of reference data's cluster, lineage, and cell type classifications
                              dims = 1:30,
                              weight.reduction = "cca")
  MapScores <- MappingScore(
    anchors = anchors@anchors,
    combined.object = anchors@object.list[[1]],
    query.neighbors =  slot(object = query[[i]], name = "neighbors")[["query_ref.nn"]],
    query.weights = Tool(object = query[[i]], slot = "TransferData")$weights.matrix,
    query.embeddings = Embeddings(object = query[[i]]),
    ref.embeddings = Embeddings(object = ref),
    nn.method = "annoy",
    # n.trees = n.trees
  )
  MappingScores[[i]] <- MapScores
  CellTypePredictions[[i]] <- predictions
} 

MappingScores <- Reduce(c,MappingScores)
MappingScores <- as.data.frame(MappingScores) %>% rownames_to_column(var='CELL')
CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)
CellTypePredictions %>% group_by(predicted.id) %>% tally()
# Save the mapping & prediction results:
# MappingScores$CellBarcodes <- rownames(MappingScores)
# CellTypePredictions$CellBarcodes <- rownames(CellTypePredictions)



tabula_blood_predictions <- CellTypePredictions %>%
  rownames_to_column(var = 'CELL') %>%
  mutate(reference='tabula_blood') %>% 
  as_tibble() %>% 
  left_join(MappingScores)

tabula_blood_predictions <- 
  tabula_blood_predictions %>%
  write_tsv('outputs/tabula_blood_predictions.tsv')

# # Incorporate cell prediction & mapping scores into original Seurat object of query data:
# ep <- AddMetaData(object = ep, 
#                   metadata = c(MappingScores, CellTypePredictions))
# 
# FeaturePlot(ep, features = 'MappingScores', 
#             reduction = 'tsne', 
#             cols = c('yellow', 'red'),
#             min.cutoff = 0.6, 
#             max.cutoff = 1)
# 
# predict <- colnames(ep@meta.data %>% select(starts_with("prediction.score."))) # extract 
# names <- sub(".*prediction.score. ", "", predict)   
# FeaturePlot(ep,
#             features = c(predict),
#             reduction = 'tsne',
#             ncol = 5) & 
#   scale_color_gradientn( colours = c('yellow', 'red'),  limits = c(0, 1), oob = squish) & 
#   NoAxes() & NoLegend() 
# names

tabula_blood_filt$cell_ontology_class %>% unique()
tabula_blood_filt$free_annotation %>% unique()

Nyquist_milk$cell_type__ontology_label %>% unique()
Nyquist_milk$General_Celltype %>% unique()



#### Nyquist milk predictions



ref <- LoadH5Seurat("reference_mapping_data/nyquist_milk.h5seurat") # load in .h5seurat file 
DefaultAssay(ref) <- 'integrated'
# Filter query data to include only one-to-one gene orthologs & convert to human gene symbols:
ref$General_Celltype %>% unique()
## Perform label transfer and mapping:
MappingScores <- list()
CellTypePredictions <- list()
for(i in 1:length(query)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query[[i]],
    reduction = "cca", # opted to use cca since the method is recommended for cross-species mapping 
    dims = 1:30, 
    normalization.method = "LogNormalize",
  ) 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = ref$General_Celltype), # predict query dataset IDs at level of reference data's cluster, lineage, and cell type classifications
                              dims = 1:30,
                              weight.reduction = "cca")
  MapScores <- MappingScore(
    anchors = anchors@anchors,
    combined.object = anchors@object.list[[1]],
    query.neighbors =  slot(object = query[[i]], name = "neighbors")[["query_ref.nn"]],
    query.weights = Tool(object = query[[i]], slot = "TransferData")$weights.matrix,
    query.embeddings = Embeddings(object = query[[i]]),
    ref.embeddings = Embeddings(object = ref),
    nn.method = "annoy",
    # n.trees = n.trees
  )
  MappingScores[[i]] <- MapScores
  CellTypePredictions[[i]] <- predictions
} 



MappingScores2 <- Reduce(c,MappingScores)

MappingScores2 <-
  as.data.frame(MappingScores2) %>%
  rownames_to_column(var='CELL') %>% 
  transmute(CELL, MappingScores=MappingScores2)

CellTypePredictions2 <- do.call(rbind, CellTypePredictions)
CellTypePredictions2 <- as.data.frame(CellTypePredictions2)
CellTypePredictions2 %>% group_by(predicted.id) %>% tally()

CellTypePredictions2 %>% rownames_to_column(var='CELL')

# Save the mapping & prediction results:
# MappingScores$CellBarcodes <- rownames(MappingScores)
# CellTypePredictions$CellBarcodes <- rownames(CellTypePredictions)

# nyquist_milk_predictions %>%
#   filter(!(predicted.id %in% c('removed', 'CSN1S1 macrophages'))) %>% 
#   # mutate(predicted=) %>% 
#   ggplot(aes(x=prediction.score.max, fill=predicted.id)) + 
#   geom_histogram() + 
#   facet_wrap(~predicted.id, nrow = 2)+ scale_y_log10()

nyquist_milk_predictions <- 
  CellTypePredictions2 %>%
  rownames_to_column(var = 'CELL') %>%
  mutate(reference='nyquist_milk') %>% 
  as_tibble() %>% 
  left_join(MappingScores2)

nyquist_milk_predictions %>%
  write_tsv('outputs/nyquist_milk_predictions.tsv')



all_predictions <- 
  bind_rows(
    nyquist_milk_predictions %>%
      select(CELL,MappingScores,predicted.id,prediction.score.max, reference),
    tabula_blood_predictions %>% 
      select(CELL,MappingScores,predicted.id,prediction.score.max, reference)
    )





prediction_summary <- 
  all_predictions %>%
  group_by(CELL) %>% 
  summarise(predicted_id=predicted.id[which.max(prediction.score.max)],
            reference=reference[which.max(prediction.score.max)], 
            prediction_score=prediction.score.max[which.max(prediction.score.max)],
            score_dif=(prediction.score.max[1] - prediction.score.max[2])**2)



prediction_summary %>%
  ggplot(aes(x=prediction_score, y=score_dif, fill=reference)) + 
  geom_point(shape=21)

prediction_summary %>%
  group_by(predicted_id) %>%
  tally() %>% 
  arrange(desc(n))

new_meta <- 
  QUERY@meta.data %>%
  rownames_to_column('CELL') %>%
  select(CELL, tissue, celltype, starts_with('integrated')) %>% 
  left_join(prediction_summary)

prediction_summary %>% arrange(desc(score_dif))

new_meta %>% group_by(tissue, reference, predicted_id) %>% tally()


new_meta %>%
  group_by(celltype, predicted_id) %>% 
  tally() %>% arrange(desc(n))


new_meta %>% 
  ggplot(aes(x=celltype, y=prediction_score)) + geom_boxplot()


LOOK <- new_meta %>% group_by(tissue, reference, predicted_id) %>% tally() %>% arrange(desc(n))

LOOK %>% ggplot(aes(x=reference, y=n, fill=predicted_id)) + geom_col()


wide_predictions <- 
  all_predictions %>%
  pivot_wider(names_from = reference, 
              values_from = c(prediction.score.max,predicted.id ))


QUERY@meta.data %>%
  rownames_to_column('CELL') %>%
  left_join(wide_predictions)


#### integrated reference predictions

ref <- LoadH5Seurat("reference_mapping_data/integrated_reference.h5seurat") # load in .h5seurat file 
DefaultAssay(ref) <- 'integrated'

# 
ref@meta.data <- 
  ref@meta.data %>%
    mutate(broad_cell_type=case_when(
    grepl('macrophage', general_cell_type,ignore.case = T) ~ 'macrophage',
    grepl('neutrophil', general_cell_type,ignore.case = T) ~ 'neutrophil',
    grepl('dendritic', general_cell_type,ignore.case = T) ~ 'dendritic',
    grepl('eosinophil', general_cell_type,ignore.case = T) ~ 'eosinophils',
    grepl('fibroblasts', general_cell_type,ignore.case = T) ~ 'fibroblasts',
    grepl('erythrocyte', general_cell_type,ignore.case = T) ~ 'erythrocyte',
    grepl('b cell', general_cell_type,ignore.case = T) ~ 'b cell',
    grepl('t cell', general_cell_type,ignore.case = T) ~ 't cell',
    grepl('nk cell', general_cell_type,ignore.case = T) ~ 'nk cell',
    grepl('plasma cell', general_cell_type,ignore.case = T) ~ 'plasma cell',
    grepl('platelet', general_cell_type,ignore.case = T) ~ 'platelet',
    grepl('stem cell', general_cell_type,ignore.case = T) ~ 'stem cell',
    grepl('monocyte', general_cell_type,ignore.case = T) ~ 'monocyte',
    grepl('basophil', general_cell_type,ignore.case = T) ~ 'basophil',
    grepl('plasmablast', general_cell_type,ignore.case = T) ~ 'plasmablast',
    grepl('^LC[1-2]$', general_cell_type,ignore.case = T) ~ 'langerhans',
    TRUE ~ 'other')) %>% 
  mutate(dataset=ifelse(!is.na(manually_annotated),'tabula_blood', 'nyquist_milk'))

ref@meta.data %>% group_by(general_cell_type) %>% tally() %>% arrange((n))
ref@meta.data %>% group_by(broad_cell_type) %>% tally() %>% arrange((n))

ref@meta.data %>% group_by(dataset,broad_cell_type) %>% tally() %>% arrange(desc(n))


## Perform label transfer and mapping:
MappingScores <- list()
CellTypePredictions <- list()
for(i in 1:length(query)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query[[i]],
    reduction = "cca", # opted to use cca since the method is recommended for cross-species mapping 
    dims = 1:30, 
    normalization.method = "LogNormalize",
  ) 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = ref$general_cell_type, 
                                             broad_type= ref$broad_cell_type), # predict query dataset IDs at level of reference data's cluster, lineage, and cell type classifications
                              dims = 1:30,
                              weight.reduction = "cca")
  # predictions <- TransferData(anchorset = anchors, 
                              # refdata = list(cell_type = 'general_cell_type', 
                                             # broad_type= 'broad_cell_type'), # predict query dataset IDs at level of reference data's cluster, lineage, and cell type classifications
                              # dims = 1:30,
                              # weight.reduction = "cca")
  MapScores <- MappingScore(
    anchors = anchors@anchors,
    combined.object = anchors@object.list[[1]],
    query.neighbors =  slot(object = query[[i]], name = "neighbors")[["query_ref.nn"]],
    query.weights = Tool(object = query[[i]], slot = "TransferData")$weights.matrix,
    query.embeddings = Embeddings(object = query[[i]]),
    ref.embeddings = Embeddings(object = ref),
    nn.method = "annoy",
    # n.trees = n.trees
  )
  MappingScores[[i]] <- MapScores
  CellTypePredictions[[i]] <- predictions
} 

length(CellTypePredictions)
CellTypePredictions %>% str()

mapping_scores <- 
  MappingScores %>%
  map(.f=~as.data.frame(.x) %>%
        rownames_to_column(var='CELL') %>% 
        transmute(CELL, mapping_score=.x)) %>% 
  bind_rows()


# the predictions list is now a list of lists
# 
Integrated_raw_cell_type_predictions <- 
  CellTypePredictions %>%
  map(pluck(1)) %>%
  bind_rows() %>%
  rownames_to_column(var='CELL') %>%
  as_tibble() %>% 
  mutate(reference='integrated_raw') %>% 
  left_join(mapping_scores) %>% 
  write_tsv('outputs/integrated_predictions_raw.tsv')

Integrated_broad_cell_type_predictions <- 
  CellTypePredictions %>%
  map(pluck(2)) %>%
  bind_rows() %>%
  rownames_to_column(var='CELL') %>% 
  mutate(reference='integrated_broad') %>% 
  left_join(mapping_scores) %>% 
  write_tsv('outputs/integrated_predictions_broad.tsv')



Integrated_broad_cell_type_predictions %>% group_by(predicted.id) %>% tally()


Integrated_broad_cell_type_predictions %>%
  # filter(!(predicted.id %in% c('removed', 'CSN1S1 macrophages'))) %>% 
  # mutate(predicted=) %>% 
  ggplot(aes(x=prediction.score.max, fill=predicted.id)) + 
  geom_histogram() + 
  facet_wrap(~predicted.id)+ scale_y_log10()+
  theme(legend.position = 'none')


Integrated_broad_cell_type_predictions %>%
  # filter(!(predicted.id %in% c('removed', 'CSN1S1 macrophages'))) %>% 
  # mutate(predicted=) %>% 
  ggplot(aes(x=mapping_score, fill=predicted.id)) + 
  geom_histogram() + 
  facet_wrap(~predicted.id)+ scale_y_log10()+
  theme(legend.position = 'none')

Integrated_broad_cell_type_predictions %>% group_by(predicted.id) %>% tally() %>% arrange(n)


Integrated_broad_cell_type_predictions %>% 
  ggplot(aes(y=mapping_score, x=prediction.score.max))+
  geom_point(aes(color=predicted.id), alpha=.025) +
  facet_wrap(~predicted.id) + 
  xlim(0,1)+
  ylim(0,1) 



