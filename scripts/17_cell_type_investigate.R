library(Seurat)
library(SeuratDisk)
library(tidyverse)


seu <- LoadH5Seurat('outputs/Integrated_classified.h5seurat')


integrated_predictions_broad <- read_tsv('outputs/integrated_predictions_broad.tsv')
integrated_predictions_raw <- read_tsv('outputs/integrated_predictions_raw.tsv')

milk_predictions <- read_tsv('outputs/nyquist_milk_predictions.tsv')
blood_predictions <- read_tsv('outputs/tabula_blood_predictions.tsv')

integrated_predictions_broad %>% colnames()
integrated_predictions_summary <- 
  integrated_predictions_broad %>% 
  transmute(CELL, broad_type=predicted.id,broad_prediction_score=prediction.score.max, integrated_mapping_score=mapping_score) %>% 
  left_join(
    integrated_predictions_raw %>% 
      transmute(CELL, raw_type=predicted.id, raw_prediction_score=prediction.score.max)
  )


rownames(seu@meta.data)
seu@meta.data$CELL <- rownames(seu@meta.data)
new_metadata <- 
  seu@meta.data %>%
  left_join(integrated_predictions_summary) %>% column_to_rownames(var='CELL')


new_metadata %>% group_by(broad_type) %>% tally() %>% arrange((n))
seu@meta.data <- new_metadata

seu_filt <- seu[,!seu@meta.data$broad_type %in% c('fibroblast', 'other')]


DimPlot(seu_filt, group.by = 'raw_type', split.by = 'tissue')


# things are a little screwy...
LOOK <- seu@meta.data %>% 
  group_by(integrated_snn_res.0.5,tissue, broad_type) %>%
  tally() %>%
  arrange(desc(n)) %>%
  slice_head(n = 2)



# Probably makes the most sense to take the blood annotations for the blood samples
# THen the milk annotations for the milk samples?



