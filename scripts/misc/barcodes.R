# investigate barcodes


library(tidyverse)


BAR <- read_lines(gzfile('cellranger_out/1312blood-1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'))
tibble(SAMPLE='1312blood-1', 
       BARCODES=BAR)


collect_barcodes <- 
  function(directory){
  filename=paste0(directory, '/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
  tibble(SAMPLE=sub('cellranger_out/(.*)/outs/filtered.*','\\1',filename), 
         BARCODES=read_lines(gzfile(filename)))
  }


collect_barcodes('cellranger_out/1312blood-1')

barcode_data <- 
  tibble(DIRS=list.dirs('cellranger_out', recursive = F),
       BARCODES=map(.x=DIRS, .f=~collect_barcodes(.x))) %>% 
  select(BARCODES) %>% 
  unnest(BARCODES) 

barcode_counts <- barcode_data %>% 
  group_by(BARCODES) %>% tally() 


four_times <- barcode_counts %>% filter(n==4) %>% pull(BARCODES)

barcode_data %>% filter(BARCODES %in% four_times) %>% group_by(SAMPLE) %>% tally()



