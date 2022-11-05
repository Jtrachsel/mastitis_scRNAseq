library(tidyverse)

# setwd('/90daydata/fsep_004/john_lippolis_milk_blood/')

cellranger_dirs <- 
  list.dirs('cellranger_out',recursive = FALSE) 

cellranger_results <- 
  tibble(
  DIR=cellranger_dirs, 
  metrics_path=paste0(DIR,'/outs/metrics_summary.csv'),
  CSV=map(.x=metrics_path, .f=~read_csv(.x))
  ) %>% 
  unnest(CSV) 

cellranger_results <- 
  cellranger_results %>% 
  mutate(individual=sub('.*/([0-9]+)([a-z]+)','\\1',DIR), 
         tissue=sub('.*/([0-9]+)([a-z]+)','\\2',DIR))

cellranger_results %>% 
  ggplot(aes(x=individual, y=`Estimated Number of Cells`, fill=tissue))+
  geom_col(position = position_dodge()) + 
  geom_text(aes(label=`Estimated Number of Cells`), size=3, position = position_dodge(width = 1))
  # facet_wrap(~individual) 

ggsave('outputs/figures/Estimated_cells.jpeg', width=7, height=5, units='in', bg='white')


cellranger_results %>% 
  ggplot(aes(x=tissue, y=`Mean Reads per Cell`, fill=tissue))+
  geom_point(shape=21,position = position_dodge(width=1), aes(size=`Estimated Number of Cells`)) +
  ylim(0,40000)+ 
  # geom_text_repel(aes(label=`Estimated Number of Cells`), position = position_dodge(width = 1), )+
  facet_wrap(~individual)

ggsave('outputs/figures/Estimated_reads_per_cell.jpeg', width=7, height=5, units='in', bg='white')

cellranger_results %>% 
  select(individual, tissue,  everything(), -DIR, -metrics_path) %>% 
  write_tsv('outputs/first_cellranger_results.tsv')




