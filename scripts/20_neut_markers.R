library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(glue)
library(future)


if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

options(future.globals.maxSize = 64000 * 1024^2)



granulocytes <- SeuratDisk::LoadH5Seurat('outputs/granulocytes.h5seurat')

granulocytes@assays

Seurat::DefaultAssay(granulocytes)

# DimPlot(granulocytes, split.by = 'tissue')

Idents(granulocytes) <- granulocytes$phyloorder

# set up all comparisonss for DE between clusters
comparisons <- 
  combn(as.character(unique(granulocytes$phyloorder)),2) %>% 
  as.data.frame() %>% t() %>% as_tibble(.name_repair = 'unique') %>% 
  transmute(ID1=...1, ID2=...2)

# function to apply that will run the DE and collect the results in a nice format
get_DE <- function(seurat_obj, ID1, ID2){
  res <- 
    FindMarkers(ident.1 = ID1, ident.2 = ID2,object = seurat_obj ) %>% 
    rownames_to_column(var='gene') %>% 
    mutate(ID1=ID1, ID2=ID2, 
           enriched_in=ifelse(avg_log2FC > 0, ID1, ID2))
  return(res)  
}

# safe version of this function that will return NULL when there's a failure
# not all the clusters exist in both tissues, this causes errors
safe_DE <- safely(get_DE)

# apply the function to all the pairwise comparisons
# blood and milk combined
all_pw_comps <- 
  comparisons %>% 
  mutate(DE=map2(.x=ID1, .y=ID2, .f=~get_DE(granulocytes, .x, .y))) #%>%
  # mutate(num_diff=map_int(.x=DE, .f=~nrow(.x))) %>% 
  # arrange(num_diff)

all_pw_comps %>% write_rds('all_pw_comps.rds')

all_pw_comps %>% 
  select(DE) %>% 
  unnest(DE) %>% 
  write_tsv('all_pw_DE.tsv')

# milk only
milk_gran <- subset(x = granulocytes, subset = tissue == "milk")

all_pw_comps_milk <- 
  comparisons %>% 
  mutate(DE=map2(.x=ID1, .y=ID2, .f=~safe_DE(milk_gran, .x, .y))) #%>%
  # mutate(num_diff=map_int(.x=DE, .f=~nrow(.x))) %>% 
  # arrange(num_diff)

all_pw_comps_milk %>% write_rds('all_pw_comps_milk.rds')

all_pw_comps_milk %>% 
  select(DE) %>% 
  unnest(DE) %>% 
  write_tsv('all_pw_milk_DE.tsv')


# blood only

blood_gran <- subset(x = granulocytes, subset = tissue == "blood")

all_pw_comps_blood <- 
  comparisons %>% 
  mutate(DE=map2(.x=ID1, .y=ID2, .f=~safe_DE(blood_gran, .x, .y))) #%>%
  # mutate(num_diff=map_int(.x=DE, .f=~nrow(.x))) %>% 
  # arrange(num_diff)

all_pw_comps_blood %>% write_rds('all_pw_comps_blood.rds')

all_pw_comps_blood %>% 
  select(DE) %>% 
  unnest(DE) %>% 
  write_tsv('all_pw_blood_DE.tsv')


#####

all_gran_markers <- FindAllMarkers(granulocytes,)

all_gran_markers %>% write_rds('all_gran_markers.rds')
all_gran_markers %>% write_tsv('all_gran_markers.tsv')


all_gran_markers_milk <- FindAllMarkers(milk_gran)
all_gran_markers_milk %>% write_rds('all_gran_markers_milk.rds')
all_gran_markers_milk %>% write_tsv('all_gran_markers_milk.tsv')

all_gran_markers_blood <- FindAllMarkers(blood_gran)
all_gran_markers_blood %>% write_rds('all_gran_markers_blood.rds')
all_gran_markers_blood %>% write_tsv('all_gran_markers_blood.tsv')


######## gene plots

genes_of_interest_plot <- 
  function(genes_of_interest, seurat_obj){
    
    split_plot <- DotPlot(seurat_obj, features = genes_of_interest, split.by = 'tissue')
    split_dat <- split_plot$data
    
    
    tst <- 
      split_dat %>% 
      filter(!features.plot %in% c('SELL', 'MPO')) %>% 
      mutate(tissue = sub('c[0-9]+_([a-z]+)','\\1',id),
             id2 = sub('(c[0-9]+)_([a-z]+)','\\1',id),
             avg_exp_scale=scale(avg.exp), 
             avg_exp_scale_lim= case_when(
               avg_exp_scale > 2.5 ~ 2.5, 
               avg_exp_scale < -2.5 ~ -2.5, 
               TRUE ~ avg_exp_scale
             )) %>% 
      mutate(id3=factor(id2, levels = levels(granulocytes@meta.data$phyloorder)))
    # 
    
    tst %>% 
      ggplot(aes(x=tissue, y=id3, fill=avg_exp_scale_lim, size=pct.exp)) + 
      geom_point(shape=21, color='grey')+
      cowplot::theme_half_open() +
      scale_fill_gradient(low = 'white', high = 'red') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab('') +
      ylab('')+
      facet_wrap(~features.plot, nrow = 1) +
      # ggtitle('priority 1 genes') +
      labs(fill='avg expression', 
           size='percent expression')
    
  }

genes_of_interest <- read_tsv('outputs/ellie_genes.tsv')

priority_1 <- genes_of_interest %>% filter(Priority ==1) %>% pull(gene_name)
priority_2 <- genes_of_interest %>% filter(Priority ==2) %>% pull(gene_name)
priority_3 <- genes_of_interest %>% filter(Priority == 3) %>% filter(!is.na(gene_name))

priority_3_plots <- 
  priority_3 %>% 
  mutate(GROUP=ceiling(row_number()/8)) %>% 
  group_by(GROUP) %>% 
  nest() %>% 
  mutate(genes_of_interest=map(.x=data, .f=~.x %>% pull(gene_name))) %>% 
  select(-data) %>% 
  mutate(plots=map(.x=genes_of_interest, .f=~genes_of_interest_plot(.x, granulocytes)))



P1_plot <- genes_of_interest_plot(genes_of_interest = priority_1, seurat_obj = granulocytes)
P2_plot <- genes_of_interest_plot(genes_of_interest = priority_2, seurat_obj = granulocytes)

P1_plot %>% ggsave(plot = .,'dotplot_priority_1_genes.jpg',device = 'jpeg', width = 9, height = 6, units = 'in')
P2_plot %>% ggsave(plot = .,'dotplot_priority_2_genes.jpg',device = 'jpeg', width = 9, height = 6, units = 'in')

priority_3_plots %>% 
  mutate(plt_save=map2(.x=plots, .y=GROUP, 
                       .f=~.x %>% 
                         ggsave(plot = ., 
                                filename = glue('dotplot_priority_3_genes_{.y}.jpg'),
                                device = 'jpeg', 
                                width = 9,
                                height = 6,
                                units = 'in')))
# #
# # granulocytes@assays$RNA@counts
# # granulocytes@assays$RNA@data
# # look <- granulocytes@commands
# 
# # look$NormalizeData.RNA
# #############
# blood <- subset(granulocytes, subset = tissue == 'blood')
# 
# milk <- subset(granulocytes, subset = tissue == 'milk')
# 
# blood_plot <- DotPlot(blood, features = priority_1)
# milk_plot <- DotPlot(milk, features = priority_1)
# 
# milk_dat <- milk_plot$data %>% mutate(tissue = 'milk')
# blood_dat <- blood_plot$data %>% mutate(tissue ='blood')
# 
# # cluster 60 not present in blood
# c60 <- milk_dat[!(milk_dat$id %in% blood_dat$id),]
# c60$avg.exp <- NA
# c60$pct.exp <- NA
# c60$avg.exp.scaled <- NA
# c60$tissue <- 'blood'
# 
# blood_dat$avg.exp.scaled %>% hist()
# split_dat$avg.exp.scaled %>% hist()
# split_dat$avg.exp %>% scale()
# 
# blood_dat <- blood_dat %>% bind_rows(c60)
# blood_dat[!(blood_dat$id %in% milk_dat$id),]
# library(cowplot)
# 
# # blood_dat <- blood_dat %>% mutate(id2=as.character(id))
# # milk_dat <- milk_dat %>% mutate(id2=as.character(id))
# 
# 
# blood_dat <- blood_dat %>% mutate(id2=factor(id, levels = levels(milk_dat$id)))
# milk_dat <- milk_dat %>% mutate(id2=(id))
# 
# 
# p_blood <- 
#   blood_dat %>%
#   ggplot(aes(x=features.plot, y=id2, fill=avg.exp.scaled, size=pct.exp)) + 
#   geom_point(shape=21, color='white')+
#   cowplot::theme_half_open() +
#   scale_fill_gradient(low = 'white', high = 'red') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   xlab('') +
#   ggtitle('Blood')
# 
# p_milk <- 
#   milk_dat %>%
#   ggplot(aes(x=features.plot, y=id2, fill=avg.exp.scaled, size=pct.exp)) + 
#   geom_point(shape=21, color='white')+
#   cowplot::theme_half_open() +
#   scale_fill_gradient(low = 'white', high = 'red')+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab('')+
#   ggtitle('Milk')
# 
# library(patchwork)
# p_blood + p_milk + plot_layout(guides = "collect")
# 
# 
# 
# 
# 
# 
# blood_dat %>% 
#   bind_rows(milk_dat) %>% 
#   mutate(avg_exp_scale=scale(avg.exp), 
#          avg_exp_scale_lim= case_when(
#            avg_exp_scale > 2.5 ~ 2.5, 
#            avg_exp_scale < -2.5 ~ -2.5, 
#            TRUE ~ avg_exp_scale
#          )) %>% 
#   ggplot(aes(x=features.plot, y=id2, fill=avg_exp_scale_lim, size=pct.exp)) + 
#   geom_point(shape=21, color='white')+
#   cowplot::theme_half_open() +
#   scale_fill_gradient(low = 'white', high = 'red') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   xlab('') +
#   facet_wrap(~tissue)
# # 
# #   ggplot(aes(x=features.plot, y=id2, fill=avg.exp.scaled, size=pct.exp)) + 
# #   geom_point(shape=21, color='white')+
# #   theme_half_open() +
# #   scale_fill_gradient(low = 'white', high = 'red') +
# #   facet_wrap(~tissue)+
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
# #   xlab('')+
# #   ggtitle('Priority 1 genes')
# 
# 
# Seurat
