library(SoupX)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(patchwork)
# need devel version for cellranger 7.0+
# remotes::install_github("constantAmateur/SoupX",ref='devel')
usethis::use_directory('SoupX')


run_soupx <- function(cellranger_directory){
  sc = load10X(cellranger_directory)
  SAMPLE_ID <- sub('.*/([0-9]+[a-z]+)/outs/','\\1',cellranger_directory)
  png(filename = paste0('SoupX/',SAMPLE_ID,'rho_.png'), width = 7, height = 5, units = 'in', res = 72)
  sc = autoEstCont(sc)
  dev.off()
  out = adjustCounts(sc, roundToInt=TRUE) # as suggested in https://github.com/satijalab/sctransform/issues/139
  
  
  THEME <- theme(legend.text = element_text(size=5), 
                 legend.title = element_text(size=6),
                 panel.grid = element_blank(), 
                 axis.text = element_blank(), 
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line = element_blank(),
                 panel.border = element_rect(color='black', fill = NA))

  B1 <- plotMarkerMap(sc, "CD3E") + ggtitle('CD3E') + THEME
  # B2 <- plotMarkerMap(sc, "ENSSSCG00000038719")
  B3 <- plotMarkerMap(sc, "CD79B")+ ggtitle("CD79B") + THEME
  B4 <- plotMarkerMap(sc, "TYROBP")+ ggtitle("TYROBP") + THEME
  B5 <- plotMarkerMap(sc, "CLEC12A")+ ggtitle('CLEC12A') + THEME
  B6 <- plotMarkerMap(sc, "BOLA-DRA")+ ggtitle('BOLA-DRA') + THEME
  B7 <- plotMarkerMap(sc, "ITGAL")+ ggtitle('ITGAL') + THEME
  combined <- (B1 +B3 + B4 ) / (B5+B6+B7)
  combined + plot_layout(guides = "collect")
  # 
  ggsave(paste0('./SoupX/',SAMPLE_ID, '_pre_soupX.jpeg'), width = 7, height = 5, units = 'in', bg = 'white')
  
  C1 <- plotChangeMap(sc,out, "CD3E") + ggtitle('CD3E') +THEME
  C3 <- plotChangeMap(sc,out, "CD79B")+ ggtitle('CD79B')+THEME
  C4 <- plotChangeMap(sc,out, "TYROBP")+ ggtitle('TYROBP')+THEME
  C5 <- plotChangeMap(sc,out, "CLEC12A")+ ggtitle('CLEC12A')+THEME
  C6 <- plotChangeMap(sc,out, "BOLA-DRA")+ ggtitle('BOLA-DRA')+THEME
  C7 <- plotChangeMap(sc,out, "ITGAL")+ ggtitle('ITGAL')+THEME
  
  combined <- (C1+C3+C4) / (C5+C6+C7)
  combined + plot_layout(guides = "collect")
  ggsave(paste0('./SoupX/',SAMPLE_ID, '_post_soupX.jpeg'), width = 7, height = 5, units = 'in', bg = 'white')
  
  
  DropletUtils::write10xCounts(paste0("SoupX/", SAMPLE_ID), out, version = '3', overwrite = TRUE)
  return(out)
}



cellranger_dirs <- 
  list.dirs('cellranger_out', recursive = FALSE)


soupx_res <- 
  tibble(cellranger_dirs=paste0(cellranger_dirs,'/outs/')) %>% 
  mutate(soupx=map(.x=cellranger_dirs, .f=~run_soupx(.x)))



