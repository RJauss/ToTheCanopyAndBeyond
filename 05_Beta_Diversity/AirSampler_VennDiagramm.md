AirSampler VennDiagramm
================

Load Data
---------

``` r
rm(list = ls())

library(ggplot2)
library(viridis)
library(ggpubr)
library(plyr)
library(ggVennDiagram)

#setwd("04_Alpha_Diversity/")
OTU_Table = as.data.frame(read.csv("../00_Data/Oomycota/05_Oomycota_OTU_Table_new_min-freq-20617_min-feat-5_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))
SampleMetadata = OTU_Table[,1:5]
Stratum = SampleMetadata$Stratum
Season = SampleMetadata$Season
OTU_Table = OTU_Table[,6:ncol(OTU_Table)]
species = as.matrix(OTU_Table)

Season_aggregated = ddply(OTU_Table, "Season", numcolwise(sum))
rownames(Season_aggregated) = Season_aggregated$Season
Season_aggregated = Season_aggregated[,-1]
Season_aggregated = ifelse(Season_aggregated > 0, T, F)

Stratum_aggregated = ddply(OTU_Table, "Stratum", numcolwise(sum))
rownames(Stratum_aggregated) = Stratum_aggregated$Stratum
Stratum_aggregated = Stratum_aggregated[,-1]
Stratum_aggregated = ifelse(Stratum_aggregated > 0, T, F)

VennSeason = 
  list(March = colnames(Season_aggregated)[Season_aggregated["March", TRUE]], 
       May = colnames(Season_aggregated)[Season_aggregated["May", TRUE]])
VennStratum = 
  list(Canopy = colnames(Stratum_aggregated)[Stratum_aggregated["Canopy", TRUE]], 
       Ground = colnames(Stratum_aggregated)[Stratum_aggregated["Ground", TRUE]])
```

Cercozoa
--------

``` r
OTU_Table_cerco = as.data.frame(read.csv("../00_Data/Cercozoa/05_Cercozoa_OTU_Table_min-freq-16922_min-feat-5_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))
SampleMetadata_cerco = OTU_Table_cerco[,1:5]
Stratum_cerco = SampleMetadata_cerco$Stratum
Season_cerco = SampleMetadata_cerco$Season
OTU_Table_cerco = OTU_Table_cerco[,6:ncol(OTU_Table_cerco)]
species_cerco = as.matrix(OTU_Table_cerco)

Season_aggregated_cerco = ddply(OTU_Table_cerco, "Season_cerco", numcolwise(sum))
rownames(Season_aggregated_cerco) = Season_aggregated_cerco$Season_cerco
Season_aggregated_cerco = Season_aggregated_cerco[,-1]
Season_aggregated_cerco = ifelse(Season_aggregated_cerco > 0, T, F)

Stratum_aggregated_cerco = ddply(OTU_Table_cerco, "Stratum_cerco", numcolwise(sum))
rownames(Stratum_aggregated_cerco) = Stratum_aggregated_cerco$Stratum_cerco
Stratum_aggregated_cerco = Stratum_aggregated_cerco[,-1]
Stratum_aggregated_cerco = ifelse(Stratum_aggregated_cerco > 0, T, F)

VennSeason_cerco = 
  list(March = colnames(Season_aggregated_cerco)[Season_aggregated_cerco["March", TRUE]], 
       May = colnames(Season_aggregated_cerco)[Season_aggregated_cerco["May", TRUE]])
VennStratum_cerco = 
  list(Canopy = colnames(Stratum_aggregated_cerco)[Stratum_aggregated_cerco["Canopy", TRUE]], 
       Ground = colnames(Stratum_aggregated_cerco)[Stratum_aggregated_cerco["Ground", TRUE]])

g_Season = ggVennDiagram(VennSeason, show.legend = F) + 
  scale_fill_gradient2(limits=c(1,75), high = "darkslategray", low = "slategray2", mid = "slategray4", midpoint = 36)
```

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
g_Stratum = ggVennDiagram(VennStratum, show.legend = F) + 
  scale_fill_gradient2(limits=c(1,75), high = "darkslategray", low = "slategray2", mid = "slategray4", midpoint = 36)
```

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
g_Season_cerco = ggVennDiagram(VennSeason_cerco, show.legend = F) + 
  #scale_fill_gradient2(limits=c(8,88), high = "indianred4", low = "#fee5d9", mid = "#fb6a4a", midpoint = 41)
  scale_fill_viridis_b(begin = 0.3, end = 0.75, direction = -1,
                       option = "magma")
```

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
g_Stratum_cerco = ggVennDiagram(VennStratum_cerco, show.legend = F) +
  #scale_fill_gradient2(limits=c(8,88), high = "indianred4", low = "#fee5d9", mid = "#fb6a4a", midpoint = 41)
  scale_fill_viridis_b(begin = 0.3, end = 0.75, direction = -1, 
                       option = "magma")
```

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
g_Season$layers[[1]]$data$y = c(5, 5)
g_Stratum$layers[[1]]$data$y = c(5, 5)
g_Season_cerco$layers[[1]]$data$y = c(5, 5)
g_Stratum_cerco$layers[[1]]$data$y = c(5, 5)

g = ggarrange(g_Season, g_Stratum, 
              ncol = 1, nrow = 2, common.legend = T, 
              labels = NULL)
g_cerco = ggarrange(g_Season_cerco, g_Stratum_cerco, 
              ncol = 1, nrow = 2, common.legend = T, 
              labels = NULL)

g
```

![](AirSampler_VennDiagramm_files/figure-markdown_github/Cerco-1.png)

``` r
g_cerco
```

![](AirSampler_VennDiagramm_files/figure-markdown_github/Cerco-2.png)

Combine
-------

``` r
combi = ggarrange(g_cerco, g, 
                  ncol = 2, nrow = 1, common.legend = T, 
                  labels = "AUTO", vjust = 4)

ggsave("AirSampler_VennDiagramCombined.png", plot = combi, 
       device = "png", dpi = 300, width = 15, height = 15, 
       units = "cm")
ggsave("AirSampler_VennDiagramCombined.pdf", plot = combi, 
       device = "pdf", dpi = 300, width = 15, height = 15, 
       units = "cm")
ggsave("AirSampler_VennDiagramCombined.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 15, height = 15, 
       units = "cm")

combi
```

![](AirSampler_VennDiagramm_files/figure-markdown_github/CombineVenn-1.png)
