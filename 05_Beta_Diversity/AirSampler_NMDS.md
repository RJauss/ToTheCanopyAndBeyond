AirSampler\_NMDS
================

Load Data
---------

``` r
rm(list = ls())

library(vegan)
library(ggplot2)
library(funrar)
library(ggpubr)


OTU_Table = read.csv("../00_Data/Oomycota/05_Oomycota_OTU_Table_new_min-freq-20617_min-feat-5_transposed_withMetadata.tsv", 
                     header=T, stringsAsFactors=TRUE, sep="\t")
species = OTU_Table[,6:ncol(OTU_Table)]
species_mat = as.matrix(species)
species_mat = make_relative(species_mat)
species_mat = log(species_mat +1)
SampleMetadata = as.data.frame(OTU_Table[,1:5])
```

Preparing the data
------------------

``` r
Dist = vegdist(species_mat, 
               diag = T, 
               na.rm = T)

OTU.NMDS.bray = metaMDS(Dist, # The distance matrix
                        distance = "bray", 
                        k=3, # How many dimensions/axes to display 
                        trymax=100, # max number of iterations
                        wascores=TRUE, # Weighted species scores
                        trace=TRUE, 
                        zerodist="add") # What to do with 0's after sqrt transformation

# Get the sample scores and add the metadata
data.scores = as.data.frame(scores(OTU.NMDS.bray))
data.scores$site = rownames(data.scores)
data.scores$Season = SampleMetadata$Season
data.scores$Stratum = SampleMetadata$Stratum
data.scores$TreeSpecies = SampleMetadata$TreeSpecies
data.scores$SampleID = SampleMetadata$X.SampleID

# Group the samples by metadata (in this case, microhabitat and tree species)

Group.Tilia = data.scores[data.scores$TreeSpecies == "Lime",][chull(data.scores[data.scores$TreeSpecies == "Lime", c("NMDS1", "NMDS2")]), ]
Group.Quercus = data.scores[data.scores$TreeSpecies == "Oak",][chull(data.scores[data.scores$TreeSpecies == "Oak", c("NMDS1", "NMDS2")]), ]
Group.Fraxinus = data.scores[data.scores$TreeSpecies == "Ash",][chull(data.scores[data.scores$TreeSpecies == "Ash", c("NMDS1", "NMDS2")]), ]
Group.Crane = data.scores[data.scores$TreeSpecies == "Crane",][chull(data.scores[data.scores$TreeSpecies == "Crane", c("NMDS1", "NMDS2")]), ]

# The hull-data will be needed by ggplot later to draw the polygons
Group.May = data.scores[data.scores$Season == "May",][chull(data.scores[data.scores$Season == "May", c("NMDS1", "NMDS2")]), ]
Group.March = data.scores[data.scores$Season == "March",][chull(data.scores[data.scores$Season == "March", c("NMDS1", "NMDS2")]), ]

Group.Canopy = data.scores[data.scores$Stratum == "Canopy",][chull(data.scores[data.scores$Stratum == "Canopy", c("NMDS1", "NMDS2")]), ]
Group.Ground = data.scores[data.scores$Stratum == "Ground",][chull(data.scores[data.scores$Stratum == "Ground", c("NMDS1", "NMDS2")]), ]

# The hull-data will be needed by ggplot later to draw the polygons

hull.data_TreeSpecies = rbind(Group.Tilia, Group.Quercus, Group.Fraxinus, Group.Crane)

hull.data_Season = rbind(Group.May, Group.March)

hull.data_Stratum = rbind(Group.Canopy, Group.Ground)
```

Plot NMDS
---------

``` r
g = ggplot() + 
  geom_polygon(data = hull.data_Stratum, 
               aes(x=NMDS1, y=NMDS2, group = Stratum, fill = Stratum), 
               alpha = 0.7, color = NA, linetype = "solid") +
  scale_fill_manual(values = c("darkolivegreen4", "burlywood3")) +
  geom_point(data = data.scores, 
             aes(x = NMDS1, y = NMDS2), 
             size = 3,
             color = "#5d5f66") + 
  geom_polygon(data = hull.data_Season, 
               aes(x=NMDS1, y=NMDS2, group = Season, color = Season), 
               alpha = 0.7, fill = NA, linetype = "dashed") +
  scale_color_manual(values = c("darkslategrey", "firebrick")) +
  geom_text(aes(x = -0.2, y = -0.5, label = as.character(paste0(OTU.NMDS.bray$ndim, "D Stress: ", round(as.numeric(OTU.NMDS.bray$stress), digits = 3)))), parse = F, color = "#5d5f66", size = 4) +
  theme_minimal() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14, face = "bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

g
```

![](AirSampler_NMDS_files/figure-markdown_github/ggplot-1.png)

Load Cerco Data
---------------

``` r
OTU_Table_cerco = read.csv("../00_Data/Cercozoa/05_Cercozoa_OTU_Table_min-freq-16922_min-feat-5_transposed_withMetadata.tsv", 
                     header=T, stringsAsFactors=TRUE, sep="\t")
species_cerco = OTU_Table_cerco[,6:ncol(OTU_Table_cerco)]
species_mat_cerco = as.matrix(species_cerco)
species_mat_cerco = make_relative(species_mat_cerco)
species_mat_cerco = log(species_mat_cerco +1)
SampleMetadata_cerco = as.data.frame(OTU_Table_cerco[,1:5])
```

Prepare Cerco Data
------------------

``` r
Dist_cerco = vegdist(species_mat_cerco, 
               diag = T, 
               na.rm = T)

OTU.NMDS.bray_cerco = metaMDS(species_mat_cerco, # The distance matrix
                        distance = "bray", 
                        k=3, # How many dimensions/axes to display 
                        trymax=100, # max number of iterations
                        wascores=TRUE, # Weighted species scores
                        trace=TRUE, 
                        zerodist="add") # What to do with 0's after sqrt transformation
```

    ## Run 0 stress 0.110583 
    ## Run 1 stress 0.109516 
    ## ... New best solution
    ## ... Procrustes: rmse 0.08763993  max resid 0.1833577 
    ## Run 2 stress 0.1095137 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005605168  max resid 0.001174718 
    ## ... Similar to previous best
    ## Run 3 stress 0.1110943 
    ## Run 4 stress 0.129178 
    ## Run 5 stress 0.1399631 
    ## Run 6 stress 0.1128075 
    ## Run 7 stress 0.1312433 
    ## Run 8 stress 0.1094977 
    ## ... New best solution
    ## ... Procrustes: rmse 0.003874541  max resid 0.009444079 
    ## ... Similar to previous best
    ## Run 9 stress 0.1139678 
    ## Run 10 stress 0.1095133 
    ## ... Procrustes: rmse 0.003787036  max resid 0.009124921 
    ## ... Similar to previous best
    ## Run 11 stress 0.1094951 
    ## ... New best solution
    ## ... Procrustes: rmse 0.007004104  max resid 0.01557613 
    ## Run 12 stress 0.1128067 
    ## Run 13 stress 0.1127909 
    ## Run 14 stress 0.1105789 
    ## Run 15 stress 0.1139999 
    ## Run 16 stress 0.1128 
    ## Run 17 stress 0.1095154 
    ## ... Procrustes: rmse 0.01119325  max resid 0.02560189 
    ## Run 18 stress 0.1095099 
    ## ... Procrustes: rmse 0.009948555  max resid 0.02323472 
    ## Run 19 stress 0.1235742 
    ## Run 20 stress 0.1128195 
    ## Run 21 stress 0.1128095 
    ## Run 22 stress 0.10952 
    ## ... Procrustes: rmse 0.01205909  max resid 0.0271065 
    ## Run 23 stress 0.1105859 
    ## Run 24 stress 0.1128042 
    ## Run 25 stress 0.1235639 
    ## Run 26 stress 0.1094925 
    ## ... New best solution
    ## ... Procrustes: rmse 0.004912314  max resid 0.01132334 
    ## Run 27 stress 0.1095138 
    ## ... Procrustes: rmse 0.005883489  max resid 0.01362403 
    ## Run 28 stress 0.129389 
    ## Run 29 stress 0.1293897 
    ## Run 30 stress 0.1218225 
    ## Run 31 stress 0.1193179 
    ## Run 32 stress 0.1128137 
    ## Run 33 stress 0.1095053 
    ## ... Procrustes: rmse 0.003988854  max resid 0.009055975 
    ## ... Similar to previous best
    ## *** Solution reached

    ## Warning in postMDS(out$points, dis, plot = max(0, plot - 1), ...): skipping
    ## half-change scaling: too few points below threshold

``` r
# Get the sample scores and add the metadata
data.scores_cerco = as.data.frame(scores(OTU.NMDS.bray_cerco))
data.scores_cerco$site = rownames(data.scores_cerco)
data.scores_cerco$Season = SampleMetadata_cerco$Season
data.scores_cerco$Stratum = SampleMetadata_cerco$Stratum
data.scores_cerco$TreeSpecies = SampleMetadata_cerco$TreeSpecies
data.scores_cerco$SampleID = SampleMetadata_cerco$X.SampleID

# Group_cerco the samples by metadata (in this case, microhabitat and tree species)

Group_cerco.Tilia = data.scores_cerco[data.scores_cerco$TreeSpecies == "Lime",][chull(data.scores_cerco[data.scores_cerco$TreeSpecies == "Lime", c("NMDS1", "NMDS2")]), ]
Group_cerco.Quercus = data.scores_cerco[data.scores_cerco$TreeSpecies == "Oak",][chull(data.scores_cerco[data.scores_cerco$TreeSpecies == "Oak", c("NMDS1", "NMDS2")]), ]
Group_cerco.Fraxinus = data.scores_cerco[data.scores_cerco$TreeSpecies == "Ash",][chull(data.scores_cerco[data.scores_cerco$TreeSpecies == "Ash", c("NMDS1", "NMDS2")]), ]
Group_cerco.Crane = data.scores_cerco[data.scores_cerco$TreeSpecies == "Crane",][chull(data.scores_cerco[data.scores_cerco$TreeSpecies == "Crane", c("NMDS1", "NMDS2")]), ]

# The hull-data will be needed by ggplot later to draw the polygons
Group_cerco.May = data.scores_cerco[data.scores_cerco$Season == "May",][chull(data.scores_cerco[data.scores_cerco$Season == "May", c("NMDS1", "NMDS2")]), ]
Group_cerco.March = data.scores_cerco[data.scores_cerco$Season == "March",][chull(data.scores_cerco[data.scores_cerco$Season == "March", c("NMDS1", "NMDS2")]), ]

Group_cerco.Canopy = data.scores_cerco[data.scores_cerco$Stratum == "Canopy",][chull(data.scores_cerco[data.scores_cerco$Stratum == "Canopy", c("NMDS1", "NMDS2")]), ]
Group_cerco.Ground = data.scores_cerco[data.scores_cerco$Stratum == "Ground",][chull(data.scores_cerco[data.scores_cerco$Stratum == "Ground", c("NMDS1", "NMDS2")]), ]

# The hull-data will be needed by ggplot later to draw the polygons

hull.data_cerco_TreeSpecies = rbind(Group_cerco.Tilia, Group_cerco.Quercus, Group_cerco.Fraxinus, Group_cerco.Crane)

hull.data_cerco_Season = rbind(Group_cerco.May, Group_cerco.March)

hull.data_cerco_Stratum = rbind(Group_cerco.Canopy, Group_cerco.Ground)
```

Plot Cerco Data
---------------

``` r
g_cerco = ggplot() + 
  geom_polygon(data = hull.data_cerco_Stratum, 
               aes(x=NMDS1, y=NMDS2, group = Stratum, fill = Stratum), 
               alpha = 0.7, color = NA, linetype = "solid") +
  scale_fill_manual(values = c("darkolivegreen4", "burlywood3")) +
  geom_point(data = data.scores_cerco, 
             aes(x = NMDS1, y = NMDS2), 
             size = 3,
             color = "#5d5f66") + 
  geom_polygon(data = hull.data_cerco_Season, 
               aes(x=NMDS1, y=NMDS2, group = Season, color = Season), 
               alpha = 0.7, fill = NA, linetype = "dashed") +
  scale_color_manual(values = c("darkslategrey", "firebrick")) +
  geom_text(aes(x = -0.2, y = -0.5, 
                label = as.character(paste0(OTU.NMDS.bray_cerco$ndim, "D Stress: ", round(as.numeric(OTU.NMDS.bray_cerco$stress), digits = 3)))), parse = F, color = "#5d5f66", size = 4) +
  theme_minimal() + 
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14, face = "bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

g_cerco
```

![](AirSampler_NMDS_files/figure-markdown_github/NMDS%20Cerco-1.png)

Combine Plots
-------------

``` r
g$labels$title = NULL
g_cerco$labels$title = NULL
combi = ggarrange(g_cerco, g, 
                  labels = c("A", "B"), 
                  ncol = 2, nrow = 1, 
                  common.legend = T, legend = "right", 
                  align = "h", vjust = 1.5) #%>%
  #annotate_figure(fig.lab = "Figure X", fig.lab.face = "bold", 
  #                fig.lab.size = 18, 
  #                top = text_grob("Non-metric multidimensional scaling", 
  #                                face = "bold", size = 20))
#ggsave("NMDSCombined.tif", plot = combi, 
#       device = "tiff", dpi = 600, width = 28, height = 11, 
#       units = "cm")
ggsave("AirSampler_NMDSCombined.png", plot = combi, 
       device = "png", dpi = 300, width = 20, height = 10, 
       units = "cm")
ggsave("AirSampler_NMDSCombined.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 20, height = 10, 
       units = "cm")
ggsave("AirSampler_NMDSCombined.pdf", plot = combi, 
       device = "pdf", dpi = 300, width = 20, height = 10, 
       units = "cm")

combi
```

![](AirSampler_NMDS_files/figure-markdown_github/CombineNMDS-1.png)
