AirSampler AlphaBoxplot
================

## Load Data

``` r
rm(list = ls())

library(ggplot2)
library(viridis)
library(vegan)
library(ggpubr)
library(ggsignif)

#setwd("04_Alpha_Diversity/")
OTU_Table = as.data.frame(read.csv("../00_Data/Oomycota/05_Oomycota_OTU_Table_new_min-freq-20617_min-feat-5_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))
SampleMetadata = OTU_Table[,1:5]
Stratum = SampleMetadata$Stratum
Timepoint = SampleMetadata$Timepoint
OTU_Table = OTU_Table[,6:ncol(OTU_Table)]
```

## Calculate Alpha Diversity

To run the diversity analyses, simply load the table and specify the
`index` - in this case: The Simpson index. Then convert it into a
dataframe and add the metadata and group.

``` r
shannon = diversity(OTU_Table, index = "shannon")
shannon = as.data.frame(shannon)
shannon$Comparison = Timepoint
rownames(shannon) = SampleMetadata$SampleID
shannon$Group = "Oomycota"
df = shannon
df$richness = specnumber(OTU_Table)
df$evenness = df$shannon/log(df$richness)
df$Type = "Timepoint"
df2 = df
df2$Comparison = Stratum
df2$Type = "Stratum"
df3 = rbind(df, df2)
df_melted = reshape2::melt(df3)
```

    ## Using Comparison, Group, Type as id variables

## Plot the Figure

Now we put the diverity measurements into a habitat specific context. It
can be easiest visualised in a boxplot:

``` r
g = ggplot(df_melted, aes(x = Type, y = value, fill = Comparison)) + 
  #stat_boxplot(geom = "errorbar", width = 0.2, show.legend = F) +
  geom_boxplot(show.legend = T) + 
  theme_minimal() + 
  labs(y = "Alpha diversity", 
       x = NULL) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        strip.text = element_text(size=14, face = "bold")) +
  #geom_signif(comparisons = list(c("March", "May")), 
  #            map_signif_level=TRUE) + 
  #geom_signif(comparisons = list(c("Canopy", "Ground")), 
  #            map_signif_level=TRUE) +
  stat_compare_means(aes(group = Comparison), label = "p.signif", 
                     size = 5, method = "wilcox.test", vjust = 0.7, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns"))) +
  scale_fill_manual(values = c("darkslategrey", "firebrick", 
                               "darkolivegreen4", "burlywood3"), 
                    limits = c("March", "May", "Canopy", "Ground"), 
                    name = "Timepoint/\nStratum") +
  facet_grid(rows = vars(variable), scales = "free", switch = "y")
```

## Cercozoa

Let’s check if the Cercozoa show a similar pattern:

``` r
Cerco_OTU_Table = as.data.frame(read.csv("../00_Data/Cercozoa/05_Cercozoa_OTU_Table_min-freq-16922_min-feat-5_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))
Cerco_SampleMetadata = Cerco_OTU_Table[,1:5]
Cerco_Timepoint = Cerco_SampleMetadata$Timepoint
Cerco_Stratum = Cerco_SampleMetadata$Stratum
Cerco_OTU_Table = Cerco_OTU_Table[,6:ncol(Cerco_OTU_Table)]

Cerco_shannon = diversity(Cerco_OTU_Table, index = "shannon")
Cerco_shannon = as.data.frame(Cerco_shannon)
rownames(Cerco_shannon) = Cerco_SampleMetadata$SampleID
colnames(Cerco_shannon) = "shannon"
Cerco_shannon$Comparison = Cerco_Timepoint
Cerco_shannon$Group = "Cercozoa"
df_cerco = Cerco_shannon
df_cerco$richness = specnumber(Cerco_OTU_Table)
df_cerco$evenness = df_cerco$shannon/log(df_cerco$richness)
df_cerco$Type = "Timepoint"
df_cerco2 = df_cerco
df_cerco2$Comparison = Cerco_Stratum
df_cerco2$Type = "Stratum"
df_cerco_3 = rbind(df_cerco, df_cerco2)

df_cerco_melted = reshape2::melt(df_cerco_3)
```

## Plot Cercozoa Figure

``` r
g_cerco = ggplot(df_cerco_melted, aes(x = Type, y = value, fill = Comparison)) + 
  #stat_boxplot(geom = "errorbar", width = 0.2, show.legend = F) +
  geom_boxplot(show.legend = T) + 
  theme_minimal() + 
  labs(y = "Alpha diversity", 
       x = NULL) +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"), 
        strip.text = element_text(size=14, face = "bold")) +
  #geom_signif(comparisons = list(c("March", "May")), 
  #            map_signif_level=TRUE) + 
  #geom_signif(comparisons = list(c("Canopy", "Ground")), 
  #            map_signif_level=TRUE) + 
  stat_compare_means(aes(group = Comparison), label = "p.signif", 
                     size = 5, method = "wilcox.test", vjust = 0.7, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns"))) +
  scale_fill_manual(values = c("darkslategrey", "firebrick", 
                               "darkolivegreen4", "burlywood3"), 
                    limits = c("March", "May", "Canopy", "Ground"), 
                    name = "Timepoint/\nStratum") +
  facet_grid(rows = vars(variable), scales = "free", switch = "y")
```

## Combine Plots

``` r
combi = ggarrange(g_cerco, g, 
                  labels = c("A", "B"), 
                  ncol = 2, nrow = 1, 
                  align = "v", common.legend = T, 
                  legend = "right")

combi = combi + theme(rect = element_rect(fill = "transparent"))#%>%
  #annotate_figure(top = text_grob("Alpha diversity of samples", 
  #                                face = "bold", size = 20), 
  #                fig.lab = "Figure X", fig.lab.face = "bold", 
  #                fig.lab.size = 18)
ggsave("AlphaBoxplotCombined.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 17.7, height = 10, 
       units = "cm")
ggsave("AlphaBoxplotCombined.png", plot = combi, 
       device = "png", dpi = 300, width = 17.7, height = 10, 
       units = "cm")
ggsave("AlphaBoxplotCombined.pdf", plot = combi, 
       device = "pdf", dpi = 300, width = 17.7, height = 10, 
       units = "cm")
ggsave("AlphaBoxplotCombined.tiff", plot = combi, 
       device = "tiff", dpi = 300, width = 17.7, height = 10, 
       units = "cm", compression = "lzw")
combi
```

![](AirSampler_AlphaBoxplot_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Account for multiple testing

The Wilcoxon Rank Sign test can be subject to multiple testing, thus we
use the general linear models in the next step for significance testing.

``` r
shannon = diversity(OTU_Table, index = "shannon")
shannon = as.data.frame(shannon)
shannon$Timepoint = Timepoint
shannon$Stratum = Stratum
rownames(shannon) = SampleMetadata$SampleID

Cerco_shannon = diversity(Cerco_OTU_Table, index = "shannon")
Cerco_shannon = as.data.frame(Cerco_shannon)
Cerco_shannon$Timepoint = Cerco_Timepoint
Cerco_shannon$Stratum = Cerco_Stratum
rownames(Cerco_shannon) = Cerco_SampleMetadata$SampleID

summary(glm(data = shannon, shannon ~ Timepoint+Stratum))
```

    ## 
    ## Call:
    ## glm(formula = shannon ~ Timepoint + Stratum, data = shannon)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.91493  -0.25953  -0.03547   0.12880   0.91657  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     0.1594     0.1254   1.271    0.213    
    ## TimepointMay    1.0070     0.1394   7.222 4.87e-08 ***
    ## StratumGround   0.1410     0.1363   1.034    0.309    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.1529632)
    ## 
    ##     Null deviance: 12.8287  on 32  degrees of freedom
    ## Residual deviance:  4.5889  on 30  degrees of freedom
    ## AIC: 36.545
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(data = Cerco_shannon, Cerco_shannon ~ Timepoint+Stratum))
```

    ## 
    ## Call:
    ## glm(formula = Cerco_shannon ~ Timepoint + Stratum, data = Cerco_shannon)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.90864  -0.21865  -0.09351   0.26620   0.86709  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     1.4173     0.2604   5.443 0.000284 ***
    ## TimepointMay    1.2357     0.3007   4.110 0.002111 ** 
    ## StratumGround  -0.4061     0.3007  -1.350 0.206642    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.2486267)
    ## 
    ##     Null deviance: 6.9400  on 12  degrees of freedom
    ## Residual deviance: 2.4863  on 10  degrees of freedom
    ## AIC: 23.388
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
richness = specnumber(OTU_Table)
richness = as.data.frame(richness)
richness$Timepoint = Timepoint
richness$Stratum = Stratum

Cerco_richness = specnumber(Cerco_OTU_Table)
Cerco_richness = as.data.frame(Cerco_richness)
Cerco_richness$Timepoint = Cerco_Timepoint
Cerco_richness$Stratum = Cerco_Stratum

summary(glm(data = richness, richness ~ Timepoint+Stratum))
```

    ## 
    ## Call:
    ## glm(formula = richness ~ Timepoint + Stratum, data = richness)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -12.843   -3.357    1.224    4.157   12.643  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     10.776      2.005   5.375 8.08e-06 ***
    ## TimepointMay    19.581      2.229   8.784 8.55e-10 ***
    ## StratumGround    3.486      2.179   1.599     0.12    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 39.0952)
    ## 
    ##     Null deviance: 4335.3  on 32  degrees of freedom
    ## Residual deviance: 1172.9  on 30  degrees of freedom
    ## AIC: 219.48
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(data = Cerco_richness, Cerco_richness ~ Timepoint+Stratum))
```

    ## 
    ## Call:
    ## glm(formula = Cerco_richness ~ Timepoint + Stratum, data = Cerco_richness)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -17.9394   -4.6364   -0.4545    3.8788   17.0606  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     34.455      6.273   5.492 0.000265 ***
    ## TimepointMay     9.182      7.244   1.268 0.233668    
    ## StratumGround   -3.515      7.244  -0.485 0.637935    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 144.2909)
    ## 
    ##     Null deviance: 1695.7  on 12  degrees of freedom
    ## Residual deviance: 1442.9  on 10  degrees of freedom
    ## AIC: 106.12
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
evenness = shannon$shannon/log(richness$richness)
evenness = as.data.frame(evenness)
evenness$Timepoint = Timepoint
evenness$Stratum = Stratum

Cerco_evenness = Cerco_shannon$Cerco_shannon/log(Cerco_richness$Cerco_richness)
Cerco_evenness = as.data.frame(Cerco_evenness)
Cerco_evenness$Timepoint = Cerco_Timepoint
Cerco_evenness$Stratum = Cerco_Stratum

summary(glm(data = evenness, evenness ~ Timepoint+Stratum))
```

    ## 
    ## Call:
    ## glm(formula = evenness ~ Timepoint + Stratum, data = evenness)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.26475  -0.08705  -0.02202   0.06808   0.31642  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    0.08091    0.04019   2.013   0.0531 .  
    ## TimepointMay   0.26780    0.04468   5.993 1.42e-06 ***
    ## StratumGround  0.02052    0.04369   0.470   0.6420    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.01570861)
    ## 
    ##     Null deviance: 1.04309  on 32  degrees of freedom
    ## Residual deviance: 0.47126  on 30  degrees of freedom
    ## AIC: -38.562
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(data = Cerco_evenness, Cerco_evenness ~ Timepoint+Stratum))
```

    ## 
    ## Call:
    ## glm(formula = Cerco_evenness ~ Timepoint + Stratum, data = Cerco_evenness)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.26439  -0.08529  -0.04488   0.09716   0.33353  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)    0.39758    0.08885   4.475  0.00119 **
    ## TimepointMay   0.31828    0.10260   3.102  0.01121 * 
    ## StratumGround -0.09320    0.10260  -0.908  0.38504   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.02894647)
    ## 
    ##     Null deviance: 0.58036  on 12  degrees of freedom
    ## Residual deviance: 0.28946  on 10  degrees of freedom
    ## AIC: -4.5683
    ## 
    ## Number of Fisher Scoring iterations: 2

These results are however the same as in the boxplot with the Wilcoxon
Test, we therefore do not need to update the figure per se.
