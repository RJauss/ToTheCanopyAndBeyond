# To the canopy and beyond: Air samples reveal wind dispersal as a driver of ubiquitous protistan pathogen assembly in tree canopies

Welcome to the **To the Canopy and Beyond** repository!

This repository is a collection of several scripts and mini-tutorials guiding you through the methods of metabarcoding analyses which were performed in the paper by [Jauss & Nowack et al., 2020](https://doi.org/10.1101/2020.11.30.405688). The data and methods are based on the paper by [Jauss & Walden et al., 2020](https://doi.org/10.22541/au.158679920.02842084) and the corresponding Github repository [FromForestSoilToCanopy](https://github.com/RJauss/FromForestSoilToCanopy). 

Most scripts here deal with the data for oomycetes, but they were also applied for Cercozoa if not stated otherwise. The raw data can be downloaded [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB37525), plots and figures were generated with the final OTU tables (not provided) and annotation files accessible in the folder [00_Data](00_Data/). 

## Table of Content
*Relevant scripts and tutorials can be found in the corresponding markdown files, which are* **linked and highlighted** *for each category.*

### 00 Data
This folder contains the final OTU table, taxonomic annotation, sample metadata as well as the oligosheet for demultiplexing the raw data. Intermediate files are not provided here, but you can generate them yourself by following the next steps.

### 01 Metabarcoding Pipeline
**[In this pipeline](01_Metabarcoding-Pipeline/Metabarcoding-Pipeline.md)**, you find the neccessary scripts to generate the (unfiltered) OTU table from raw .fastq files.

### 02 Taxonomic Annotation and Visualisation
This is a collection of several scripts neccessary for the taxonomic annotation of our OTUs. What we do first is to **[BLAST against the NCBI nt database to remove contaminants](02_Taxonomic_Annotation_and_Visualisation/BLAST-against-NCBI-nt-Database.md)**. After that, we **[download and process public oomycete/cercozoan sequences](02_Taxonomic_Annotation_and_Visualisation/Downloading-&-Processing-ITS-Sequences.md)**, which we use as a reference database for our **[taxonomic annotation with `vsearch`](02_Taxonomic_Annotation_and_Visualisation/Annotate-with-vsearch-and-the-ITS1-reference-database.md)**.

The visualisation of the taxonomy then includes a diagram showing the **[total taxonomic composition](02_Taxonomic_Annotation_and_Visualisation/AirSampler_TaxonomyOverview.md)**.
What's more, we also visualised the **[Taxonomy per sampling season in a sankey diagram](02_Taxonomic_Annotation_and_Visualisation/AirSampler_Sankey.md)**

### 03 Postprocessing the OTU Table
This section provides scripts on how to **[import, explore and filter the OTU table with `Qiime`](03_Postprocessing_OTU-Table/Importing-and-Filtering-OTU-Table.md)**, how to **[extract sequences from the filtered table](03_Postprocessing_OTU-Table/Postprocessing-the-OTU-Table.md#Extract-Sequences-from-Filtered-Table)** and last but not least how to **[paste the filtered metadata into the filtered table](03_Postprocessing_OTU-Table/Postprocessing-the-OTU-Table.md#Paste-Filtered-OTU-Table-and-Filtered-Metadata)**.

### 04 Determining Alpha Diversity
Here we deal with the methods of how to **[plot alpha diversity indices in a boxplot](04_Alpha_Diversity/AirSampler_AlphaBoxplot.md)** including pairwise significance testing. Also, we performed **[rarefaction analyses](04_Alpha_Diversity/AirSampler_RarefactionCurves.md)** to test if a higher sequencing depth would have revealed more OTUs.

### 05 Exploring Beta Diversity
One of the most straightforward methods of visualising beta diversity is an NMDS plot, the script is provided **[here](05_Beta_Diversity/AirSampler_NMDS.md)**. But we also plot the number of **[shared OTUs in a Venn diagram](05_Beta_Diversity/AirSampler_VennDiagramm.md)**. 

