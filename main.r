##main file
library("vegan")
library("FactoMineR")
library("ggplot2")
library("ggrepel")
library("ggiraph")
library("edgeR")
library('gsubfn')
source('nmds_all_genus_level.R')

## data download and initial summary tables
source("download_summary_fcts.r")
meta_data = meta_data_download()
sum_asv_count_genus = asv_data_download()
overdisp_plot = overdisp_plot(sum_asv_count_genus)
#print(overdisp_plot)
summary_table = summary_table(sum_asv_count_genus)
#print(summary_table)


## 4 Results of Analysis 
### 4.1 Explorative Analysis of Piglet Microbiome Data
#### 4.1.1 Nonmetric Multidimensional Scaling (NMDS)

data_microbiome <- microbiome_data()

ind_S14B = as.logical(meta_data$Animal == "S" & meta_data$Time_sampling == "-14")
ind_S7B = as.logical(meta_data$Animal == "S" & meta_data$Time_sampling == "7")
ind_F7B = as.logical(meta_data$Animal == "F" & meta_data$Time_sampling == "7"& meta_data$Time_points == "B")
ind_F7W = as.logical(meta_data$Animal == "F" & meta_data$Time_sampling == "7"& meta_data$Time_points == "W")

## NMDS on phylum level
source('nmds_phylum_level.R')
nmds_phylum(data_microbiome, ind_S14B)
nmds_phylum(data_microbiome, ind_S7B)
nmds_phylum(data_microbiome, ind_F7B)
nmds_phylum(data_microbiome, ind_F7W)


## NMDS on genus level
genus_data <- genus_raw_and_tss_data()
genus_tss_data <- genus_data[[2]]
print(genus_data)
nmds_genus(genus_tss_data, ind_S14B, 'ind_S14B')
nmds_genus(genus_tss_data, ind_S7B, 'ind_S7B')
nmds_genus(genus_tss_data, ind_F7B, 'ind_F7B')
nmds_genus(genus_tss_data, ind_F7W, 'ind_F7W')

## NMDS of all 4 groups + hierarchical clustering (genus level)
genus_raw <- genus_data[[1]]
genus_tmm = tmm_norm(genus_raw)
nmds_hc_four_groups(genus_tmm, ind_S14B, ind_S7B, ind_F7B, ind_F7W)
