##main file
library("vegan")
library("FactoMineR")
library("ggplot2")
library("ggrepel")
library("ggiraph")
library("edgeR")
library('gsubfn')
library("factoextra")
library("utils")
library("BBmisc")
library("mefa4")
library("tibble")
library("DataCombine")
library("dplyr")
library("mixOmics")
library("apcluster") 
library("propr")
library("gap")
library("fdrtool")
library("EnhancedVolcano")
library("DESeq2")


source('differential_testing_deseq2.R')
source("download_summary_fcts.r")
source('nmds_phylum_genus_level_groupwise.R')
source('nmds_all_genus_level.R')
source('family_observations_exploratory.R')


## data download and initial summary tables
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


## NMDS + hierarchical clustering family observations genus level
raw_counts_family_matrix <- family_data_raw()[[1]]
family_data_tmm <- tmm_norm_family(raw_counts_family_matrix)

##nmds 
tmm_mds_family <- nmds_family(family_data_tmm)
cluster_tmm_family <- hc_family(family_data_tmm)

## 3 clusters division
clust_tmm_family3 <- cluster_division(cluster_tmm_family, 3)

scores.tmm.nmds_family3 <- nmds_hc_plts(tmm_mds_family, clust_tmm_family3, family_data_tmm)

## farm scores for 3 clusters
farm_scores3 <- farm_scores(scores.tmm.nmds_family3, 3)
print(farm_scores3)


## 2 clusters
clust_tmm_family2 <- cluster_division(cluster_tmm_family, 2)
scores.tmm.nmds_family2 <- nmds_hc_plts(tmm_mds_family, clust_tmm_family2, family_data_tmm)

## farm scores for 2 clusters
farm_scores2 <- farm_scores(scores.tmm.nmds_family2, 2)
print(farm_scores2)



### Affinity
raw_counts_family <- family_data_raw()[[3]]
affinity_results <- affinity_prop_clust(raw_counts_family)

## differential testing results

##DT of ap clusters
scores.tmm.nmds_family2 <- affinity_results[[2]]
cluster_division_ap = scores.tmm.nmds_family2$apclusters 
dt_deseq2_results_ap <- differential_testing(raw_counts_family, scores.tmm.nmds_family2, cluster_division_ap)
significant_genera_ap <- dt_deseq2_results_ap[[2]]
print(significant_genera_ap)

#p-value re-estimation with 'fdrtool'
res_ap = dt_deseq2_results_ap[[1]] 
pvalue_diagostics(res_ap, 'ap')

##volcano_plots ap
volcano_plots_ap(res_ap)

## DT of hierarchical clusters
cluster_division_hc = scores.tmm.nmds_family2$cluster
dt_deseq2_results_hc <- differential_testing(raw_counts_family, scores.tmm.nmds_family2, cluster_division_hc)
significant_genera_hc <- dt_deseq2_results_hc[[2]]
print(significant_genera_hc)

##p-value diagnostics of hc
res_hc = dt_deseq2_results_hc[[1]] 
res_hc_corr <- pvalue_diagostics(res_hc, 'hc')

##volcano plots hc
volcano_plots_hc(res_hc_corr)
