tmm_norm_family <-function(raw_counts_family_matrix){
  ## tmm normalization
  dge_family = DGEList(raw_counts_family_matrix)
  dge_family <- calcNormFactors(dge_family, method = "TMM")
  #dge2$samples
  
  pseudo_TMM_family <- log2(cpm(dge_family) + 1)
  empty_rows = apply(pseudo_TMM_family, 1, function(x) sum(x)==0)
  pseudo_TMM_family <- pseudo_TMM_family[!empty_rows,]
  
  pseudo_TMM_family_nmds = as.data.frame(t(pseudo_TMM_family))
  pseudo_TMM_family_nmds$farm = as.factor(substr(row.names(pseudo_TMM_family_nmds), 1, 1))
  return(pseudo_TMM_family_nmds)
}


nmds_family <- function(family_data_tmm){
  set.seed(23)
  tmm_mds_family = metaMDS(comm = family_data_tmm[,-1169], k=2, distance = "bray", trymax = 100) 
  return(tmm_mds_family)
}
  
  
hc_family <- function(family_data_tmm){
  bc_dist_family <- vegdist(family_data_tmm[,-1169], method = "bray")
  cluster_tmm_family <- hclust(bc_dist_family, method = 'ward.D2')
  #plot(cluster_tmm_family) 
  
  fviz_nbclust(family_data_tmm[,-1169], FUNcluster=hcut, diss=bc_dist_family, method = "silhouette")
  #ggsave("optimal_n_clusters_genus_family.eps", device = "eps", width=9, height=5.5)
  return(cluster_tmm_family)
}  


cluster_division <- function(cluster_tmm_family, k){
  clust_tmm_family <- cutree(cluster_tmm_family, k = k)
  clust_tmm_family = as.data.frame(clust_tmm_family)
  return(clust_tmm_family)
}

  

nmds_hc_plts <-function(tmm_mds_family, clust_tmm_family, family_data_tmm){
  ## plots inside, returns df with nmds scores
  
  ## nmds scores for plots
  scores.tmm.nmds_family = as.data.frame(scores(tmm_mds_family))
  scores.tmm.nmds_family$site = rownames(scores.tmm.nmds_family)
  scores.tmm.nmds_family$farm = as.factor(family_data_tmm$farm)
  scores.tmm.nmds_family$cluster = clust_tmm_family$clust_tmm_family
  scores.tmm.nmds_family$cluster = as.factor(scores.tmm.nmds_family$cluster)
  
  genus.tmm.nmds_family = as.data.frame(scores(tmm_mds_family, "species"))
  genus.tmm.nmds_family$genus = rownames(genus.tmm.nmds_family)
  
  gg_cluster = ggplot() + 
    geom_text(data = scores.tmm.nmds_family, aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.025), size = 2) + 
    geom_point(data = scores.tmm.nmds_family, aes(x = NMDS1, y = NMDS2, colour = cluster, fill = cluster), size = 2, shape = 21) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_equal() +
    theme_bw()
  
  gg_cluster
  #ggsave("nmds_hc_genus_family.eps", device = "eps", width=9, height=6)
  
  
  gg_cluster_zoom = ggplot() + 
    geom_text(data = subset(scores.tmm.nmds_family, scores.tmm.nmds_family$NMDS1 < 0.4 & scores.tmm.nmds_family$NMDS1 > -0.4 & scores.tmm.nmds_family$NMDS2 < 0.25 & scores.tmm.nmds_family$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.013), size = 2.4) + 
    geom_point(data = subset(scores.tmm.nmds_family, scores.tmm.nmds_family$NMDS1 < 0.4 & scores.tmm.nmds_family$NMDS1 > -0.4 & scores.tmm.nmds_family$NMDS2 < 0.25 & scores.tmm.nmds_family$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, colour = cluster, fill = cluster), size = 2, shape = 21) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_cartesian(ylim=c(- 0.25, 0.25), xlim = c(- 0.4, 0.4)) +
    theme_bw()
  gg_cluster_zoom
  #ggsave("nmds_hc_genus_family_zoom.eps", device = "eps", width=9, height=6)
  
  gg_farm = ggplot() + 
    geom_text(data = scores.tmm.nmds_family, aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.025), size = 2) + 
    geom_point(data = scores.tmm.nmds_family, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, shape = cluster), size = 2) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_equal() +
    theme_bw()
  gg_farm
  #ggsave("nmds_hc_genus_family_farmcolour.eps", device = "eps", width=9, height=6)
  
  gg_farm_zoom = ggplot() + 
    geom_text(data = subset(scores.tmm.nmds_family, scores.tmm.nmds_family$NMDS1 < 0.4 & scores.tmm.nmds_family$NMDS1 > -0.4 & scores.tmm.nmds_family$NMDS2 < 0.25 & scores.tmm.nmds_family$NMDS1 > -0.3), aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.013), size = 2.4) + 
    geom_point(data = subset(scores.tmm.nmds_family, scores.tmm.nmds_family$NMDS1 < 0.4 & scores.tmm.nmds_family$NMDS1 > -0.4 & scores.tmm.nmds_family$NMDS2 < 0.25 & scores.tmm.nmds_family$NMDS1 > -0.3), aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, shape = cluster), size = 2) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_cartesian(ylim=c(- 0.3, 0.25), xlim = c(- 0.4, 0.4))  +
    theme_bw()
  gg_farm_zoom
  
  return(scores.tmm.nmds_family)
}  
  

##calculating farm scores for each cluster 
farm_scores <- function(scores.tmm.nmds_family, k) {
  cluster1_family = scores.tmm.nmds_family[which(scores.tmm.nmds_family$cluster==1),]
  cluster2_family = scores.tmm.nmds_family[which(scores.tmm.nmds_family$cluster==2),]
  
  scores1_family = as.data.frame(summary(cluster1_family$farm))
  scores1_family$scores1 = prop.table(scores1_family$`summary(cluster1_family$farm)`)
  #scores1_family$`summary(cluster1_family$farm)`
  colnames(scores1_family) = c("count", "scores")
  
  scores2_family = as.data.frame(summary(cluster2_family$farm))
  scores2_family$scores2 = prop.table(scores2_family$`summary(cluster2_family$farm)`)
  #scores2_family$`summary(cluster2_family$farm)`
  colnames(scores2_family) = c("count", "scores")
  
  if (k==3){
    cluster3_family = scores.tmm.nmds_family[which(scores.tmm.nmds_family$cluster==3),]
    scores3_family = as.data.frame(summary(cluster3_family$farm))
    scores3_family$scores3 = prop.table(scores3_family$`summary(cluster3_family$farm)`)
    #scores3_family$`summary(cluster3_family$farm)`
    colnames(scores3_family) = c("count", "scores")
    farm_scores_family <- data.frame("farm" = levels(cluster1_family$farm), 
                                     "Scores_cluster1" = round(scores1_family$scores, 2), 
                                     "Scores_cluster2" = round(scores2_family$scores, 2), 
                                     "Scores_cluster3" = round(scores3_family$scores, 2))
    return(farm_scores_family) 
  } else{
    farm_scores_family <- data.frame("farm" = levels(cluster1_family$farm), 
                                     "Scores_cluster1" = round(scores1_family$scores, 2), 
                                     "Scores_cluster2" = round(scores2_family$scores, 2))
    return(farm_scores_family) 
  } 
}




affinity_prop_clust <- function(raw_counts_family){
  #rownames(raw_counts_family) = with(raw_counts_family, paste0(farm, individual))
  #raw_counts_family = raw_counts_family[,-c(1,2)]
  raw_counts_family_pseudo = raw_counts_family+1
  
  
  ## log-ratio transformation
  clr_trans_data = propr(as.matrix(raw_counts_family_pseudo))
  clr_data = clr_trans_data@logratio
  
  ## filter those genera, which sum count less than 5
  genus_counts <- apply(raw_counts_family, 2, sum)
  ind_ok_counts = as.logical(genus_counts > 5) ## 566 filtered
  #keep <- rowSums(cpm(y)>5) >= 2
  raw_counts_family_filtered = raw_counts_family[,ind_ok_counts]
  
  
  # AP clustering itself
  ## first Bray-Curtis dissimilarity is calculated separately
  ### BC dissimilarity on full data
  dissim_bray = vegdist(raw_counts_family, method="bray")
  dissim_bray_matr = as.matrix(dissim_bray)
  
  ### BC dissimilarity on filtered data 
  dissim_bray_filtered = vegdist(raw_counts_family_filtered, method="bray")
  dissim_bray_matr_filtered = as.matrix(dissim_bray_filtered) #... stored as matrix
  
  
  ## APC on full data
  apclust = apcluster(dissim_bray_matr, details = TRUE)
  #sink(file = "ap_output.txt")
  apclust
  #sink(file = NULL)
  
  
  apclust_filtered = apcluster(dissim_bray_matr_filtered, details = TRUE)
  #sink(file = "ap_output_filtered.txt")
  apclust_filtered
  #sink(file = NULL)

  
  cluster1 = as.data.frame(summary(as.factor(substr(rownames(as.data.frame(apclust_filtered@clusters[[1]])),1,1))))
  colnames(cluster1) = c("count")
  cluster1$farm = rownames(cluster1)
  
  farms <- as.factor(LETTERS[1:20])
  for(i in farms){
    if (i %notin% cluster1$farm){
      cluster1 = add_row(cluster1, count = 0, farm = i)
    }
  }
  
  cluster1$scores = prop.table(cluster1$count)
  rownames(cluster1) = cluster1$farm
  cluster1 = cluster1 %>% arrange(as.factor(farm))
  rownames(cluster1) = cluster1$farm    
  
  
  cluster2 = as.data.frame(summary(as.factor(substr(rownames(as.data.frame(apclust_filtered@clusters[[2]])),1,1))))
  colnames(cluster2) = c("count")
  cluster2$scores = prop.table(cluster2$count)
  
  farm_scores_bycluster <- data.frame("farm" = levels(as.factor(rownames(cluster2))), "Count_cluster1" = cluster1$count, "Count_cluster2" = cluster2$count)
  rownames(farm_scores_bycluster) = farm_scores_bycluster$farm
  
  farm_scores_bycluster_norm = decostand(farm_scores_bycluster[-1], method = "total", MARGIN = 2)
  farm_scores_bycluster_norm = apply(farm_scores_bycluster_norm, 2, function(x) round(x, 2))
  
  
  scores.tmm.nmds_family2$apclusters = 2
  first_ckust = names(apclust_filtered@clusters[[1]])
  for (j in scores.tmm.nmds_family2$site){
    if(j %in% first_ckust){scores.tmm.nmds_family2$apclusters[which(scores.tmm.nmds_family2$site==j)]=1}
  }
  scores.tmm.nmds_family2$apclusters = as.factor(scores.tmm.nmds_family2$apclusters)
  
  
  gg_cluster2_ap = ggplot() + 
    geom_text(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.025), size = 2) + 
    geom_point(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, colour = apclusters, fill = apclusters), size = 2, shape = 21) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_equal() +
    theme_bw()
  
  gg_cluster2_ap
  #ggsave("nmds_ap.eps", device = "eps", width=9, height=6)
  
  gg_cluster2_zoom_ap = ggplot() + 
    geom_text(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.013), size = 2.4) + 
    geom_point(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, colour = apclusters, fill = apclusters), size = 2, shape = 21) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_cartesian(ylim=c(- 0.25, 0.25), xlim = c(- 0.4, 0.4)) +
    theme_bw()
  gg_cluster2_zoom_ap
  #ggsave("nmds_zoom_ap.eps", device = "eps", width=9, height=6)
  
  gg_farm_ap = ggplot() + 
    geom_text(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.025), size = 2) + 
    geom_point(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, shape = apclusters), size = 2) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_equal() +
    theme_bw()
  gg_farm_ap
  #ggsave("nmds_farm_ap.eps", device = "eps", width=9, height=6)
  
  gg_farm_zoom_ap = ggplot() + 
    geom_text(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.013), size = 2.4) + 
    geom_point(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, shape = apclusters), size = 2) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_cartesian(ylim=c(- 0.25, 0.25), xlim = c(- 0.4, 0.4))  +
    theme_bw()
  gg_farm_zoom_ap
  #ggsave("nmds_farm_zoom_ap.eps", device = "eps", width=9, height=6)
  
  return(list(farm_scores_bycluster_norm, scores.tmm.nmds_family2))
}