tmm_norm <-function(genus_data){
  dge <- DGEList(genus_data)
  dge <- calcNormFactors(dge, method = "TMM")
  logCPM <- cpm(dge, log=TRUE) #cpm() produces (log2-)transformed normalised counts per million (log2 CPM)
  #plot(logCPM)
  counts.norm = as.data.frame(t(logCPM))
  counts.norm$farm = meta_data$Farm
  return(counts.norm)
}

nmds_hc_four_groups <- function(tmm_norm, animal_code1, animal_code2, animal_code3, animal_code4){
  set.seed(43)
  tmm_mds = metaMDS(comm = tmm_norm[,-417], k=2, distance = "bray", trymax=100) 
  bc_dist <- vegdist(tmm_norm[,-417], method = "bray")
  cluster_tmm <- hclust(bc_dist, method = 'ward.D2')
  
  ##choice of optimal number of clusters
  fviz_nbclust(tmm_norm[,-417], FUNcluster=hcut, diss=bc_dist, method = "silhouette")
  clust_tmm <- cutree(cluster_tmm, k = 2)
  clust_tmm = as.data.frame(clust_tmm)
  
  ## extract scores 
  scores.tmm.nmds = as.data.frame(scores(tmm_mds))
  scores.tmm.nmds$site = rownames(scores.tmm.nmds)
  scores.tmm.nmds$farm = as.factor(tmm_norm$farm)
  scores.tmm.nmds$cluster = clust_tmm$clust_tmm
  scores.tmm.nmds$cluster = as.factor(scores.tmm.nmds$cluster)
  scores.tmm.nmds$animal_group = NA
  scores.tmm.nmds$animal_group[animal_code1] = "Sows_14B"
  scores.tmm.nmds$animal_group[animal_code2] = "Sows_7B"
  scores.tmm.nmds$animal_group[animal_code3] = "Piglets_7B"
  scores.tmm.nmds$animal_group[animal_code4] = "Piglets_7W"
  scores.tmm.nmds$animal_group = as.factor(scores.tmm.nmds$animal_group)
  
  #head(scores.tmm.nmds)
  
  genus.tmm.nmds = as.data.frame(scores(tmm_mds, "species"))
  genus.tmm.nmds$genus = rownames(genus.tmm.nmds)
  #head(genus.tmm.nmds)
  
  g_clusters = ggplot() + 
    #  geom_text(data = genus.tmm.nmds, aes(x = NMDS1, y = NMDS2, label = genus), colour = "darkcyan") + 
    geom_point(data = scores.tmm.nmds, aes(x = NMDS1, y = NMDS2, colour = animal_group, fill = animal_group, shape = cluster), size = 2) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_equal() +
    theme_bw()
  
  g_groups = ggplot() + 
    # geom_text(data = genus.tmm.nmds, aes(x = NMDS1, y = NMDS2, label = genus), colour = "darkcyan") + 
    geom_point(data = scores.tmm.nmds, aes(x = NMDS1, y = NMDS2, colour = animal_group, fill = animal_group), size = 2, shape = 21) + # add the point markers
    #  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
    #  scale_colour_manual(values = gP) +
    coord_equal() +
    theme_bw()
  
  plot(g_groups)
  #ggsave("nmds_genus_allgroups.eps", device = "eps", width=9, height=5.5)
  
  fviz_nbclust(tmm_norm[,-417], FUNcluster=hcut, diss=bc_dist, method = "silhouette")
  #ggsave("hc_numberofclusters_genus_allgroups.eps", device = "eps", width=9, height=5.5)
  
  plot(g_clusters)
  #ggsave("nmds_hierarchicalclustering_genus_allgroups.eps", device = "eps", width=9, height=5.5)
}