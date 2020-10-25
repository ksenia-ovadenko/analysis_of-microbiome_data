

nmds_phylum <- function(data_microbiome, animal_ind){
  micro_data = data_microbiome[animal_ind,]
  micro_data_norm <- decostand(x = micro_data[,-24], method="total", MARGIN = 1)
  set.seed(45)
  micro.mds = metaMDS(comm = micro_data[, - 24], k=2, distance = "bray", trymax = 100) 
  
  
  data.scores.micro = as.data.frame(scores(micro.mds))
  data.scores.micro$site = rownames(data.scores.micro)
  data.scores.micro$farm = as.factor(micro_data$farm)
  
  species.scores.micro = as.data.frame(scores(micro.mds, "species"))
  species.scores.micro$species = rownames(species.scores.micro)
  
  data.scores.micro$samplesAndfarms = apply(data.scores.micro, 1, function(x) paste(x["site"], x["farm"], sep = ":"))
  
  int_micro = ggplot() + 
    geom_point_interactive(data = data.scores.micro, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) +
    geom_point(data = species.scores.micro, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
    geom_text_repel(data = species.scores.micro, aes(x = NMDS1, y = NMDS2, label = species)) +
    coord_equal() +
    theme_bw()
  ggiraph(ggobj = int_micro)
}


nmds_genus <- function(genus_norm_data, animal_ind, animal_code){
  genus_data = genus_norm_data[animal_ind,]
  genus_data_numeric = as.data.frame(apply(genus_data[,-417], 2, function(x) as.numeric(x)))
  genus_data_numeric$farm = genus_data$farm
  set.seed(45)
  micro.mds.genus = metaMDS(comm = genus_data_numeric[,-417], k=2) 
  
  data.scores.micro.genus = as.data.frame(scores(micro.mds.genus))
  data.scores.micro.genus$site = rownames(data.scores.micro.genus)
  data.scores.micro.genus$farm = as.factor(genus_data$farm)
  
  species.scores.micro.genus = as.data.frame(scores(micro.mds.genus, "species"))
  species.scores.micro.genus$species = rownames(species.scores.micro.genus)
  
  data.scores.micro.genus$samplesAndfarms = apply(data.scores.micro.genus, 1, function(x) paste(x["site"], x["farm"], sep = ":"))
  
  int_micro_genus = ggplot() + 
    #  geom_text(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
    #  geom_text_repel(data = species.scores.micro.S14Bgenus, aes(x = NMDS1, y = NMDS2, label = species)) +
    geom_point_interactive(data = data.scores.micro.genus, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
    geom_text(data = data.scores.micro.genus, aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.05, size = 2.5) +
    # geom_point(data = species.scores.micro.S14Bgenus, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
    coord_equal() +
    theme_bw()
  ggiraph(ggobj = int_micro_genus)
  #ggsave("nmds_genus_S14B.eps", device = "eps", width=9, height=5.5)
  
  if(animal_code == 'ind_S7B'){
    int_micro_S7Bgenus_zoom = ggplot() + 
      geom_point_interactive(data = subset(data.scores.micro.genus, data.scores.micro.genus$NMDS1 < 0.5 & data.scores.micro.genus$NMDS2 < 0.5 & data.scores.micro.genus$NMDS2 > -0.5), 
                             aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
      geom_text(data = subset(data.scores.micro.genus, data.scores.micro.genus$NMDS1 < 0.5 & data.scores.micro.genus$NMDS2 < 0.5 & data.scores.micro.genus$NMDS2 > -0.5), 
                aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.025, size = 2.5) +
      coord_cartesian(ylim=c(- 0.5, 0.5), xlim = c(- 0.5, 0.5)) +
      theme_bw()
    ggiraph(ggobj = int_micro_S7Bgenus_zoom)
  } else if(animal_code == 'ind_F7W'){
    int_micro_F7Wgenus_zoom = ggplot() + 
      geom_point_interactive(data = subset(data.scores.micro.genus, data.scores.micro.genus$NMDS1 > -1 & data.scores.micro.genus$NMDS1 < 1 & data.scores.micro.genus$NMDS2 > -1 & data.scores.micro.genus$NMDS2 < 1), 
                             aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
      geom_text(data = subset(data.scores.micro.genus, data.scores.micro.genus$NMDS1 > -1 & data.scores.micro.genus$NMDS1 < 1 & data.scores.micro.genus$NMDS2 > -1 & data.scores.micro.genus$NMDS2 < 1), 
                aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.05, size = 2.5) +
      coord_cartesian(ylim=c(- 1, 1), xlim = c(- 1, 1)) +
      theme_bw()
    ggiraph(ggobj = int_micro_F7Wgenus_zoom)
  }
}