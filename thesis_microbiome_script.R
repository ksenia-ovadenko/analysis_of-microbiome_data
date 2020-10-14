knitr::opts_chunk$set(echo = TRUE)
default_set = par(no.readonly = T)


library("ggplot2")
library("readxl")
#library("fastDummies")
library("vegan")
library("edgeR")
library("factoextra")


meta_data <- as.data.frame(read_xlsx("~/Desktop/Thesis/Optibiom/Optibiom_metadata_perSample_AL10022020.xlsx", 2))

meta_data[["Probensammlung"]] <- as.POSIXct(meta_data[["Probensammlung"]]*(60*60*24), origin = "1899-12-30")

meta_data[meta_data==-99] <- NA


## important meta variables
## medical initiation of birth
#meta_data_sows7B_med_einl <-  na.omit(meta_data[which(meta_data$Animal == "S" & meta_data$Time_sampling == "7"), c("Farm", "MedEinl")], na.action="omit")
#meta_data$MedEinl <- as.factor(meta_data$MedEinl)
t#able(meta_data_sows7B_med_einl$Farm, meta_data_sows7B_med_einl$MedEinl)

## asv classification data 
asv = as.data.frame(read_xlsx("~/Desktop/Thesis/Optibiom/Optibiom_ASV_table_AL_10022020.xlsx", 2))
asv.genus = as.factor(asv[,"Genus"])

#sink(file = "all_genera.txt")
#levels(asv.genus)
#sink(file = NULL)



# row.names(data_m) <- meta_data[,1]
# data_microbiome = data_microbiome[,-1]
# data_microbiome$farm = meta_data$Farm


## section 2.3 Modelling Microbiome Count Data
## download asv count data
asv_count <- as.data.frame(read_xlsx("~/Desktop/Thesis/asv/Optibiom_SchweinASV&taxonomy05_2020 short.xlsx"))

asv_count_genus <- asv_count[,-c(1:5)]
colnames(asv_count_genus) <-  sub("sample.", "", colnames(asv_count_genus))
groups <- asv_count_genus$Genus
sum_asv_count_genus <- as.data.frame(rowsum(x = asv_count_genus[,-1], group = groups))


## overdispersion plot 
raw_counts <- sum_asv_count_genus
df <- data.frame(mean = apply(raw_counts, 1, mean),
                 var = apply(raw_counts, 1, var))
#df <- df[df$mean <= 50, ]
p <- ggplot(data=df, aes(x = mean, y = var))
p <- p + geom_point(colour = "orange")
p <- p + theme_bw()
p <- p + geom_abline(aes(intercept=0, slope=1))
p <- p + ggtitle("Variance vs. Mean in Microbiome Counts") + ylab("variance")
p <- p + coord_cartesian(ylim=c(0, 1500), xlim=c(0, 50))
print(p)

#ggsave("overdispersion_raw_counts.eps", device = "eps", width=9, height=5.5)

ncol(sum_asv_count_genus)
asv_summary = sum_asv_count_genus
asv_summary$mean = round(apply(asv_summary, 1, mean), 2)
asv_summary = asv_summary[order(-asv_summary$mean),]
asv_summary = asv_summary[1:30,]
asv_summary$median = apply(asv_summary, 1, median)
asv_summary$variance = round(apply(asv_summary, 1, var), 2)
asv_summary$zero = apply(asv_summary[,-c(803,804,805)], 1, function(x){ round(length(which(x==0))/length(x)*100,2)})

summary_genera = as.data.frame(asv_summary[, c(803:806)])
rownames(summary_genera) = gsub("_", ".", rownames(summary_genera))
colnames(summary_genera) = c("Mean", "Median", "Variance", "Zero.pct")
summary_genera
#write.csv(summary_genera, file="summary_gen.csv", quote=FALSE)

library("vegan")
library("FactoMineR")
library("ggplot2")
library("ggrepel")
library("ggiraph")

## download microbiome data (phylum level)
data_microbiome <- as.data.frame(read_xlsx("~/Desktop/Thesis/Optibiom/Optibiom_tax_profile_AL_10022020.xlsx"), header = TRUE, row.names = 1) ##phylum level
row.names(data_microbiome) <- meta_data[,1]
data_microbiome = data_microbiome[,-1]
data_microbiome$farm = meta_data$Farm


micro_data.S14B = data_microbiome[which(meta_data$Animal == "S" & meta_data$Time_sampling == "-14"),]
micro_data.S14B_norm <- decostand(x = micro_data.S14B[,-24], method="total", MARGIN = 1)
sum(micro_data.S14B_norm[1,])
set.seed(45)
micro.mds.S14B = metaMDS(comm = micro_data.S14B[, - 24], k=2, distance = "bray", trymax = 100) 


data.scores.micro.S14B = as.data.frame(scores(micro.mds.S14B))
data.scores.micro.S14B$site = rownames(data.scores.micro.S14B)
data.scores.micro.S14B$farm = as.factor(micro_data.S14B$farm)
#head(data.scores.micro.S14B)

species.scores.micro.S14B = as.data.frame(scores(micro.mds.S14B, "species"))
species.scores.micro.S14B$species = rownames(species.scores.micro.S14B)
#head(species.scores.micro.S14B)

data.scores.micro.S14B$samplesAndfarms = apply(data.scores.micro.S14B, 1, function(x) paste(x["site"], x["farm"], sep = ":"))

int_micro_S14B = ggplot() + 
#  geom_text(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +
  geom_point_interactive(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + # add the point markers
  geom_point(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
  geom_text_repel(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species)) +
  coord_equal() +
  theme_bw()
ggiraph(ggobj = int_micro_S14B)

#ggsave("int_micro_S14B.eps", device = "eps", width=9, height=5.5)


# int_micro_S14B_zoom = ggplot() + 
# #  geom_text(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
#   geom_point_interactive(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + # add the point markers
#   geom_point(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
#   geom_text_repel(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species)) +
#   coord_cartesian(ylim=c(- 1, 1), xlim = c(- 1, 1)) +
#   theme_bw()
# 
# ggiraph(ggobj = int_micro_S14B_zoom)
# ggsave("int_micro_S14B_zoom.eps", device = "eps", width=9, height=5.5)



micro_data.S7B = data_microbiome[which(meta_data$Animal == "S" & meta_data$Time_sampling == "7"),]
micro_data.S7B_norm <- decostand(x = micro_data.S7B[, -24], method="total", MARGIN = 1)
set.seed(45)
micro.mds.S7B = metaMDS(comm = micro_data.S7B[, - 24], k=2, distance = "bray", trymax = 100) 

## extract scores 
data.scores.micro.S7B = as.data.frame(scores(micro.mds.S7B))
data.scores.micro.S7B$site = rownames(data.scores.micro.S7B)
data.scores.micro.S7B$farm = as.factor(micro_data.S7B$farm)

species.scores.micro.S7B = as.data.frame(scores(micro.mds.S7B, "species"))
species.scores.micro.S7B$species = rownames(species.scores.micro.S7B)

data.scores.micro.S7B$samplesAndfarms = apply(data.scores.micro.S7B, 1, function(x) paste(x["site"], x["farm"], sep = ":"))

plot_int_S7B = ggplot() + 
  geom_point_interactive(data = data.scores.micro.S7B, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text_repel(data = species.scores.micro.S7B, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point(data = species.scores.micro.S7B, aes(x = NMDS1, y = NMDS2), size = 2) +
  coord_equal() +
  theme_bw()

ggiraph(ggobj = plot_int_S7B)
#ggsave("int_micro_S7B.eps", device = "eps", width=9, height=5.5)

# 
# #zoomed
# plot_int_S7Bz = ggplot() + 
#   geom_point_interactive(data = data.scores.micro.S7B, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
#   geom_text_repel(data = species.scores.micro.S7B, aes(x = NMDS1, y = NMDS2, label = species)) +
#   geom_point(data = species.scores.micro.S7B, aes(x = NMDS1, y = NMDS2), size = 2) +
#   coord_cartesian(ylim=c(- 1, 0.75), xlim = c(- 1, 1)) +
#   theme_bw()
# 
# ggiraph(ggobj = plot_int_S7Bz)


micro_data.F7B = data_microbiome[which(meta_data$Animal == "F" & meta_data$Time_sampling == "7"& meta_data$Time_points == "B"),]
micro_data.F7B_norm <- decostand(x = micro_data.F7B[, -24], method="total", MARGIN = 1)
set.seed(45)
micro.mds.F7B = metaMDS(comm = micro_data.F7B[, - 24], k=2, distance = "bray", trymax = 100) 

## extract scores 
data.scores.micro.F7B = as.data.frame(scores(micro.mds.F7B))
data.scores.micro.F7B$site = rownames(data.scores.micro.F7B)
data.scores.micro.F7B$farm = as.factor(micro_data.F7B$farm)

species.scores.micro.F7B = as.data.frame(scores(micro.mds.F7B, "species"))
species.scores.micro.F7B$species = rownames(species.scores.micro.F7B)

data.scores.micro.F7B$samplesAndfarms = apply(data.scores.micro.F7B, 1, function(x) paste(x["site"], x["farm"], sep = ":"))


##plot
plot_int_F7B = ggplot() + 
  geom_point_interactive(data = data.scores.micro.F7B, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text_repel(data = species.scores.micro.F7B, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point(data = species.scores.micro.F7B, aes(x = NMDS1, y = NMDS2), size = 2) +
  coord_equal() +
  theme_bw()

ggiraph(ggobj = plot_int_F7B)
#ggsave("int_micro_F7B.eps", device = "eps", width=9, height=5.5)

# #zoomed
# plot_int_F7Bz = ggplot() +
#   geom_point_interactive(data = subset(data.scores.micro.F7B, data.scores.micro.F7B$NMDS1 < 1 & data.scores.micro.F7B$NMDS2 < 1), aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) +
#   geom_text_repel(data = subset(species.scores.micro.F7B, species.scores.micro.F7B$NMDS1 < 1 & species.scores.micro.F7B$NMDS2 < 1), aes(x = NMDS1, y = NMDS2, label = species)) +
#   geom_point(data = subset(species.scores.micro.F7B, species.scores.micro.F7B$NMDS1 < 1& species.scores.micro.F7B$NMDS2 < 1), aes(x = NMDS1, y = NMDS2), size = 2) +
#   coord_cartesian(ylim=c(- 0.75, 0.75), xlim = c(- 0.75, 1)) +
#   theme_bw()
# 
# ggiraph(ggobj = plot_int_F7Bz)


micro_data.F7W = data_microbiome[which(meta_data$Animal == "F" & meta_data$Time_sampling == "7"& meta_data$Time_points == "W"),]
micro_data.F7W_norm <- decostand(x = micro_data.F7W[, -24], method="total", MARGIN = 1)
set.seed(45)
micro.mds.F7W = metaMDS(comm = micro_data.F7W[, - 24], k=2, trymax = 100) 

## extract scores 
data.scores.micro.F7W = as.data.frame(scores(micro.mds.F7W))
data.scores.micro.F7W$site = rownames(data.scores.micro.F7W)
data.scores.micro.F7W$farm = as.factor(micro_data.F7W$farm)

species.scores.micro.F7W = as.data.frame(scores(micro.mds.F7W, "species"))
species.scores.micro.F7W$species = rownames(species.scores.micro.F7W)

data.scores.micro.F7W$samplesAndfarms = apply(data.scores.micro.F7W, 1, function(x) paste(x["site"], x["farm"], sep = ":"))


##plot
plot_int_F7W = ggplot() + 
  geom_point_interactive(data = data.scores.micro.F7W, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text_repel(data = species.scores.micro.F7W, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point(data = species.scores.micro.F7W, aes(x = NMDS1, y = NMDS2), size = 2) +
  coord_equal() +
  theme_bw()

ggiraph(ggobj = plot_int_F7W)
#ggsave("int_micro_F7W.eps", device = "eps", width=9, height=5.5)

#zoomed
# plot_int_F7Bz = ggplot() +
#   geom_point_interactive(data = subset(data.scores.micro.F7W, data.scores.micro.F7W$NMDS1 < 1 & data.scores.micro.F7W$NMDS2 < 1), aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) +
#   geom_text_repel(data = subset(species.scores.micro.F7W, species.scores.micro.F7W$NMDS1 < 1 & species.scores.micro.F7W$NMDS2 < 1), aes(x = NMDS1, y = NMDS2, label = species)) +
#   geom_point(data = subset(species.scores.micro.F7W, species.scores.micro.F7W$NMDS1 < 1& species.scores.micro.F7W$NMDS2 < 1), aes(x = NMDS1, y = NMDS2), size = 2) +
#   coord_cartesian(ylim=c(- 0.75, 0.75), xlim = c(- 0.75, 1)) +
#   theme_bw()
# 
# ggiraph(ggobj = plot_int_F7Bz)


## asv data (genus level)
asv_count <- as.data.frame(read_xlsx("~/Desktop/Thesis/asv/Optibiom_SchweinASV&taxonomy05_2020 short.xlsx"))
asv_count_genus <- asv_count[,-c(1:5)]
#colnames(asv_count_genus) <-  sub("sample.", "", colnames(asv_count_genus))
#colnames(asv_count_genus)
#asv_count_genus_norm <- decostand(x = asv_count_genus[,-1], method="total", MARGIN = 2)
#asv_count_genus_norm$Genus <- as.factor(asv_count_genus$Genus)
#groups <- asv_count_genus_norm$Genus
#sum_asv_count_genus_norm <- as.data.frame(rowsum(x = asv_count_genus_norm[,-803], group = groups))

asv_count_genus$Genus <- as.factor(asv_count_genus$Genus)
groups <-asv_count_genus$Genus
sum_asv_count_genus <- as.data.frame(rowsum(x = asv_count_genus[,-1], group = groups))


sum_asv_count_genus_norm <- decostand(x = sum_asv_count_genus, method="total", MARGIN = 2)
dim(sum_asv_count_genus_norm)

sum_asv_count_genus_norm$Genus <- levels(asv_count_genus$Genus)
rownames(sum_asv_count_genus_norm) = sum_asv_count_genus_norm$Genus
sum_asv_count_genus_norm = sum_asv_count_genus_norm[,-803]

sum_asv_count_genus_norm <- as.data.frame(t(sum_asv_count_genus_norm))
sum_asv_count_genus_norm$farm = meta_data$Farm
#write.csv(sum_asv_count_genus_norm, file="genera_norm.csv", quote=FALSE)

## nmds based on genus level
genus_data.S14B = sum_asv_count_genus_norm[which(meta_data$Animal == "S" & meta_data$Time_sampling == "-14"),]
genus_data.S14B_numeric = as.data.frame(apply(genus_data.S14B[,-417], 2, function(x) as.numeric(x)))
genus_data.S14B_numeric$farm = genus_data.S14B$farm
set.seed(45)
micro.mds.S14Bgenus = metaMDS(comm = genus_data.S14B_numeric[,-417], k=2) 

#apply(genus_data.S14B_numeric, 2, function(x) is.numeric(x))

data.scores.micro.S14Bgenus = as.data.frame(scores(micro.mds.S14Bgenus))
data.scores.micro.S14Bgenus$site = rownames(data.scores.micro.S14Bgenus)
data.scores.micro.S14Bgenus$farm = as.factor(genus_data.S14B$farm)
#head(data.scores.micro.S14B)

species.scores.micro.S14Bgenus = as.data.frame(scores(micro.mds.S14Bgenus, "species"))
species.scores.micro.S14Bgenus$species = rownames(species.scores.micro.S14Bgenus)
#head(species.scores.micro.S14B)

data.scores.micro.S14Bgenus$samplesAndfarms = apply(data.scores.micro.S14Bgenus, 1, function(x) paste(x["site"], x["farm"], sep = ":"))

int_micro_S14Bgenus = ggplot() + 
#  geom_text(data = species.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
#  geom_text_repel(data = species.scores.micro.S14Bgenus, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point_interactive(data = data.scores.micro.S14Bgenus, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text(data = data.scores.micro.S14Bgenus, aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.05, size = 2.5) +
 # geom_point(data = species.scores.micro.S14Bgenus, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
  coord_equal() +
  theme_bw()
ggiraph(ggobj = int_micro_S14Bgenus)
#ggsave("nmds_genus_S14B.eps", device = "eps", width=9, height=5.5)





genus_data.S7B = sum_asv_count_genus_norm[which(meta_data$Animal == "S" & meta_data$Time_sampling == "7"),]
genus_data.S7B_numeric = as.data.frame(apply(genus_data.S7B[,-417], 2, function(x) as.numeric(x)))
genus_data.S7B_numeric$farm = genus_data.S7B$farm
set.seed(45)
micro.mds.S7Bgenus = metaMDS(comm = genus_data.S7B_numeric[,-417], k=2) 

#apply(genus_data.S7B_numeric, 2, function(x) is.numeric(x))

data.scores.micro.S7Bgenus = as.data.frame(scores(micro.mds.S7Bgenus))
data.scores.micro.S7Bgenus$site = rownames(data.scores.micro.S7Bgenus)
data.scores.micro.S7Bgenus$farm = as.factor(genus_data.S7B$farm)
#head(data.scores.micro.S7B)

species.scores.micro.S7Bgenus = as.data.frame(scores(micro.mds.S7Bgenus, "species"))
species.scores.micro.S7Bgenus$species = rownames(species.scores.micro.S7Bgenus)
#head(species.scores.micro.S7B)

data.scores.micro.S7Bgenus$samplesAndfarms = apply(data.scores.micro.S7Bgenus, 1, function(x) paste(x["site"], x["farm"], sep = ":"))

int_micro_S7Bgenus = ggplot() + 
#  geom_text(data = species.scores.micro.S7B, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
#  geom_text_repel(data = species.scores.micro.S7Bgenus, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point_interactive(data = data.scores.micro.S7Bgenus, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text(data = data.scores.micro.S7Bgenus, aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.05, size = 2.5) +
 # geom_point(data = species.scores.micro.S7Bgenus, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
  coord_equal() +
  theme_bw()
ggiraph(ggobj = int_micro_S7Bgenus)
#ggsave("nmds_genus_S7B.eps", device = "eps", width=9, height=5.5)

int_micro_S7Bgenus_zoom = ggplot() + 
#  geom_text(data = species.scores.micro.S7B, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
#  geom_text_repel(data = species.scores.micro.S7Bgenus, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point_interactive(data = subset(data.scores.micro.S7Bgenus, data.scores.micro.S7Bgenus$NMDS1 < 0.5 & data.scores.micro.S7Bgenus$NMDS2 < 0.5 & data.scores.micro.S7Bgenus$NMDS2 > -0.5), aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text(data = subset(data.scores.micro.S7Bgenus, data.scores.micro.S7Bgenus$NMDS1 < 0.5 & data.scores.micro.S7Bgenus$NMDS2 < 0.5 & data.scores.micro.S7Bgenus$NMDS2 > -0.5), aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.025, size = 2.5) +
 # geom_point(data = species.scores.micro.S7Bgenus, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
  coord_cartesian(ylim=c(- 0.5, 0.5), xlim = c(- 0.5, 0.5)) +
  theme_bw()
ggiraph(ggobj = int_micro_S7Bgenus_zoom)
#ggsave("nmds_genus_S7B_zoom.eps", device = "eps", width=9, height=5.5)






genus_data.F7B = sum_asv_count_genus_norm[which(meta_data$Animal == "F" & meta_data$Time_sampling == "7" & meta_data$Time_points == "B"),]
genus_data.F7B_numeric = as.data.frame(apply(genus_data.F7B[,-417], 2, function(x) as.numeric(x)))
genus_data.F7B_numeric$farm = genus_data.F7B$farm
set.seed(45)
micro.mds.F7Bgenus = metaMDS(comm = genus_data.F7B_numeric[,-417], k=2) 

#apply(genus_data.F7B_numeric, 2, function(x) is.numeric(x))

data.scores.micro.F7Bgenus = as.data.frame(scores(micro.mds.F7Bgenus))
data.scores.micro.F7Bgenus$site = rownames(data.scores.micro.F7Bgenus)
data.scores.micro.F7Bgenus$farm = as.factor(genus_data.F7B$farm)
#head(data.scores.micro.F7B)

species.scores.micro.F7Bgenus = as.data.frame(scores(micro.mds.F7Bgenus, "species"))
species.scores.micro.F7Bgenus$species = rownames(species.scores.micro.F7Bgenus)
#head(species.scores.micro.F7B)

data.scores.micro.F7Bgenus$samplesAndfarms = apply(data.scores.micro.F7Bgenus, 1, function(x) paste(x["site"], x["farm"], sep = ":"))

int_micro_F7Bgenus = ggplot() + 
#  geom_text(data = species.scores.micro.F7B, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
#  geom_text_repel(data = species.scores.micro.F7Bgenus, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point_interactive(data = data.scores.micro.F7Bgenus, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text(data = data.scores.micro.F7Bgenus, aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.05, size = 2.5) +
 # geom_point(data = species.scores.micro.F7Bgenus, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
  coord_equal() +
  theme_bw()
ggiraph(ggobj = int_micro_F7Bgenus)
#ggsave("nmds_genus_F7B.eps", device = "eps", width=9, height=5.5)




genus_data.F7W = sum_asv_count_genus_norm[which(meta_data$Animal == "F" & meta_data$Time_sampling == "7" & meta_data$Time_points == "W"),]
genus_data.F7W_numeric = as.data.frame(apply(genus_data.F7W[,-417], 2, function(x) as.numeric(x)))
genus_data.F7W_numeric$farm = genus_data.F7W$farm
set.seed(45)
micro.mds.F7Wgenus = metaMDS(comm = genus_data.F7W_numeric[,-417], k=2) 

#apply(genus_data.F7W_numeric, 2, function(x) is.numeric(x))

data.scores.micro.F7Wgenus = as.data.frame(scores(micro.mds.F7Wgenus))
data.scores.micro.F7Wgenus$site = rownames(data.scores.micro.F7Wgenus)
data.scores.micro.F7Wgenus$farm = as.factor(genus_data.F7W$farm)
#head(data.scores.micro.F7W)

species.scores.micro.F7Wgenus = as.data.frame(scores(micro.mds.F7Wgenus, "species"))
species.scores.micro.F7Wgenus$species = rownames(species.scores.micro.F7Wgenus)
#head(species.scores.micro.F7W)

data.scores.micro.F7Wgenus$samplesAndfarms = apply(data.scores.micro.F7Wgenus, 1, function(x) paste(x["site"], x["farm"], sep = ":"))

int_micro_F7Wgenus = ggplot() + 
#  geom_text(data = species.scores.micro.F7W, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
#  geom_text_repel(data = species.scores.micro.F7Wgenus, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point_interactive(data = data.scores.micro.F7Wgenus, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text(data = data.scores.micro.F7Wgenus, aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.05, size = 2.5) +
 # geom_point(data = species.scores.micro.F7Wgenus, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
  coord_equal() +
  theme_bw()
ggiraph(ggobj = int_micro_F7Wgenus)
#ggsave("nmds_genus_F7W.eps", device = "eps", width=9, height=5.5)


int_micro_F7Wgenus_zoom = ggplot() + 
#  geom_text(data = species.scores.micro.F7W, aes(x = NMDS1, y = NMDS2, label = species), colour = "darkcyan") +  
#  geom_text_repel(data = species.scores.micro.F7Wgenus, aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point_interactive(data = subset(data.scores.micro.F7Wgenus, data.scores.micro.F7Wgenus$NMDS1 > -1 & data.scores.micro.F7Wgenus$NMDS1 < 1 & data.scores.micro.F7Wgenus$NMDS2 > -1 & data.scores.micro.F7Wgenus$NMDS2 < 1), aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, tooltip = samplesAndfarms), size = 3, shape = 21) + 
  geom_text(data = subset(data.scores.micro.F7Wgenus, data.scores.micro.F7Wgenus$NMDS1 > -1 & data.scores.micro.F7Wgenus$NMDS1 < 1 & data.scores.micro.F7Wgenus$NMDS2 > -1 & data.scores.micro.F7Wgenus$NMDS2 < 1), aes(x = NMDS1, y = NMDS2, label = farm), nudge_y = 0.05, size = 2.5) +
 # geom_point(data = species.scores.micro.F7Wgenus, aes(x = NMDS1, y = NMDS2, label = species), size = 2) +
  coord_cartesian(ylim=c(- 1, 1), xlim = c(- 1, 1)) +
  theme_bw()
ggiraph(ggobj = int_micro_F7Wgenus_zoom)
#ggsave("nmds_genus_F7W_zoom.eps", device = "eps", width=9, height=5.5)

library("edgeR")
library("ggplot2")
library("ggrepel")

dge <- DGEList(sum_asv_count_genus)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log=TRUE) #cpm() produces (log2-)transformed normalised counts per million (log2 CPM)
plot(logCPM)

counts.norm = as.data.frame(t(logCPM))
counts.norm$farm = meta_data$Farm


set.seed(43)
tmm_mds = metaMDS(comm = counts.norm[,-417], k=2, distance = "bray", trymax=100) 
#tmm_mds_asv = metaMDS(comm = sum_asv_count_genus_norm[,-417], k=2, distance = "bray") 
#ordiplot(tmm_mds, display = "sites", type = "point")
#ordiplot(tmm_mds_asv, display = "sites", type = "point")
bc_dist <- vegdist(counts.norm[,-417], method = "bray")
cluster_tmm <- hclust(bc_dist, method = 'ward.D2')
#plot(cluster_tmm) 

fviz_nbclust(counts.norm[,-417], FUNcluster=hcut, diss=bc_dist, method = "silhouette")
#fviz_nbclust(counts.norm[,-417], FUNcluster=hcut, diss=bc_dist, method = "wss") ## elbow approx at k=3
#fviz_nbclust(counts.norm[,-417], FUNcluster=hcut, diss=bc_dist, method = "gap_stat")#

colnames(counts.norm)
clust_tmm <- cutree(cluster_tmm, k = 2)
clust_tmm = as.data.frame(clust_tmm)


ind_S14B = as.logical(meta_data$Animal == "S" & meta_data$Time_sampling == "-14")
ind_S7B = as.logical(meta_data$Animal == "S" & meta_data$Time_sampling == "7")
ind_F7B = as.logical(meta_data$Animal == "F" & meta_data$Time_sampling == "7"& meta_data$Time_points == "B")
ind_F7W = as.logical(meta_data$Animal == "F" & meta_data$Time_sampling == "7"& meta_data$Time_points == "W")




## extract scores 
scores.tmm.nmds = as.data.frame(scores(tmm_mds))
scores.tmm.nmds$site = rownames(scores.tmm.nmds)
scores.tmm.nmds$farm = as.factor(counts.norm$farm)
scores.tmm.nmds$cluster = clust_tmm$clust_tmm
scores.tmm.nmds$cluster = as.factor(scores.tmm.nmds$cluster)
scores.tmm.nmds$animal_group = NA
scores.tmm.nmds$animal_group[ind_S14B] = "Sows_14B"
scores.tmm.nmds$animal_group[ind_S7B] = "Sows_7B"
scores.tmm.nmds$animal_group[ind_F7B] = "Piglets_7B"
scores.tmm.nmds$animal_group[ind_F7W] = "Piglets_7W"
scores.tmm.nmds$animal_group = as.factor(scores.tmm.nmds$animal_group)

#head(scores.tmm.nmds)

genus.tmm.nmds = as.data.frame(scores(tmm_mds, "species"))
genus.tmm.nmds$genus = rownames(genus.tmm.nmds)
#head(genus.tmm.nmds)



##plot
#library("randomcoloR")
#cols = distinctColorPalette(20)


#library(RColorBrewer)
#getPalette = colorRampPalette(brewer.pal(12, "Set3"))
#gP = getPalette(nlevels(samples.scores$farm))
#names(gP) = levels(samples.scores$farm)


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

fviz_nbclust(counts.norm[,-417], FUNcluster=hcut, diss=bc_dist, method = "silhouette")
#ggsave("hc_numberofclusters_genus_allgroups.eps", device = "eps", width=9, height=5.5)

plot(g_clusters)
#ggsave("nmds_hierarchicalclustering_genus_allgroups.eps", device = "eps", width=9, height=5.5)

library("edgeR")
library("vegan")
library("factoextra")
library("utils")

df = read.delim(file = "~/Desktop/Thesis/stacked data/stacked_data.txt", sep = "\t")
df_na = df[complete.cases(df), ]

raw_counts_family = df_na
#sapply(df_na, function(x) sum(is.na(x)))

rownames(raw_counts_family) = with(raw_counts_family, paste0(farm, individual))
raw_counts_family = raw_counts_family[,-c(1,2)]
raw_counts_family_matrix = as.matrix(t(raw_counts_family))

raw_counts_wo_normalization = as.data.frame(raw_counts_family_matrix)
## tmm normalization
dge_family = DGEList(raw_counts_family_matrix)
dge_family <- calcNormFactors(dge_family, method = "TMM")
#dge2$samples

pseudo_TMM_family <- log2(cpm(dge_family) + 1)
empty_rows = apply(pseudo_TMM_family, 1, function(x) sum(x)==0)
pseudo_TMM_family <- pseudo_TMM_family[!empty_rows,]

raw_counts_wo_normalization = raw_counts_wo_normalization[!empty_rows,]

pseudo_TMM_family_nmds = as.data.frame(t(pseudo_TMM_family))
pseudo_TMM_family_nmds$farm = as.factor(substr(row.names(pseudo_TMM_family_nmds), 1, 1))
#dim(pseudo_TMM_family_nmds)
#sapply(raw_counts_family, function(x) sum(x)==0)
set.seed(23)
tmm_mds_family = metaMDS(comm = pseudo_TMM_family_nmds[,-1169], k=2, distance = "bray", trymax = 100) 


## HClust
bc_dist_family <- vegdist(pseudo_TMM_family_nmds[,-1169], method = "bray")
cluster_tmm_family <- hclust(bc_dist_family, method = 'ward.D2')
#plot(cluster_tmm_family) 

fviz_nbclust(pseudo_TMM_family_nmds[,-1169], FUNcluster=hcut, diss=bc_dist_family, method = "silhouette") ## k=3
#ggsave("optimal_n_clusters_genus_family.eps", device = "eps", width=9, height=5.5)
#fviz_nbclust(pseudo_TMM_family_nmds[,-1169], FUNcluster=hcut, diss=bc_dist_family, method = "wss") 
#fviz_nbclust(pseudo_TMM_family_nmds[,-1169], FUNcluster=hcut, diss=bc_dist_family, method = "gap_stat")

clust_tmm_family <- cutree(cluster_tmm_family, k = 3)
clust_tmm_family = as.data.frame(clust_tmm_family)



## nmds scores for plots
scores.tmm.nmds_family = as.data.frame(scores(tmm_mds_family))
scores.tmm.nmds_family$site = rownames(scores.tmm.nmds_family)
scores.tmm.nmds_family$farm = as.factor(pseudo_TMM_family_nmds$farm)
scores.tmm.nmds_family$cluster = clust_tmm_family$clust_tmm_family
scores.tmm.nmds_family$cluster = as.factor(scores.tmm.nmds_family$cluster)


genus.tmm.nmds_family = as.data.frame(scores(tmm_mds_family, "species"))
genus.tmm.nmds_family$genus = rownames(genus.tmm.nmds_family)
head(genus.tmm.nmds_family)


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
#ggsave("nmds_hc_genus_family_farmcolour_zoom.eps", device = "eps", width=9, height=6)

cluster1_family = scores.tmm.nmds_family[which(scores.tmm.nmds_family$cluster==1),]
cluster2_family = scores.tmm.nmds_family[which(scores.tmm.nmds_family$cluster==2),]
cluster3_family = scores.tmm.nmds_family[which(scores.tmm.nmds_family$cluster==3),]

library("BBmisc")
scores1_family = as.data.frame(summary(cluster1_family$farm))
scores1_family$scores1 = prop.table(scores1_family$`summary(cluster1_family$farm)`)
#scores1_family$`summary(cluster1_family$farm)`
colnames(scores1_family) = c("count", "scores")

scores2_family = as.data.frame(summary(cluster2_family$farm))
scores2_family$scores2 = prop.table(scores2_family$`summary(cluster2_family$farm)`)
#scores2_family$`summary(cluster2_family$farm)`
colnames(scores2_family) = c("count", "scores")

scores3_family = as.data.frame(summary(cluster3_family$farm))
scores3_family$scores3 = prop.table(scores3_family$`summary(cluster3_family$farm)`)
#scores3_family$`summary(cluster3_family$farm)`
colnames(scores3_family) = c("count", "scores")

farm_scores_family <- data.frame("farm" = levels(cluster1_family$farm), "Scores_cluster1" = round(scores1_family$scores, 2), "Scores_cluster2" = round(scores2_family$scores, 2), "Scores_cluster3" = round(scores3_family$scores, 2))

farm_scores_family 

clust_tmm_family2 <- cutree(cluster_tmm_family, k = 2)
clust_tmm_family2 = as.data.frame(clust_tmm_family2)

plot(cluster_tmm_family)

## nmds scores for plots
scores.tmm.nmds_family2 = as.data.frame(scores(tmm_mds_family))
scores.tmm.nmds_family2$site = rownames(scores.tmm.nmds_family)
scores.tmm.nmds_family2$farm = as.factor(pseudo_TMM_family_nmds$farm)
scores.tmm.nmds_family2$cluster = clust_tmm_family2$clust_tmm_family2
scores.tmm.nmds_family2$cluster = as.factor(scores.tmm.nmds_family2$cluster)


genus.tmm.nmds_family2 = as.data.frame(scores(tmm_mds_family, "species"))
genus.tmm.nmds_family2$genus = rownames(genus.tmm.nmds_family2)



gg_cluster2 = ggplot() + 
  geom_text(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.025), size = 2) + 
  geom_point(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, colour = cluster, fill = cluster), size = 2, shape = 21) + # add the point markers
#  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
#  scale_colour_manual(values = gP) +
  coord_equal() +
  theme_bw()

gg_cluster2
#ggsave("nmds_hc_genus_family2.eps", device = "eps", width=9, height=6)



gg_cluster2_zoom = ggplot() + 
  geom_text(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.013), size = 2.4) + 
  geom_point(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, colour = cluster, fill = cluster), size = 2, shape = 21) + # add the point markers
#  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
#  scale_colour_manual(values = gP) +
  coord_cartesian(ylim=c(- 0.25, 0.25), xlim = c(- 0.4, 0.4)) +
  theme_bw()
gg_cluster2_zoom
#ggsave("nmds_hc_genus_family_zoom2.eps", device = "eps", width=9, height=6)

gg_farm2 = ggplot() + 
 geom_text(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.025), size = 2) + 
  geom_point(data = scores.tmm.nmds_family2, aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, shape = cluster), size = 2) + # add the point markers
#  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
#  scale_colour_manual(values = gP) +
  coord_equal() +
  theme_bw()
gg_farm2
#ggsave("nmds_hc_genus_family_farmcolour2.eps", device = "eps", width=9, height=6)

gg_farm2_zoom = ggplot() + 
 geom_text(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, label = site), position = position_nudge(y = -0.013), size = 2.4) + 
  geom_point(data = subset(scores.tmm.nmds_family2, scores.tmm.nmds_family2$NMDS1 < 0.4 & scores.tmm.nmds_family2$NMDS1 > -0.4 & scores.tmm.nmds_family2$NMDS2 < 0.25 & scores.tmm.nmds_family2$NMDS1 > -0.25), aes(x = NMDS1, y = NMDS2, colour = farm, fill = farm, shape = cluster), size = 2) + # add the point markers
#  geom_text(data = data.scores.micro.S14B, aes(x = NMDS1, y = NMDS2, label = site), size = 6, vjust = 0) +  # add the site labels
#  scale_colour_manual(values = gP) +
  coord_cartesian(ylim=c(- 0.25, 0.25), xlim = c(- 0.4, 0.4))  +
  theme_bw()
gg_farm2_zoom
#ggsave("nmds_hc_genus_family_farmcolour_zoom2.eps", device = "eps", width=9, height=6)

library("dendextend")
#cluster_tmm_family <- hclust(bc_dist_family, method = 'ward.D2')
#clust_tmm_family2 <- cutree(cluster_tmm_family, k = 2)
#clust_tmm_family2 = as.data.frame(clust_tmm_family2)



dend = as.dendrogram(cluster_tmm_family)
# color the branches based on the clusters
dend3 = dend
dend3 = color_branches(dend3, k=3, groupLabels = c("Group 1", "Group 3", "Group 2"))
dend3 = set(dend3, "labels_col", "white")
#setEPS()
postscript("tree3.eps", height = 5, width = 9)
#plot(dend3, main = "Tree of Hierarchical Clustering for 3 Clusters")
#dev.off()



dend2 = dend
dend2 = color_branches(dend2, k=2, groupLabels = c("Group 1", "Group 2"))
dend2 = set(dend2, "labels_col", "white")
#setEPS()
#postscript("tree2.eps", height = 5, width = 9)
plot(dend2, main = "Tree of Hierarchical Clustering for 2 Clusters")
#dev.off()

cluster1_family2 = scores.tmm.nmds_family[which(scores.tmm.nmds_family2$cluster==1),]
cluster2_family2 = scores.tmm.nmds_family[which(scores.tmm.nmds_family2$cluster==2),]

library("BBmisc")
scores1_family2 = as.data.frame(summary(cluster1_family2$farm))
scores1_family2$scores1 = prop.table(scores1_family2$`summary(cluster1_family2$farm)`)
#scores1_family$`summary(cluster1_family$farm)`
colnames(scores1_family2) = c("count", "scores")

scores2_family2 = as.data.frame(summary(cluster2_family2$farm))
scores2_family2$scores2 = prop.table(scores2_family2$`summary(cluster2_family2$farm)`)
#scores2_family$`summary(cluster2_family$farm)`
colnames(scores2_family2) = c("count", "scores")


farm_scores_family2 <- data.frame("farm" = levels(cluster1_family2$farm), "Scores_cluster1" = round(scores1_family2$scores, 2), "Scores_cluster2" = round(scores2_family2$scores, 2))


farm_scores_family2 

library("mixOmics")
library("factoextra")
library("apcluster") 
library("propr")
library("vegan")



raw_counts_family = df_na
rownames(raw_counts_family) = with(raw_counts_family, paste0(farm, individual))
raw_counts_family = raw_counts_family[,-c(1,2)]
raw_counts_family_pseudo = raw_counts_family+1


## log-ratio transformation
library("propr")
clr_trans_data = propr(as.matrix(raw_counts_family_pseudo))
clr_data = clr_trans_data@logratio

#lib.size <- apply(df_na[,-c(1,2)], 1, sum)
# barplot(lib.size)



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


library("mefa4")
library("tibble")
library("DataCombine")
library("dplyr")

cluster1 = as.data.frame(summary(as.factor(substr(rownames(as.data.frame(apclust_filtered@clusters[[1]])),1,1))))
colnames(cluster1) = c("count")
cluster1$farm = rownames(cluster1)

for(i in levels(as.factor(df_na$farm))){
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
apply(farm_scores_bycluster_norm, 2, function(x) round(x, 2))



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

library("DESeq2")
#head(otu_tab)
#dim(raw_counts_family)


otu_tab = raw_counts_family

countData <- as(otu_tab, "matrix")
#head(countData)

countData<-(t(countData)) #DESeq2 need taxa(genes=rows) by samples(=columns) format
#head(countData)

## meta data
cluster = scores.tmm.nmds_family2$apclusters 
#cluster[cluster == 1] = 0
#cluster[cluster == 2] = 1
group <- as.factor(cluster)


## build DESeq object
metaData<-data.frame(row.names=colnames(countData), group=group) 


dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = ~group)

## filter
dds <- dds[rowSums(counts(dds)) > 0,]

## normalization
dds <- estimateSizeFactors(dds)

## estimate dispersion
dds<- estimateDispersions(dds)

## test the differential abundance
dds$group <- relevel(dds$group, "1")
#dds$group <- factor(dds$group, levels = c("NonSmoker", "Smoker"))
dds <- DESeq(dds)
res <- results(dds)
#sink(file = "deseq2_results.txt")
res
#sink(file = NULL)



sum(res$pvalue < 0.01, na.rm=TRUE) ## 78 genera were found differently abundant on 0.01-level without correction
##Next we would like to see how many genera would show different abundance with Benjamini-Hochberg correction.
#Default values for FDR correction in DESeq2 is 10%.
table(res[,"padj"] < 0.01) ## 29

#The following table shows the strongest down-regulated genera in cluster 2 comparing to cluster 1:
res_Sig <- res[which(res$padj < 0.01 ),]
#sink(file = "deseq2_downreg.txt")
head(res_Sig[order(res_Sig$log2FoldChange),])
#sink(file = NULL)

#sink(file = "deseq2_upreg.txt")
head(res_Sig[order(-res_Sig$log2FoldChange),]) ##overubandant genera
#sink(file = NULL)

# res_Sig = res_Sig[order(-res_Sig$log2FoldChange),]
# res_ap_genera = as.data.frame(res_Sig[, c(2, 6)])
# rownames(res_ap_genera) = gsub("_", ".", rownames(res_ap_genera))
# colnames(res_ap_genera) = c("log2FoldChange", "BH adjusted p-value")
# write.csv(res_ap_genera, file="DA_ap_BH.csv", quote=FALSE)



#sink(file = "da_ap.txt")
#res_Sig[order(-res_Sig$log2FoldChange),] ##overubandant genera
#sink(file = NULL)
#res_sig_rounded = res_Sig
#res_sig_rounded = apply(res_sig_rounded, 2, function(x) round(x, 5))
#write.csv(as.data.frame(res_sig_rounded), file="DA_AP_rd.csv")


#setEPS()
#postscript("thr.eps", height = 5, width = 9)
plot(metadata(res)$filterNumRej, type="b", ylab="number of rejections", xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red") 
abline(v=metadata(res)$filterTheta)
#dev.off()

#setEPS()
#postscript("MA.eps", height = 5, width = 9)
plotMA(res)
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
#dev.off()


#setEPS()
#postscript("disp_est.eps", height = 5, width = 9)
plotDispEsts(dds, ylim = c(1e-2, 1e3))
#dev.off()

library("EnhancedVolcano")

## test the differential abundance
#dds$group <- relevel(dds$group, "1")
##dds$group <- factor(dds$group, levels = c("NonSmoker", "Smoker"))
#dds <- DESeq(dds)
#res <- results(dds)


EnhancedVolcano(res, 
                lab = rownames(res), 
                x = 'log2FoldChange', 
                y = 'padj', 
                pCutoff = 0.01,
                ylim = c(0,6),
                legendLabSize = 6.5,
                legendIconSize = 3.0,
                legendLabels=c('NS',
                               expression(paste('Log'[2], 'FC')),
                               'p-value',
                               expression(paste('p-v. & Log'[2], 'FC'))))


EnhancedVolcano(res, 
                lab = rownames(res), 
                x = 'log2FoldChange', 
                y = 'padj',
                pCutoff = 0.01,
                xlim = c(-5,5),
                ylim = c(0,6),
                legendLabSize = 6.5,
                legendIconSize = 3.0,
                legendLabels=c('NS',
                               expression(paste('Log'[2], 'FC')),
                               'p-value',
                               expression(paste('p-v. & Log'[2], 'FC'))))

vst <-varianceStabilizingTransformation(dds)
topVarGenes <- head(order(rowVars(assay(vst)), decreasing=TRUE), 10)
assay(vst)[topVarGenes,0]

library("gplots")
library("RColorBrewer")
library("genefilter")
library("SummarizedExperiment")
par("mar")
heatmap.2(assay(vst)[topVarGenes,], scale="row", trace="none", dendrogram="none", col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

library("gap")
#setEPS()
#postscript("p_values_hist_na_removed", height = 5, width = 9)
hist(res$pvalue, breaks = 0:20/20, col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Independent Filtering", xlab = "p-values")
#dev.off()

#setEPS()
#postscript("p-value_low_means_removed.eps", height = 5, width = 9)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Filter of Low Mean Normalized Counts", xlab = "p-values")
#dev.off()




library("fdrtool")
res <- res[ !is.na(res$padj),]
res <- res[ !is.na(res$pvalue),]
res <- res[, -which(names(res) == "padj")]

res_fdr <- fdrtool(res$stat, statistic= "normal", plot = T)
#res_fdr$param[1, "sd"]
res[,"padj"] <- p.adjust(res_fdr$pval, method = "BH")

#setEPS()
#postscript("p_values_hist_corr", height = 5, width = 9)
hist(res_fdr$pval,breaks = 0:20/20,
     col = "grey50", border = "white", main = "Cluster 1 vs. Cluster 2, corrected p-values distribution", xlab = "Corrected p-values")
#dev.off()


setEPS()
postscript("pv_hist_line.eps", height = 5, width = 9)
hist(res_fdr$pval, breaks = 0:20/20, col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Histogram of Corrected p-values", xlab = "Corrected p-values")
abline(h=41, col=c("dodgerblue"))
dev.off()



# After adjusting p-values with BH with fdrtool(), it appeared that 584 genera were differently abundant in cluster 2.
table(res[,"padj"] < 0.1)



## meta data
cluster2 = scores.tmm.nmds_family2$cluster 
#cluster[cluster == 2] = 1
#cluster2
group2 <- as.factor(cluster2)


## build DESeq object
metaData.alt2<-data.frame(row.names=colnames(countData),group=group2) 
#head(metaData)


dds_alt2 <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData.alt2, 
                              design = ~group)

## filter
dds_alt2 <- dds_alt2[rowSums(counts(dds_alt2)) > 0,]

## normalization
dds_alt2 <- estimateSizeFactors(dds_alt2)

## estimate dispersion
dds_alt2 <- estimateDispersions(dds_alt2)

## test the differential abundance
dds_alt2$group <- relevel(dds_alt2$group, "1")
#dds$group <- factor(dds$group, levels = c("NonSmoker", "Smoker"))
dds_alt2 <- DESeq(dds_alt2)
res_alt2 <- results(dds_alt2)

#sink(file = "deseq2_results_alt2.txt")
res_alt2
#sink(file = NULL)


sum(res_alt2$pvalue < 0.01, na.rm=TRUE) ## 96 genera were found differently abundant on 0.01-level without correction
##Next we would like to see how many genera would show different abundance with Benjamini-Hochberg correction.
#Default values for FDR correction in DESeq2 is 10%.
table(res_alt2[,"padj"] < 0.01) ## 65


#The following table shows the strongest down-regulated genera in cluster 2 comparing to cluster 1:
res_Sig_alt2 <- res_alt2[which(res_alt2$padj < 0.01 ),]
#sink(file = "deseq2_results_down_alt2.txt")
head(res_Sig_alt2[order(res_Sig_alt2$log2FoldChange),])
#sink(file = NULL)

#sink(file = "deseq2_results_up_alt2.txt")
head(res_Sig_alt2[order(-res_Sig_alt2$log2FoldChange),])
#sink(file = NULL)

# res_Sig_hc = res_Sig_alt2[order(-res_Sig_alt2$log2FoldChange),]
# res_hc_genera = as.data.frame(res_Sig_hc[, c(2, 6)])
# rownames(res_hc_genera) = gsub("_", ".", rownames(res_hc_genera))
# colnames(res_hc_genera) = c("log2FoldChange", "BH adjusted p-value")
# write.csv(res_hc_genera, file="DA_hc_BH.csv", quote=FALSE)

#sink(file = "da_hc.txt")
#res_Sig_alt2[order(-res_Sig_alt2$log2FoldChange),] ##overubandant genera
#sink(file = NULL)
#res_sig_rounded_alt2 = res_Sig_alt2
#res_sig_rounded_alt2 = apply(res_sig_rounded, 2, function(x) round(x, 5))
#write.csv(as.data.frame(res_Sig_alt2), file="da_hc.csv")
#dim(res_Sig_alt2)


#setEPS()
#postscript("ma_plot_HC.eps", height = 5, width = 9)
topGene_alt2 <- rownames(res_alt2)[which.min(res_alt2$padj)]
plotMA(res_alt2)
with(res[topGene_alt2, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene_alt2, pos=2, col="dodgerblue")
})
#dev.off()

#setEPS()
#postscript("disp_est_HC.eps", height = 5, width = 9)
#plotDispEsts(dds_alt2, ylim = c(1e-2, 1e3))
#dev.off()

#topGene2_alt <- rownames(res_alt2)[which.min(res_alt2$padj)]
#plotCounts(dds_alt2, gene = topGene2_alt, intgroup = c("group"))

#topGene2 <- rownames(res_alt2)[which.max(res_alt2$log2FoldChange)]
#plotCounts(dds_alt2, gene = topGene2, intgroup = c("group"))

 
library("ggbeeswarm")
geneCounts <- plotCounts(dds_alt2, gene = topGene2, intgroup = c("group"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = group, y = count)) + 
         scale_y_log10() + 
         geom_beeswarm(cex = 3)
ggplot(geneCounts, aes(x = group, y = count, group = group)) +
  scale_y_log10() + 
  geom_point(size = 3) + 
  geom_line()

#Diagnostics with p-value

#setEPS()
#postscript("pv_alt2_indep_filt.eps", height = 5, width = 9)
hist(res_alt2$pvalue, breaks = 0:20/20, col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Independent Filtering", xlab = "p-values")
#dev.off()

#setEPS()
#postscript("pv_alt2_low_removed.eps", height = 5, width = 9)
hist(res_alt2$pvalue[res_alt2$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Filter of Low Mean Normalized Counts", xlab = "p-values")
#dev.off()

hist(res_alt2$pvalue, breaks=20, col="grey", main = "1st cluster vs. 2nd cluster", xlab = "p-values")


library("fdrtool")

#res
dim(res_alt2)
res_alt2 <- res_alt2[ !is.na(res_alt2$padj),]
res_alt2 <- res_alt2[ !is.na(res_alt2$pvalue),]

res_alt2 <- res_alt2[, -which(names(res_alt2) == "padj")]
#res_fdr$param[1, "sd"]

res_fdr_alt2 <- fdrtool(res_alt2$stat, statistic= "normal", plot = T)
res_fdr_alt2$param[1, "sd"]
res_alt2[,"padj"] <- p.adjust(res_fdr_alt2$pval, method = "BH")


#setEPS()
#postscript("pv_hist_line_alt2.eps", height = 5, width = 9)
hist(res_fdr_alt2$pval, breaks = 0:20/20, col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Histogram of Corrected p-values", xlab = "Corrected p-values")
abline(h=41, col=c("dodgerblue"))
#dev.off()


# After adjusting p-values with BH with fdrtool(), it appeared that 54 genera were fount to be differently abundant in cluster 2.
sum(res_alt2[,"padj"] < 0.01)
res_Sig_hc_corr = res_alt2[which(res_alt2$padj < 0.01 ),]

# res_Sig_hc_corr = res_Sig_hc_corr[order(-res_Sig_hc_corr$log2FoldChange),]
# res_hc_genera_corr = as.data.frame(res_Sig_hc_corr[, c(2, 6)])
# rownames(res_hc_genera_corr) = gsub("_", ".", rownames(res_hc_genera_corr))
# colnames(res_hc_genera_corr) = c("log2FoldChange", "BH adjusted p-value")
# write.csv(res_hc_genera_corr, file="DA_hc_corr.csv", quote=FALSE)

## Volcano Plot 
EnhancedVolcano(res_alt2,
    lab = rownames(res_alt2),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.001)
EnhancedVolcano(res_alt2,
    lab = rownames(res_alt2),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.01)

plot(res_alt2$baseMean+1, -log10(res_alt2$padj),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     cex=.4, col=rgb(0,0,0,.3))

library("edgeR")
# data and groups for creating DGEList type object
counts = as.data.frame(t(raw_counts_family))
group = df_na$farm

y = DGEList(counts=counts,group=group)
y_full <- y
y_full2 = y_full

## keeping genera which have count of more than 5 in at least 2 samples
keep <- rowSums(y$counts>5) #>= 2
y <- y[keep,]

## keeping genera which have cpm of more than 100 in at least 2 samples
keep_cpm <- rowSums(cpm(y) > 100) >= 2
y_cpm <- y_full[keep_cpm,]


## reset library size
y$samples$lib.size <- colSums(y$counts)


## calculate TMM normalization factors fro effective library size
y <- calcNormFactors(y)

## calculate effective library sizes
y$samples$lib.size*y$samples$norm.factors

#(1) estimate the common dispersion:
y1 <- estimateCommonDisp(y, verbose=T)

#(2) Fit a trended model to get a genus-wise dispersion.
#A plot the tag-wise biological coefficient of variation (square root of #dispersions) against log2-CPM:

y1 <- estimateTagwiseDisp(y1)
plotBCV(y1)

#(3) Fit a generalized linear model to estimate the genus-wise dispersion.
#use a generalized linear model to estimate the dispersion
cluster = as.numeric(scores.gmds$apclusters)
#cluster[cluster == 1] = 0
#cluster[cluster == 2] = 1
#cluster
design <- model.matrix(~cluster)
rownames(design) <- colnames(y)
#design

library("statmod")
y2 <- estimateDisp(y, design, robust=TRUE)
#y2$common.dispersion
plotBCV(y2)

# Estimate and visualize the quasi-likelihood dispersions (based on the fitted quasi-likelihood GLM):
fit <- glmQLFit(y2, design, robust=TRUE)
plotQLDisp(fit)


## Test differential abundances
## 2 clusters are taken from finite apclustering
y_clusters = DGEList(counts=counts,group=as.factor(scores.gmds$apclusters))
y1_clusters <- estimateCommonDisp(y_clusters, verbose=T)
y1_clusters <- estimateTagwiseDisp(y1_clusters)

#The following table includes the list of top-10 differently abundant genera (BH corrected p-value = column "FDR"):
et <- exactTest(y1_clusters, pair = c("1", "2"))
topTags(et)

design2 = model.matrix(~cluster)
rownames(design2) <- colnames(y)
#design2

fit <- glmQLFit(y1, design2)
qlf <- glmQLFTest(fit, contrast=c(-1,1))
topTags(qlf)

FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05) ## 515

qlf_lrt <- glmLRT(fit, contrast=c(-1,1))
topTags(qlf_lrt)

#another design
y$samples$cluster = scores.gmds$apclusters
#y$samples$cluster[y$samples$cluster == "1"] = 0
#y$samples$cluster[y$samples$cluster == "2"] = 1
y$samples$cluster = as.factor(y$samples$cluster)

design1 <- model.matrix(~0+cluster, data=y$samples)
colnames(design1) <- levels(as.factor(y1$clust))
design1

fit1 <- glmQLFit(y1, design1)
qlf1 <- glmQLFTest(fit1, contrast=c(-1,1))
topTags(qlf1)

FDR1 <- p.adjust(qlf1$table$PValue, method="BH")
sum(FDR1 < 0.05) ## 287
topTags(qlf1, n=15)

qlf1_lrt <- glmLRT(fit1, contrast=c(-1,1))
topTags(qlf1_lrt)

da = decideTestsDGE(et, p.value = 0.1)
da_OTUs = rownames(y1)[as.logical(da)]
plotSmear(et , de.tags = da_OTUs, cex = 0.5)
abline(h = c(-2, 2), col = "blue") 

tab = data.frame(logFC = et$table[, 1], negLogPval = -log10(et$table[, 3] ))
head(tab)
## plotting

par(mar = c(5, 4, 4, 4))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))

# Log2 fold change and p-value cutoffs
lfc = 2
pval = 0.1

# Selecting OTUs of interest
sig_OTUs = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
# Identifying the selected OTUs
points(tab[sig_OTUs, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
