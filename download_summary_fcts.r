meta_data_download <- function(){
  meta_data <- as.data.frame(read_xlsx("~/Desktop/Thesis/Optibiom/Optibiom_metadata_perSample_AL10022020.xlsx", 2))
  meta_data[["Probensammlung"]] <- as.POSIXct(meta_data[["Probensammlung"]]*(60*60*24), origin = "1899-12-30")
  meta_data[meta_data==-99] <- NA
  return(meta_data)
}
 

asv_data_download <-function(){
  ## asv classification data
  asv = as.data.frame(read_xlsx("~/Desktop/Thesis/Optibiom/Optibiom_ASV_table_AL_10022020.xlsx", 2))
  asv.genus = as.factor(asv[,"Genus"])
  
  ## download asv count data
  asv_count <- as.data.frame(read_xlsx("~/Desktop/Thesis/asv/Optibiom_SchweinASV&taxonomy05_2020 short.xlsx"))
  asv_count_genus <- asv_count[,-c(1:5)]
  colnames(asv_count_genus) <-  sub("sample.", "", colnames(asv_count_genus))
  groups <- asv_count_genus$Genus
  sum_asv_count_genus <- as.data.frame(rowsum(x = asv_count_genus[,-1], group = groups))
  return(sum_asv_count_genus)
} 


overdisp_plot <-function(sum_asv_count_genus){
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
  #print(p)
  return(p)
}


summary_table <- function(sum_asv_count_genus){
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
  #write.csv(summary_genera, file="summary_gen.csv", quote=FALSE)
  return(summary_genera)
}

## microbiome data on phylum level
microbiome_data <- function(){
  data_microbiome <- as.data.frame(read_xlsx("~/Desktop/Thesis/Optibiom/Optibiom_tax_profile_AL_10022020.xlsx"), header = TRUE, row.names = 1) ##phylum level
  row.names(data_microbiome) <- meta_data[,1]
  data_microbiome = data_microbiome[,-1]
  data_microbiome$farm = meta_data$Farm
  return(data_microbiome)
}

## asv data (genus level)

genus_raw_and_tss_data <- function(){
  ## outputs raw microbiome counts on genus level
  asv_count <- as.data.frame(read_xlsx("~/Desktop/Thesis/asv/Optibiom_SchweinASV&taxonomy05_2020 short.xlsx"))
  asv_count_genus <- asv_count[,-c(1:5)]
  asv_count_genus$Genus <- as.factor(asv_count_genus$Genus)
  groups <- asv_count_genus$Genus
  sum_asv_count_genus <- as.data.frame(rowsum(x = asv_count_genus[,-1], group = groups))
  
  ## total sum normalization
  sum_asv_count_genus_norm <- decostand(x = sum_asv_count_genus, method="total", MARGIN = 2)
  
  sum_asv_count_genus_norm$Genus <- levels(asv_count_genus$Genus)
  rownames(sum_asv_count_genus_norm) = sum_asv_count_genus_norm$Genus
  sum_asv_count_genus_norm = sum_asv_count_genus_norm[,-803]
  
  sum_asv_count_genus_norm <- as.data.frame(t(sum_asv_count_genus_norm))
  sum_asv_count_genus_norm$farm = meta_data$Farm
  #write.csv(sum_asv_count_genus_norm, file="genera_norm.csv", quote=FALSE)
  return(list(sum_asv_count_genus, sum_asv_count_genus_norm))
}

