differential_testing <- function(raw_counts_family, scores.tmm.nmds_family2, cluster_division){
  
  otu_tab = raw_counts_family
  countData <- as(otu_tab, "matrix")
  countData<-(t(countData)) #DESeq2 needs taxa(genes=rows) by samples(=columns) format
  
  ## meta data
  group <- as.factor(cluster_division)
  
  
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
  
  return(list(res, res_Sig))
}


pvalue_diagostics <- function(res, clustering_method){
  #setEPS()
  #postscript("p_values_hist_na_removed", height = 5, width = 9)
  hist(res$pvalue, breaks = 0:20/20, col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Independent Filtering", xlab = "p-values")

  #postscript("p-value_low_means_removed.eps", height = 5, width = 9)
  hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
       col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Filter of Low Mean Normalized Counts", xlab = "p-values")
  res <- res[ !is.na(res$padj),]
  res <- res[ !is.na(res$pvalue),]
  res <- res[, -which(names(res) == "padj")]
  
  res_fdr <- fdrtool(res$stat, statistic= "normal", plot = T)
  #res_fdr$param[1, "sd"]
  res[,"padj"] <- p.adjust(res_fdr$pval, method = "BH")
  
  #postscript("p_values_hist_corr", height = 5, width = 9)
  hist(res_fdr$pval,breaks = 0:20/20,
       col = "grey50", border = "white", main = "Cluster 1 vs. Cluster 2, corrected p-values distribution", xlab = "Corrected p-values")

  #postscript("pv_hist_line.eps", height = 5, width = 9)
  hist(res_fdr$pval, breaks = 0:20/20, col = "grey50", border = "white", main = "1st cluster vs. 2nd cluster: Histogram of Corrected p-values", xlab = "Corrected p-values")
  abline(h=41, col=c("dodgerblue"))
  #dev.off()
  if (clustering_method == 'ap'){
    return(table(res[,"padj"] < 0.01))
  } else if (clustering_method == 'hc') {
    res_hc_corr = res[which(res$padj < 0.01),]
    return(res_hc_corr)
  }

}


volcano_plots_ap <- function(res_ap){
  EnhancedVolcano(res_ap, 
                  lab = rownames(res_ap), 
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
  
  EnhancedVolcano(res_ap, 
                  lab = rownames(res_ap), 
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
}


volcano_plots_hc <- function(res_hc){
  EnhancedVolcano(res_hc, 
                  lab = rownames(res_hc), 
                  x = 'log2FoldChange', 
                  y = 'padj', 
                  pCutoff = 0.01,
                  legendLabSize = 6.5,
                  legendIconSize = 3.0,
                  legendLabels=c('NS',
                                 expression(paste('Log'[2], 'FC')),
                                 'p-value',
                                 expression(paste('p-v. & Log'[2], 'FC'))))
  
  
  EnhancedVolcano(res_hc, 
                  lab = rownames(res_hc), 
                  x = 'log2FoldChange', 
                  y = 'padj',
                  pCutoff = 0.01,
                  xlim = c(0,8),
                  legendLabSize = 6.5,
                  legendIconSize = 3.0,
                  legendLabels=c('NS',
                                 expression(paste('Log'[2], 'FC')),
                                 'p-value',
                                 expression(paste('p-v. & Log'[2], 'FC'))))
}