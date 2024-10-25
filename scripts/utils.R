library(EnhancedVolcano)

customVolcano <- function(top.table, phenotype, DEGs){
  volcano.col <- ifelse(top.table$logFC > 1, 'red', 'black')
  names(volcano.col)[volcano.col == 'red'] <- 'Up-regulated'
  names(volcano.col)[volcano.col == 'black'] <- 'Non-significant'
  EnhancedVolcano(top.table,
                  title = paste0(phenotype, " vs all other phenotypes"),
                  subtitle = NULL,
                  lab = rownames(top.table),
                  selectLab = DEGs,
                  x = 'logFC',
                  y = 'adj.P.Val', #P.Value
                  pCutoff = 0.01, #0.05
                  drawConnectors = TRUE,
                  widthConnectors = 0.75,
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  colCustom = volcano.col)
}

wilcoxontestBarplot <- function(paris.to.plot, signatures, thres.barplot){
  barplot.pval <- data.frame()
  # Loop over each pair
  for (pair in pairs.to.plot) {
    # Get the Phenotype and Signature from the pair
    phenotype <- pair[2]
    selected_signature <- pair[1]
    
    # Subset the dataframe for the current Phenotype
    df_phenotype <- thres.barplot[thres.barplot$Phenotype == phenotype,]
    
    # Subset the dataframe for the selected Signature
    df_selected <- df_phenotype[df_phenotype$Signature == selected_signature,]
    
    # Loop over each other Signature
    for (signature in signatures) {
      if (signature != selected_signature) {
        # Subset the dataframe for the current Signature
        df_other <- df_phenotype[df_phenotype$Signature == signature,]
        
        # Perform the Wilcoxon test
        test_result <- wilcox.test(df_selected$Expression, df_other$Expression)
        
        # Store the p-value in the result dataframe
        barplot.pval <- rbind(barplot.pval, data.frame(Phenotype = phenotype, Signature = signature, P_Value = test_result$p.value))
      }
    }
  }
  
  # Reshape the result dataframe to the desired format
  barplot.pval <- barplot.pval %>% spread(key = Signature, value = P_Value)
  return(barplot.pval)
}

