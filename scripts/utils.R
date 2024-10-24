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