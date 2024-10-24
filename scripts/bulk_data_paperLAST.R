library(tidyverse)
library(limma)
library(edgeR)
library(pheatmap)
library(umap)
library(dendextend)
library(ggridges)
library(EnhancedVolcano)
library(RColorBrewer)

theme_set(theme_classic())

# root_path <- "M:/DER/LABOS/IMMUNODER/Antoine/Data/"

# Get the directory of the current script
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the root of your project
setwd(file.path(script_dir, ".."))

## Load raw data
library(nanostringr)
rcc_path <- "M:/DER/LABOS/IMMUNODER/Antoine/Data/NCOMMS-24-24396A/."
biopsy.mat <- read_rcc(rcc_path)
biopsy.mat <- as.data.frame(biopsy.mat$raw)
rownames(biopsy.mat) <- biopsy.mat$Name
biopsy.mat <- as.data.frame(t(biopsy.mat[,c(-1:-3)]))

biopsy.meta <- read.csv("data/biopsy_metaLAST17.csv", row.names = "X")
dupi.evol <- read.csv("data/testdata_dupi.csv", row.names = "X")

## Gene panels

hkg <- read.csv("data/housekeepinggenes.csv")$x
pos.probes <- colnames(biopsy.mat)[grepl("POS", colnames(biopsy.mat))]
neg.probes <- colnames(biopsy.mat)[grepl("NEG", colnames(biopsy.mat))][1:6]

signature.names <- c("Th1", "Th2", "Th17", "Neutro", "Macro", "Eosino", "IFN")
sign.col <- c("#FFCC66", "#99CCFF", "#CCFF99", "#FFCCFF", "#FF6666", "#6699CC", "#FFFF99")
names(sign.col) <- signature.names

der.sign <- read.csv("data/der_signature.csv")
nsf.sign <- read.csv("data/NSForest_markers_alldiseases.csv")

## Normalization
mpgm <- 1247
mhgm <- 1405
posgeomean <- apply(biopsy.mat[,pos.probes], 1, function(x) exp(mean(log(x))))
norm.factor <- mpgm/posgeomean
biopsy.norm.mat <- norm.factor*biopsy.mat
hkggeomean <- apply(biopsy.norm.mat[,hkg], 1, function(x) exp(mean(log(x))))
norm.factor <- mhgm/hkggeomean
biopsy.norm.mat <- norm.factor*biopsy.norm.mat
biopsy.norm.mat[biopsy.norm.mat < 1] <- 1
biopsy.norm.mat <- log2(biopsy.norm.mat)

## Scaling on sentinels

phenotypes <- c("LP", "AD", "PsO", "NeuD", "Wells", "CLE", "Healthy")
pheno.col <- c("#FFCC66", "#99CCFF", "#CCFF99", "#FFCCFF", "#6699CC", "#FFFF99", "#CCCCCC")
names(pheno.col) <- phenotypes
pheno2 <- c("BP", "DHR", "PRG", "Erythroderma", "undetermined rashes")
pheno.col2 <-  c("#66FFFF", "#66FF66","#FF9900", "#FF6666", "#FFFC33")
names(pheno.col2) <- pheno2
pheno.col.all <- c(pheno.col, pheno.col2)

sent1 <- biopsy.meta[biopsy.meta$Sentinel &
                        biopsy.meta$Phenotype %in% phenotypes, "Biopsy"]
sentinels.mat <- biopsy.norm.mat[sent1, der.sign$genes]
sentinels.meta <- biopsy.meta[sent1,]
sentinels.meta$Phenotype <- factor(sentinels.meta$Phenotype, levels = phenotypes)

testsample <- biopsy.meta[!(biopsy.meta$Biopsy %in% sent1), "Biopsy"]
test.mat <- biopsy.norm.mat[testsample, der.sign$genes]
test.meta <- biopsy.meta[testsample,]

m <- colMeans(sentinels.mat)
sd <- apply(sentinels.mat, 2, function(x) sd(x))
scaled.mat <- rbind(t((t(sentinels.mat) - m) / sd), 
                        t((t(test.mat) - m) / sd))
scaled.mat <- scaled.mat[order(row.names(scaled.mat)),]
scaled.mat <- scaled.mat[rownames(biopsy.meta),]
for (signature in levels(factor(der.sign$signature))){
  genes <- der.sign[der.sign$signature %in% signature, ]$genes
  if (signature %in% colnames(biopsy.meta)){
    biopsy.meta[, signature] <- apply(scaled.mat[, genes],1 , mean)
  } else{
    biopsy.meta[, ncol(biopsy.meta) + 1] <- apply(scaled.mat[, genes],1 , mean)
    colnames(biopsy.meta)[ncol(biopsy.meta)] <- signature
  }
}

thres <- c(0.56,0.33,0.44,0.53,0.5,1.05,0.86)
names(thres) <- signature.names
scaled.thres <- t(t(biopsy.meta[, signature.names]) - thres)
scaled.thres <- 1/(1+exp(3*(-scaled.thres)))

sentinels.meta[, signature.names] <- scaled.thres[rownames(sentinels.meta), signature.names]
sentinels.meta$dominant.module <- signature.names[max.col(sentinels.meta[,signature.names])]

#### FIG1A  DEGs

## Limma

mm <- model.matrix(~0 + sentinels.meta$Phenotype)
d <- t(biopsy.mat[sent1,])
d <- DGEList(counts = d, group = sentinels.meta$Phenotype)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
levels <- colnames(coef(fit))
contr <- diag(length(levels))
contr[contr == 0] <- -1/(length(levels)-1)
rownames(contr) <- levels
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

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


pheno.ind <- 1
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pdf("output/fig1LP.pdf", width = 7, height = 7)
customVolcano(top.table, phenotypes[pheno.ind], der.sign$genes[der.sign$signature %in% "Th1"])
dev.off()

pheno.ind <- pheno.ind +1
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pdf("output/fig1AD.pdf", width = 7, height = 7)
customVolcano(top.table, phenotypes[pheno.ind], der.sign$genes[der.sign$signature %in% "Th2"])
dev.off()

pheno.ind <- pheno.ind +1
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pdf("output/fig1PSO.pdf", width = 7, height = 7)
customVolcano(top.table, phenotypes[pheno.ind], der.sign$genes[der.sign$signature %in% "Th17"])
dev.off()

pheno.ind <- pheno.ind +1
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pdf("output/fig1ND.pdf", width = 7, height = 7)
customVolcano(top.table, phenotypes[pheno.ind], der.sign$genes[der.sign$signature %in% c("Neutro", "Macro")])
dev.off()

pheno.ind <- pheno.ind +1
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pdf("output/fig1WELLS.pdf", width = 7, height = 7)
customVolcano(top.table, phenotypes[pheno.ind], der.sign$genes[der.sign$signature %in% "Eosino"])
dev.off()

pheno.ind <- pheno.ind +1
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pdf("output/fig1LUPUS.pdf", width = 7, height = 7)
customVolcano(top.table, phenotypes[pheno.ind], der.sign$genes[der.sign$signature %in% "IFN"])
dev.off()

pheno.ind <- pheno.ind +1
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pdf("output/fig1HD.pdf", width = 7, height = 7)
customVolcano(top.table, phenotypes[pheno.ind], der.sign$genes)
dev.off()

#### FIG1B - UMAP

sentinels.umap <- umap(sentinels.mat, transform_seed=42)
sentinels.meta$UMAP1 <- sentinels.umap$layout[,1]
sentinels.meta$UMAP2 <- sentinels.umap$layout[,2]
pdf("output/fig1UMAP.pdf", width = 5, height = 3.5)
sentinels.meta %>% ggplot(aes(x=UMAP1, y=UMAP2, color=Phenotype))+
  geom_point(size=2, shape=21, col="black", aes(fill=Phenotype))+
  scale_fill_manual(values=pheno.col)+
  theme(legend.position = "top")+
  labs(x = "UMAP1",
       y = "UMAP2")
dev.off()

#### FIG1C - Heatmap

annot_col <- data.frame(Phenotype=sentinels.meta$Phenotype)
rownames(annot_col) <- rownames(sentinels.mat)

pdf("output/fig1Heatmap.pdf", width = 15, height = 15)
pheatmap <- pheatmap(t(na.omit(sentinels.mat)), cluster_rows = FALSE, scale = "row", 
                     clustering_distance_cols = "correlation", 
                     clustering_method = "complete", 
                     cellheight = 8, cellwidth = 8,fontsize = 12,
                     fontsize_col =10, fontsize_row = 8,
                     breaks = seq(-2.5, 2.5, by = 0.05), border_color = "black" ,
                     color = colorRampPalette(rev(c("red", "black", "green")))(100),
                     annotation_col = annot_col,
                     annotation_colors = list(Phenotype=pheno.col))
dev.off()

#### FIG1D - Cluster

dist <- as.dist(1 - cor(na.omit(t(scale(sentinels.mat))), method = "pearson"))
dend <- hclust(dist, method = "complete")
k <- length(unique(sentinels.meta$Phenotype))
dend <- as.dendrogram(dend) 
labels(dend) <- "" 
pdf("output/fig1Sign.pdf", width = 7, height = 4)
plot(dend)
colored_bars(pheno.col[sentinels.meta$Phenotype], dend  = dend)
rect.dendrogram(k=k, border = 8, lty = 5, lwd = 1, tree = dend)
dev.off()

sentinels.meta$cluster <- as.factor(cutree(dend, k=k))
table <- table(sentinels.meta[, c("cluster", "Phenotype")])
levels(sentinels.meta$cluster) <- colnames(table)[apply(table, 1, which.max)]
sentinels.meta$cluster <- factor(sentinels.meta$cluster, levels=levels(sentinels.meta$Phenotype))
perf.sent <- mltest::ml_test(sentinels.meta$cluster, sentinels.meta$Phenotype)
as.data.frame(perf.sent[c("precision", "recall", "specificity")])
FM_index(sentinels.meta$cluster, sentinels.meta$Phenotype, assume_sorted_vectors = TRUE)

sent.all.mat <- biopsy.norm.mat[sent1,]

dist <- as.dist(1 - cor(na.omit(t(scale(sent.all.mat))), method = "pearson"))
dend <- hclust(dist, method = "complete")
k <- length(unique(sentinels.meta$Phenotype))
dend <- as.dendrogram(dend) 
labels(dend) <- "" 
pdf("output/fig1All.pdf", width = 7, height = 4)
plot(dend)
colored_bars(pheno.col[sentinels.meta$Phenotype], dend  = dend)
rect.dendrogram(k=k, border = 8, lty = 5, lwd = 1, tree = dend)
dev.off()

sentinels.meta$clust.all <- as.factor(cutree(dend, k=k))
table <- table(sentinels.meta[, c("clust.all", "Phenotype")])
levels(sentinels.meta$clust.all) <- colnames(table)[apply(table, 1, which.max)]
sentinels.meta$clust.all <- factor(sentinels.meta$clust.all, levels=levels(sentinels.meta$Phenotype))
perf.sent <- mltest::ml_test(sentinels.meta$clust.all, sentinels.meta$Phenotype)
as.data.frame(perf.sent[c("precision", "recall", "specificity")])
FM_index(sentinels.meta$clust.all, sentinels.meta$Phenotype, assume_sorted_vectors = TRUE)

sent.nsf.mat <- biopsy.norm.mat[sent1, nsf.sign$markerGene]

dist <- as.dist(1 - cor(na.omit(t(scale(sent.nsf.mat))), method = "pearson"))
dend <- hclust(dist, method = "complete")
k <- length(unique(sentinels.meta$Phenotype))
dend <- as.dendrogram(dend) 
labels(dend) <- ""
pdf("output/fig1NSF.pdf", width = 7, height = 4)
plot(dend)
colored_bars(pheno.col[sentinels.meta$Phenotype], dend  = dend)
rect.dendrogram(k=k, border = 8, lty = 5, lwd = 1, tree = dend)
dev.off()

sentinels.meta$clust.nsf <- as.factor(cutree(dend, k=k))
table <- table(sentinels.meta[, c("clust.nsf", "Phenotype")])
levels(sentinels.meta$clust.nsf) <- colnames(table)[apply(table, 1, which.max)]
sentinels.meta$clust.nsf <- factor(sentinels.meta$clust.nsf, levels=levels(sentinels.meta$Phenotype))
perf.sent <- mltest::ml_test(sentinels.meta$clust.nsf, sentinels.meta$Phenotype)
as.data.frame(perf.sent[c("precision", "recall", "specificity")])
FM_index(sentinels.meta$clust.nsf, sentinels.meta$Phenotype, assume_sorted_vectors = TRUE)

#### FIG2A - Density plot sentinels/signature

violin <- biopsy.meta[sent1,]
violin <- violin %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
violin$thres <- thres[violin$Signature]
pdf("output/fig2Density.pdf", width = 6, height = 18)
ggplot(violin, aes(x=Expression,y=Phenotype, fill=Signature))+
  geom_density_ridges(jittered_points = TRUE) +
  geom_vline(aes(xintercept = thres),data=violin) +
  facet_wrap(~Signature, nrow = 7) +
  theme(legend.position = "top")+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG2B - Boxplot logit transform sentinels/signature

violin <- as.data.frame(scaled.thres[sent1,])
violin$Phenotype <- biopsy.meta[sent1, "Phenotype"]
violin <- violin %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig2Boxplot.pdf", width = 6, height = 18)
ggplot(violin, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2)) +
  facet_wrap(~Phenotype, nrow = 7)+
  theme(legend.position = "top")+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG2C - Boxplot logit transform challengers/signature

test.meta <- as.data.frame(scaled.thres[biopsy.meta$Challenger,])
test.meta$Phenotype <- biopsy.meta[biopsy.meta$Challenger, "Phenotype"]
test.meta$dominant.module <- signature.names[max.col(test.meta[,signature.names])]
violin <- test.meta %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig2Boxplot2.pdf", width = 6, height = 18)
ggplot(violin, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2)) +
  facet_wrap(~Phenotype, nrow = 7)+
  theme(legend.position = "top")+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG2D - UMAP

test.mat <- biopsy.norm.mat[biopsy.meta$Challenger, der.sign$genes]
test.umap <- as.data.frame(predict(sentinels.umap, test.mat))
test.meta$UMAP1 <- test.umap[rownames(test.meta),1]
test.meta$UMAP2 <- test.umap[rownames(test.meta),2]
plot.umap <- rbind(test.meta[, c("Phenotype", "UMAP1", "UMAP2")], 
                   sentinels.meta[, c("Phenotype", "UMAP1", "UMAP2")])
plot.umap[rownames(plot.umap) %in% rownames(sentinels.meta), "Biospy"] <- "Sentinel"
plot.umap[is.na(plot.umap$Biospy), "Biospy"] <- "Challenger"
pdf("output/fig2UMAP.pdf", width = 5, height = 3.5)
plot.umap %>% ggplot(aes(x=UMAP1, y=UMAP2))+
  geom_point(size=2, shape=21, aes(fill=Phenotype, col=Biospy))+
  scale_fill_manual(values=pheno.col)+
  scale_color_manual(values=c("black", NA))+
  theme(legend.position = "top")+
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

test.meta$cluster <- test.meta$Phenotype
test.meta[c("CLE_003", "AD_047"), "cluster"] <- c("LP", "Healthy") #"SW_002","PsO"
test.meta$cluster <- factor(test.meta$cluster, levels=levels(sentinels.meta$Phenotype))
test.meta$Phenotype <- factor(test.meta$Phenotype, levels=levels(sentinels.meta$Phenotype))
perf.chal <- mltest::ml_test(test.meta$cluster, test.meta$Phenotype)
as.data.frame(perf.chal[c("precision", "recall", "specificity")])
FM_index(test.meta$cluster, test.meta$Phenotype, assume_sorted_vectors = TRUE)

#### SuppTable1 - Significant pathway sentinels

sent.thres <- as.data.frame(scaled.thres[sent1,])
sent.thres$Phenotype <- biopsy.meta[sent1, "Phenotype"]
sent.thres <- sent.thres %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)

# Create an empty dataframe to store the results
sent.pval <- data.frame()

# Get unique Phenotypes and Signatures
# phenotypes <- unique(sent.thres$Phenotype)
signatures <- unique(sent.thres$Signature)
pairs.to.plot <- list(c("Th1", "LP"),
                      c("Th2", "AD"),
                      c("Th17", "PsO"),
                      c("Neutro", "NeuD"),
                      c("Eosino", "Wells"),
                      c("IFN", "CLE"))

# Loop over each pair
for (pair in pairs.to.plot) {
  # Get the Phenotype and Signature from the pair
  phenotype <- pair[2]
  selected_signature <- pair[1]
  
  # Subset the dataframe for the current Phenotype
  df_phenotype <- sent.thres[sent.thres$Phenotype == phenotype,]
  
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
      sent.pval <- rbind(sent.pval, data.frame(Phenotype = phenotype, Signature = signature, P_Value = test_result$p.value))
    }
  }
}

# Reshape the result dataframe to the desired format
sent.pval <- sent.pval %>% spread(key = Signature, value = P_Value)

#### Significant pathway challenger

chal.thres <- test.meta
chal.thres$Phenotype <- biopsy.meta[rownames(chal.thres), "Phenotype"]
chal.thres <- chal.thres %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)

# Create an empty dataframe to store the results
chal.pval <- data.frame()

# Get unique Phenotypes and Signatures
# phenotypes <- unique(chal.thres$Phenotype)
signatures <- unique(chal.thres$Signature)
pairs.to.plot <- list(c("Th1", "LP"),
                      c("Th2", "AD"),
                      c("Th17", "PsO"),
                      c("Neutro", "NeuD"),
                      c("IFN", "CLE"))

# Loop over each pair
for (pair in pairs.to.plot) {
  # Get the Phenotype and Signature from the pair
  phenotype <- pair[2]
  selected_signature <- pair[1]
  
  # Subset the dataframe for the current Phenotype
  df_phenotype <- chal.thres[chal.thres$Phenotype == phenotype,]
  
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
      chal.pval <- rbind(chal.pval, data.frame(Phenotype = phenotype, Signature = signature, P_Value = test_result$p.value))
    }
  }
}

# Reshape the result dataframe to the desired format
chal.pval <- chal.pval %>% spread(key = Signature, value = P_Value)

### FIG3 - Boxplot logit transform PB, TOX-MP, PRG

pheno2 <- c(phenotypes, "BP", "DHR")
sent2 <- biopsy.meta[biopsy.meta$Sentinel &
                       biopsy.meta$Phenotype %in% pheno2, "Biopsy"]
subset.mat <- scaled.mat[sent2, der.sign$genes]
subset.meta <- as.data.frame(scaled.thres[rownames(subset.mat),])
subset.meta$Phenotype <- biopsy.meta[sent2, "Phenotype"]
subset.meta$dominant.module <- signature.names[max.col(subset.meta[,signature.names])]

subset.meta$threshold <- 0.7 * apply(subset.meta, 1, function(row){
  column_name <- row["dominant.module"]
  return(as.numeric(row[column_name]))
  })

subset.meta$codominant <- apply(subset.meta, 1, function(row) {
  threshold <- max(as.numeric(row["threshold"]), 0.3)
  return(paste(names(which(row[signature.names] > threshold)), collapse = ","))
})

# subset.meta$Classification <- as.factor(Mclust(as.data.frame(biopsy.meta[rownames(subset.meta), signature.names]), 3)$classification)
violin <- subset.meta
violin <- violin %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig3others.pdf", width = 6, height = 10)
ggplot(violin, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "top")+
  facet_wrap(~Phenotype, nrow = 4)+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG3 - UMAP

sent.other.umap <- umap(subset.mat)
sent.other.meta <- subset.meta
sent.other.meta$UMAP1 <- sent.other.umap$layout[,1]
sent.other.meta$UMAP2 <- sent.other.umap$layout[,2]
pdf("output/fig3otherUMAP.pdf", width = 5, height = 3.5)
sent.other.meta %>% ggplot(aes(x=UMAP1, y=UMAP2))+
  geom_point(size=2, shape=21, aes(fill=Phenotype))+
  theme(legend.position = "top")+
  scale_fill_manual(values=pheno.col.all)
dev.off()

#### FIG3 - Heatmap PB/TOX/

pheno2 <- c("BP", "DHR")
sent2 <- biopsy.meta[biopsy.meta$Sentinel & biopsy.meta$Phenotype %in% pheno2,
                     "Biopsy"]

merge.mat <- scaled.mat[c(sent1, sent2), ]
merge.meta <- biopsy.meta[c(sent1, sent2), ]
merge.meta$Phenotype <- as.character(merge.meta$Phenotype)

annot_col <- data.frame(Phenotype=merge.meta$Phenotype)
rownames(annot_col) <- rownames(merge.mat)

pdf("output/fig3othersHM.pdf", width = 18, height = 15)
pheatmap <- pheatmap(t(na.omit(merge.mat)), cluster_rows = FALSE, 
                     clustering_distance_cols = "correlation", 
                     clustering_method = "complete", 
                     cellheight = 8, cellwidth = 8,fontsize = 12,
                     fontsize_col =10, fontsize_row = 8,
                     breaks = seq(-2.5, 2.5, by = 0.05), border_color = "black" ,
                     color = colorRampPalette(rev(c("red", "black", "green")))(100),
                     annotation_col = annot_col,
                     annotation_colors = list(Phenotype=pheno.col.all)
                     )
dev.off()

### Fig3 - Hierachical clustering

pheno2 <- c("BP", "DHR", "PRG")
sent2 <- biopsy.meta[biopsy.meta$Sentinel & biopsy.meta$Phenotype %in% pheno2,
                     "Biopsy"]

merge.mat <- scaled.mat[c(sent1, sent2), ]
merge.meta <- biopsy.meta[c(sent1, sent2), ]
merge.meta$Phenotype <- as.character(merge.meta$Phenotype)
merge.col <- pheno.col.all[merge.meta$Phenotype]

dist <- as.dist(1 - cor(na.omit(t(merge.mat)), method = "pearson"))
dend <- hclust(dist, method = "complete")
k <- length(unique(merge.meta$Phenotype))
dend <- as.dendrogram(dend) 
labels(dend) <- "" 
pdf("output/fig3Sign.pdf", width = 7, height = 4)
plot(dend)
colored_bars(merge.col, dend  = dend)
rect.dendrogram(k=k, border = 8, lty = 5, lwd = 1, tree = dend)
dev.off()

merge.meta$cluster <- as.factor(cutree(dend, k=k))
FM_index(merge.meta$cluster, merge.meta$Phenotype, assume_sorted_vectors = TRUE)

merge.all.mat <- biopsy.norm.mat[rownames(merge.mat),]

dist <- as.dist(1 - cor(na.omit(t(scale(merge.all.mat))), method = "pearson"))
dend <- hclust(dist, method = "complete")
k <- length(unique(merge.meta$Phenotype))
dend <- as.dendrogram(dend) 
labels(dend) <- "" 
pdf("output/fig3bAll.pdf", width = 7, height = 4)
plot(dend)
colored_bars(merge.col, dend  = dend)
rect.dendrogram(k=k, border = 8, lty = 5, lwd = 1, tree = dend)
dev.off()

merge.meta$clust.all <- as.factor(cutree(dend, k=k))
FM_index(merge.meta$clust.all, merge.meta$Phenotype, assume_sorted_vectors = TRUE)

merge.nsf.mat <- biopsy.norm.mat[rownames(merge.mat), nsf.sign$markerGene]

dist <- as.dist(1 - cor(na.omit(t(scale(merge.nsf.mat))), method = "pearson"))
dend <- hclust(dist, method = "complete")
k <- length(unique(merge.meta$Phenotype))
dend <- as.dendrogram(dend) 
labels(dend) <- ""
pdf("output/fig3NSF.pdf", width = 7, height = 4)
plot(dend)
colored_bars(merge.col, dend  = dend)
rect.dendrogram(k=k, border = 8, lty = 5, lwd = 1, tree = dend)
dev.off()

merge.meta$clust.nsf <- as.factor(cutree(dend, k=k))
FM_index(merge.meta$clust.nsf, merge.meta$Phenotype, assume_sorted_vectors = TRUE)

#### FIG4 Erythroderma & undetermined rashes

# Erythrodermia

subset.mat <- scaled.mat[biopsy.meta$erythro, der.sign$genes]
subset.meta <- as.data.frame(scaled.thres[biopsy.meta$erythro,])
subset.meta$Phenotype <- "Erythroderma"
subset.meta$dominant.module <- signature.names[max.col(subset.meta[,signature.names])]
# subset.meta[erythro$Patient, 9:20] <- erythro[,-4]

pheno.erythro <- c("AD", "PsO", "DHR")
sent.erythro <- biopsy.meta[biopsy.meta$Sentinel & 
                              biopsy.meta$Phenotype %in% pheno.erythro,
                            "Biopsy"]

merge.mat <- scaled.mat[c(sent.erythro, rownames(subset.mat)), ]
merge.meta <- biopsy.meta[rownames(merge.mat),]
merge.meta$Phenotype <- as.character(merge.meta$Phenotype)
merge.meta[rownames(subset.mat), "Phenotype"] <- "Erythroderma"

sent.erythro.umap <- umap(merge.mat[sent.erythro,])
sent.erythro.meta <- merge.meta[sent.erythro,]
sent.erythro.meta$UMAP1 <- sent.erythro.umap$layout[,1]
sent.erythro.meta$UMAP2 <- sent.erythro.umap$layout[,2]

subset.umap <- predict(sent.erythro.umap, subset.mat)
subset.meta$UMAP1 <- subset.umap[,1]
subset.meta$UMAP2 <- subset.umap[,2]
plot.umap <- rbind(subset.meta[, c("Phenotype", "UMAP1", "UMAP2")], 
                   sent.erythro.meta[, c("Phenotype", "UMAP1", "UMAP2")])
plot.umap$Sentinels <- rownames(plot.umap) %in% rownames(sent.erythro.meta)
pdf("output/fig4erythroUMAP.pdf", width = 5, height = 3.5)
plot.umap %>% ggplot(aes(x=UMAP1, y=UMAP2))+
  geom_point(size=2, shape=21, aes(fill=Phenotype))+
  theme(legend.position = "top")+
  scale_fill_manual(values=c(pheno.col, pheno.col2))
dev.off()

subset.meta[subset.meta$UMAP1 < 1 & subset.meta$UMAP2 < 0, "Prediction"] <- "AD"
subset.meta[subset.meta$UMAP1 > 1 & subset.meta$UMAP2 > -4, "Prediction"] <- "PsO"
subset.meta[subset.meta$UMAP1 < 1 & subset.meta$UMAP2 > 4, "Prediction"] <- "DHR"

violin <- subset.meta
violin <- violin %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig4erythroUMAPpred.pdf", width = 6, height = 10)
ggplot(violin, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "top")+
  facet_wrap(~Prediction, ncol = 1)+
  scale_fill_manual(values=sign.col)
dev.off()

# Undetermined rashes

subset.mat <- scaled.mat[biopsy.meta$advspso, der.sign$genes]
subset.meta <- as.data.frame(scaled.thres[biopsy.meta$advspso,])
subset.meta$Phenotype <- "undetermined rashes"

pheno.advspso <- c("AD", "PsO", "DHR", "LP", "CLE", "BP")
sent.advspso <- biopsy.meta[biopsy.meta$Sentinel & biopsy.meta$Phenotype %in% pheno.advspso,
                            "Biopsy"]

merge.mat <- scaled.mat[c(sent.advspso, rownames(subset.mat)), ]
merge.meta <- biopsy.meta[rownames(merge.mat),]
merge.meta$Phenotype <- as.character(merge.meta$Phenotype)
merge.meta[rownames(subset.mat), "Phenotype"] <- "undetermined rashes"

sent.advspso.umap <- umap(merge.mat[sent.advspso,])
sent.advspso.meta <- merge.meta[sent.advspso,]
sent.advspso.meta$UMAP1 <- sent.advspso.umap$layout[,1]
sent.advspso.meta$UMAP2 <- sent.advspso.umap$layout[,2]

subset.umap <- predict(sent.advspso.umap, subset.mat)
subset.meta$UMAP1 <- subset.umap[,1]
subset.meta$UMAP2 <- subset.umap[,2]
plot.umap <- rbind(subset.meta[, c("Phenotype", "UMAP1", "UMAP2")], 
                   sent.advspso.meta[, c("Phenotype", "UMAP1", "UMAP2")])
plot.umap$Sentinels <- rownames(plot.umap) %in% rownames(sent.advspso.meta)
pdf("output/fig4rashesUMAP.pdf", width = 5, height = 3.5)
plot.umap %>% ggplot(aes(x=UMAP1, y=UMAP2))+
  geom_point(size=2, shape=21, aes(fill=Phenotype))+
  theme(legend.position = "top")+
  scale_fill_manual(values=pheno.col.all)
dev.off()

subset.meta[subset.meta$UMAP1 > 0 & subset.meta$UMAP2 > 0, "Prediction"] <- "PsO"
subset.meta[subset.meta$UMAP1 < 0 & subset.meta$UMAP1 > -6, "Prediction"] <- "AD"
subset.meta[subset.meta$UMAP1 < 0 & subset.meta$UMAP2 < 0, "Prediction"] <- "LP"

pdf("output/fig4rashesUMAPpred.pdf", width = 6, height = 10)
ggplot(violin, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "top")+
  facet_wrap(~Prediction, ncol = 1)+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG6 - DUPI non-responder

subset.mat <- scaled.mat[dupi.evol$Biopsy, der.sign$genes]
subset.meta <- as.data.frame(scaled.thres[dupi.evol$Biopsy,])
subset.meta$Treatment <- factor(dupi.evol$Treatment, 
                                   levels=c("Pre-treatment", "Post-treatment"))
subset.meta$Patient <- dupi.evol$Patient
selected.pathway <- c("Th1", "Th2", "Th17")

violin <- subset.meta
violin <- violin %>% gather(Signature, Expression, selected.pathway, factor_key = TRUE)
violin$Expression <- round(100*violin$Expression)
pdf("output/fig6dupi.pdf", width = 6, height = 4)
ggplot(data=violin, aes(x=Patient, y=Expression, fill=Signature, label=Expression)) +
  geom_bar(width = 1, colour = "black", stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = sign.col) +
dev.off()

#### SUP2 - PCA/UMAP

all.umap <- umap(biopsy.norm.mat[rownames(sentinels.mat),])
sentinels.meta$UMAP1 <- all.umap$layout[,1]
sentinels.meta$UMAP2 <- all.umap$layout[,2]
pdf("output/supp2UMAPall.pdf", width = 5, height = 3.5)
sentinels.meta %>% ggplot(aes(x=UMAP1, y=UMAP2, color=Phenotype))+
  geom_point(size=2, shape=21, col="black", aes(fill=Phenotype))+
  scale_fill_manual(values=pheno.col)+
  labs(x = "UMAP1",
       y = "UMAP2")+
  geom_text(data=subset(sentinels.meta, Phenotype =="NeuD"), aes(label=Biopsy))
dev.off()

#### SUP5 - reproducibility

repro <- read.csv(paste0(root_path, "repro2.csv"), row.names = "X")
repro.meta <- biopsy.meta[rownames(repro)[c(1,2,5,6)], ] 
# repro.meta <- biopsy.meta[rownames(repro)[c(3,4,7,8)], ]
repro.mat <- biopsy.norm.mat[rownames(repro.meta), der.sign$genes]
repro.umap <- predict(sentinels.umap, repro.mat)
repro.meta$UMAP1 <- repro.umap[,1]
repro.meta$UMAP2 <- repro.umap[,2]
plot.umap <- rbind(repro.meta[, c("Phenotype", "UMAP1", "UMAP2")], 
                   sentinels.meta[, c("Phenotype", "UMAP1", "UMAP2")])
plot.umap[rownames(plot.umap) %in% rownames(sentinels.meta), "Biospy"] <- "Sentinel"
plot.umap[is.na(plot.umap$Biospy), "Biospy"] <- "Repro"
plot.umap %>% ggplot(aes(x=UMAP1, y=UMAP2))+
  geom_point(size=2, shape=21, aes(fill=Phenotype, col=Biospy))+
  scale_fill_manual(values=pheno.col)+
  scale_color_manual(values=c("black", NA))+
  theme(legend.position = "top")+
  labs(x = "UMAP1", y = "UMAP2")

pheno.repro <- c("PsO", "AD", "LP")
sent.repro <- biopsy.meta[biopsy.meta$Phenotype %in% pheno.repro, "Biopsy"]
# sent.repro <- sent.repro[-c(58, 99, 109, 113)]

repro.mat <- scaled.mat[sent.repro, ]
repro.meta <- biopsy.meta[sent.repro,]
repro.meta$Phenotype <- as.factor(repro.meta$Phenotype)
repro.meta[, signature.names] <- scaled.thres[sent.repro, signature.names]

# repro.meta[is.na(repro.meta$Phenotype), "Phenotype"] <- "PsO"

repro.umap <- umap(repro.mat)
repro.meta$UMAP1 <- repro.umap$layout[,1]
repro.meta$UMAP2 <- repro.umap$layout[,2]
repro.meta %>% ggplot(aes(x=UMAP1, y=UMAP2, color=Phenotype))+
  geom_point(size=2, shape=21, col="black", aes(fill=Phenotype))+
  scale_fill_manual(values=pheno.col.all)+
  labs(x = "UMAP1",
       y = "UMAP2")

# repro.meta[grepl("abdo|bras|body|dos|cou|clavière|uisse|épaule|fesse|flanc|genou|jambe|joue|lombaire|main|embre|pied|poignet|ronc", repro.meta$Location), "Location2"] <- "Body"
# repro.meta[grepl("bucc|jugal|langue", repro.meta$Location), "Location2"] <- "Oral"
# repro.meta[grepl("paume|plantaire|Palm", repro.meta$Location), "Location2"] <- "Palmoplantar"
# repro.meta[grepl("poplitée|axillaire|pli", repro.meta$Location), "Location2"] <- "Intertriginous"
# repro.meta[repro.meta$Phenotype %in% "LP" & repro.meta$Location2 %in% "Intertriginous", "Location2"] <- "Body"

repro.meta[grepl("Extremities|Head|Trunk", repro.meta$Location), "Location"] <- "Body"
repro.meta["65_AD_2022.Bio.070", "Location"] <- "Intertriginous"
repro.meta["68_LP_2022.Bio.087", "Location"] <- "Body"
ggplot(subset(repro.meta, !is.na(Location)), aes(x=Location, y=Th2)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Phenotype, ncol = 3, scales = "free")+
  scale_fill_manual(values=sign.col)+
  ylim(0,1)

#### SUP6 - DEGs BP & DHR

pheno2 <- c("BP", "DHR")
sent2 <- biopsy.meta[biopsy.meta$Sentinel & biopsy.meta$Phenotype %in% pheno2,
                     "Biopsy"]

merge.mat <- scaled.mat[c(sent1, sent2), ]
merge.meta <- biopsy.meta[c(sent1, sent2), ]

merge.meta$Phenotype <- as.factor(merge.meta$Phenotype)
mm <- model.matrix(~0 + merge.meta$Phenotype)
d <- t(biopsy.mat[rownames(merge.mat),])
d <- DGEList(counts = d, group = merge.meta$Phenotype)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
levels <- colnames(coef(fit))
contr <- diag(length(levels))
contr[contr == 0] <- -1/(length(levels)-1)
rownames(contr) <- levels
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

pheno.ind <- 2
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
pb.dge <- top.table %>% filter(logFC > 1)
pb.sign <- der.sign[der.sign$genes %in% rownames(pb.dge), ]
pb.dge[pb.sign$genes, "Signature"] <- pb.sign$signature

pdf("output/sup6DEGsBP.pdf", width = 7, height = 7)
customVolcano(top.table, levels(merge.meta$Phenotype)[pheno.ind], rownames(pb.dge))
dev.off()

pheno.ind <- 4
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
toxmp.dge <- top.table %>% filter(logFC > 1)
toxmp.sign <- der.sign[der.sign$genes %in% rownames(toxmp.dge), ]
toxmp.dge[toxmp.sign$genes, "Signature"] <- toxmp.sign$signature

pdf("output/sup6DEGsDHR.pdf", width = 7, height = 7)
customVolcano(top.table, levels(merge.meta$Phenotype)[pheno.ind], rownames(toxmp.dge))
dev.off()

#### SUP7 - Barplot example mismatch

patient <- "NR_004"
patient <- "NR_007"

data.frame(Signature=signature.names, Score=scaled.thres[patient, signature.names]) %>% 
  ggplot(aes(x=Signature, y=Score, fill=Signature))+
  geom_col(width=0.9, color="black")+
  scale_fill_manual(values=sign.col)+
  scale_x_discrete(limits = signature.names)+
  ylim(0,1)+
  theme_classic()+
  ggtitle(patient)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

#### SUP8 - Check distribution of singature genes raw and normalized

pdf("output/suppfig8genexp_raw.pdf", width = 25, height = 5)
boxplot(biopsy.mat[sent1, der.sign$genes], las=2)
dev.off()
pdf("output/suppfig8genexp_norm.pdf", width = 25, height = 5)
boxplot(sentinels.mat, las=2)
dev.off()
pdf("output/suppfig8genexp_scaled.pdf", width = 25, height = 5)
boxplot(scaled.mat[sent1, der.sign$genes], las=2)
dev.off()

#### SupTable5 - Basics

biopsy.meta$cat <- biopsy.meta$Phenotype
biopsy.meta[biopsy.meta$erythro, "cat"] <- "erythro"
biopsy.meta[biopsy.meta$advspso, "cat"] <- "rashes"
biopsy.meta[biopsy.meta$post.treatment, "cat"] <- "NR"
table(biopsy.meta[, "cat"])

median(biopsy.meta$Age, na.rm = T)
min(biopsy.meta$Age, na.rm = T)
max(biopsy.meta$Age, na.rm = T)
table(biopsy.meta$Gender)/length(biopsy.meta$Gender)*100
table(biopsy.meta[, "Phenotype"])
table(biopsy.meta[biopsy.meta$advspso, "Phenotype"])
sum(is.na(biopsy.meta[biopsy.meta$advspso, "Phenotype"]))
table(biopsy.meta[biopsy.meta$erythro, "Phenotype"])
table(biopsy.meta[biopsy.meta$pre.treatment, "Phenotype"])
sum(is.na(biopsy.meta[biopsy.meta$pre.treatment, "Phenotype"]))
table(biopsy.meta[biopsy.meta$post.treatment, "Phenotype"])
sum(is.na(biopsy.meta[biopsy.meta$post.treatment, "Phenotype"]))

table(biopsy.meta$Location)/sum(!is.na((biopsy.meta$Location)))*100

#### NSF forest

library(reticulate)
use_python("c:/ProgramData/Anaconda3/python.exe")
nsf <- import("nsforest")
sc <- import("scanpy")
sent.sc <- sc$AnnData(biopsy.norm.mat[rownames(merge.mat),])
sent.sc$obs$Phenotype <- merge.meta$Phenotype
nsf.genes <- nsf$NSForest(sent.sc, "Phenotype", output_folder = "M:/DER/LABOS/IMMUNODER/Antoine/Script/nsforestReticulate/")
