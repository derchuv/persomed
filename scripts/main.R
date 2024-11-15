library(tidyverse)
library(limma)
library(edgeR)
library(pheatmap)
library(umap)
library(dendextend)
library(ggridges)
library(RColorBrewer)
library(nanostringr)

theme_set(theme_classic())

# Get the directory of the current script
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the root of your project
setwd(file.path(script_dir, ".."))

source("scripts/utils.R")

## Load raw data
rcc_path <- "data/rcc/."
biopsy.mat <- read_rcc(rcc_path)
biopsy.mat <- as.data.frame(biopsy.mat$raw)
rownames(biopsy.mat) <- biopsy.mat$Name
biopsy.mat <- as.data.frame(t(biopsy.mat[,c(-1:-3)]))

biopsy.meta <- read.csv("data/biopsy.csv", row.names = "X")
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
pheno2 <- c("BP", "DHR", "Erythroderma", "undetermined rashes")
pheno.col2 <-  c("#66FFFF", "#66FF66", "#FF6666", "#FFFC33")
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
biopsy.scaled.mat <- rbind(t((t(sentinels.mat) - m) / sd), 
                        t((t(test.mat) - m) / sd))
biopsy.scaled.mat <- biopsy.scaled.mat[order(row.names(biopsy.scaled.mat)),]
biopsy.scaled.mat <- biopsy.scaled.mat[rownames(biopsy.meta),]
for (signature in levels(factor(der.sign$signature))){
  genes <- der.sign[der.sign$signature %in% signature, ]$genes
  if (signature %in% colnames(biopsy.meta)){
    biopsy.meta[, signature] <- apply(biopsy.scaled.mat[, genes],1 , mean)
  } else{
    biopsy.meta[, ncol(biopsy.meta) + 1] <- apply(biopsy.scaled.mat[, genes],1 , mean)
    colnames(biopsy.meta)[ncol(biopsy.meta)] <- signature
  }
}

thres <- c(0.56,0.33,0.44,0.53,0.5,1.05,0.86)
names(thres) <- signature.names
biopsy.thres.mat <- t(t(biopsy.meta[, signature.names]) - thres)
biopsy.thres.mat <- 1/(1+exp(3*(-biopsy.thres.mat)))

sentinels.meta[, signature.names] <- biopsy.thres.mat[rownames(sentinels.meta), signature.names]
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

toplot.meta <- biopsy.meta[sent1,]
toplot.meta <- toplot.meta %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
toplot.meta$thres <- thres[toplot.meta$Signature]
pdf("output/fig2Density.pdf", width = 6, height = 18)
ggplot(toplot.meta, aes(x=Expression,y=Phenotype, fill=Signature))+
  geom_density_ridges(jittered_points = TRUE) +
  geom_vline(aes(xintercept = thres),data=toplot.meta) +
  facet_wrap(~Signature, nrow = 7) +
  theme(legend.position = "top")+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG2B - Boxplot logit transform sentinels/signature

toplot.meta <- as.data.frame(biopsy.thres.mat[sent1,])
toplot.meta$Phenotype <- biopsy.meta[sent1, "Phenotype"]
toplot.meta <- toplot.meta %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig2Boxplot.pdf", width = 6, height = 18)
ggplot(toplot.meta, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2)) +
  facet_wrap(~Phenotype, nrow = 7)+
  theme(legend.position = "top")+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG2C - Boxplot logit transform challengers/signature

subset.meta <- as.data.frame(biopsy.thres.mat[biopsy.meta$Challenger,])
subset.meta$Phenotype <- biopsy.meta[biopsy.meta$Challenger, "Phenotype"]
subset.meta$dominant.module <- signature.names[max.col(subset.meta[,signature.names])]
toplot.meta <- subset.meta %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig2Boxplot2.pdf", width = 6, height = 18)
ggplot(toplot.meta, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2)) +
  facet_wrap(~Phenotype, nrow = 7)+
  theme(legend.position = "top")+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG2D - UMAP

subset.mat <- biopsy.norm.mat[biopsy.meta$Challenger, der.sign$genes]
subset.umap <- as.data.frame(predict(sentinels.umap, subset.mat))
subset.meta$UMAP1 <- subset.umap[rownames(subset.meta),1]
subset.meta$UMAP2 <- subset.umap[rownames(subset.meta),2]
plot.umap <- rbind(subset.meta[, c("Phenotype", "UMAP1", "UMAP2")], 
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

subset.meta$cluster <- subset.meta$Phenotype
subset.meta[c("CLE_003", "AD_047"), "cluster"] <- c("LP", "Healthy") #"SW_002","PsO"
subset.meta$cluster <- factor(subset.meta$cluster, levels=levels(sentinels.meta$Phenotype))
subset.meta$Phenotype <- factor(subset.meta$Phenotype, levels=levels(sentinels.meta$Phenotype))
perf.chal <- mltest::ml_test(subset.meta$cluster, subset.meta$Phenotype)
as.data.frame(perf.chal[c("precision", "recall", "specificity")])
FM_index(subset.meta$cluster, subset.meta$Phenotype, assume_sorted_vectors = TRUE)

#### SuppTable1 - Significant pathway sentinels and chalengers

sent.thres <- as.data.frame(biopsy.thres.mat[sent1,])
sent.thres$Phenotype <- biopsy.meta[sent1, "Phenotype"]
sent.thres <- sent.thres %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)

# Get unique Phenotypes and Signatures
# phenotypes <- unique(sent.thres$Phenotype)
signatures <- unique(sent.thres$Signature)
pairs.to.plot <- list(c("Th1", "LP"),
                      c("Th2", "AD"),
                      c("Th17", "PsO"),
                      c("Neutro", "NeuD"),
                      c("Eosino", "Wells"),
                      c("IFN", "CLE"))
sent.pval <- wilcoxontestBarplot(pairs.to.plot, signatures, sent.thres)

chal.thres <- as.data.frame(biopsy.thres.mat[biopsy.meta$Challenger,])
chal.thres$Phenotype <- biopsy.meta[rownames(chal.thres), "Phenotype"]
chal.thres <- chal.thres %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)

# Get unique Phenotypes and Signatures
# phenotypes <- unique(chal.thres$Phenotype)
signatures <- unique(chal.thres$Signature)
pairs.to.plot <- list(c("Th1", "LP"),
                      c("Th2", "AD"),
                      c("Th17", "PsO"),
                      c("Neutro", "NeuD"),
                      c("IFN", "CLE"))
chal.pval <- wilcoxontestBarplot(pairs.to.plot, signatures, chal.thres)

#### FIG3 - Boxplot logit transform BP and DHR

pheno.other <- c("BP", "DHR")
sent.other <- biopsy.meta[biopsy.meta$Sentinel &
                       biopsy.meta$Phenotype %in% pheno.other, "Biopsy"]
subset.mat <- biopsy.scaled.mat[c(sent1, sent.other), der.sign$genes]
subset.meta <- as.data.frame(biopsy.thres.mat[rownames(subset.mat),])
subset.meta$Phenotype <- biopsy.meta[rownames(subset.mat), "Phenotype"]
subset.meta$dominant.module <- signature.names[max.col(subset.meta[,signature.names])]

subset.meta$threshold <- 0.7 * apply(subset.meta, 1, function(row){
  column_name <- row["dominant.module"]
  return(as.numeric(row[column_name]))
  })

subset.meta$codominant <- apply(subset.meta, 1, function(row) {
  threshold <- max(as.numeric(row["threshold"]), 0.3)
  return(paste(names(which(row[signature.names] > threshold)), collapse = ","))
})

toplot.meta <- subset.meta
toplot.meta <- toplot.meta %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig3others.pdf", width = 6, height = 10)
ggplot(toplot.meta, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot(outlier.shape = NA)+
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

merge.mat <- biopsy.scaled.mat[c(sent1, sent.other), ]
merge.meta <- biopsy.meta[c(sent1, sent.other), ]
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

merge.mat <- biopsy.scaled.mat[c(sent1, sent.other), ]
merge.meta <- biopsy.meta[c(sent1, sent.other), ]
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

subset.mat <- biopsy.scaled.mat[biopsy.meta$erythro, der.sign$genes]
subset.meta <- as.data.frame(biopsy.thres.mat[biopsy.meta$erythro,])
subset.meta$Phenotype <- "Erythroderma"
subset.meta$dominant.module <- signature.names[max.col(subset.meta[,signature.names])]
pheno.erythro <- c("AD", "PsO", "DHR")
sent.erythro <- biopsy.meta[biopsy.meta$Sentinel & 
                              biopsy.meta$Phenotype %in% pheno.erythro,
                            "Biopsy"]

merge.mat <- biopsy.scaled.mat[c(sent.erythro, rownames(subset.mat)), ]
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

# Manual settings of the limits for each category
subset.meta[subset.meta$UMAP1 > 0 & subset.meta$UMAP2 < 1, "Prediction"] <- "AD"
subset.meta[subset.meta$UMAP1 < 0 & subset.meta$UMAP2 < 0, "Prediction"] <- "PsO"
subset.meta[subset.meta$UMAP1 < 3 & subset.meta$UMAP2 > 4, "Prediction"] <- "DHR"

toplot.meta <- subset.meta
toplot.meta <- toplot.meta %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig4erythroUMAPpred.pdf", width = 6, height = 10)
ggplot(toplot.meta, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "top")+
  facet_wrap(~Prediction, ncol = 1)+
  scale_fill_manual(values=sign.col)
dev.off()

# Undetermined rashes

subset.mat <- biopsy.scaled.mat[biopsy.meta$rashes, der.sign$genes]
subset.meta <- as.data.frame(biopsy.thres.mat[biopsy.meta$rashes,])
subset.meta$Phenotype <- "undetermined rashes"

pheno.rashes <- c("AD", "PsO", "DHR", "LP", "CLE", "BP")
sent.rashes <- biopsy.meta[biopsy.meta$Sentinel & biopsy.meta$Phenotype %in% pheno.rashes,
                            "Biopsy"]

merge.mat <- biopsy.scaled.mat[c(sent.rashes, rownames(subset.mat)), ]
merge.meta <- biopsy.meta[rownames(merge.mat),]
merge.meta$Phenotype <- as.character(merge.meta$Phenotype)
merge.meta[rownames(subset.mat), "Phenotype"] <- "undetermined rashes"

sent.rashes.umap <- umap(merge.mat[sent.rashes,])
sent.rashes.meta <- merge.meta[sent.rashes,]
sent.rashes.meta$UMAP1 <- sent.rashes.umap$layout[,1]
sent.rashes.meta$UMAP2 <- sent.rashes.umap$layout[,2]

subset.umap <- predict(sent.rashes.umap, subset.mat)
subset.meta$UMAP1 <- subset.umap[,1]
subset.meta$UMAP2 <- subset.umap[,2]
plot.umap <- rbind(subset.meta[, c("Phenotype", "UMAP1", "UMAP2")], 
                   sent.rashes.meta[, c("Phenotype", "UMAP1", "UMAP2")])
plot.umap$Sentinels <- rownames(plot.umap) %in% rownames(sent.rashes.meta)
pdf("output/fig4rashesUMAP.pdf", width = 5, height = 3.5)
plot.umap %>% ggplot(aes(x=UMAP1, y=UMAP2))+
  geom_point(size=2, shape=21, aes(fill=Phenotype))+
  theme(legend.position = "top")+
  scale_fill_manual(values=pheno.col.all)
dev.off()

## Manual settings of the limits for each category
subset.meta[subset.meta$UMAP1 > 0 & subset.meta$UMAP2 > 0, "Prediction"] <- "PsO"
subset.meta[subset.meta$UMAP1 < 0 & subset.meta$UMAP2 > 3, "Prediction"] <- "AD"
subset.meta[subset.meta$UMAP1 < 0 & subset.meta$UMAP2 < 0, "Prediction"] <- "LP"

toplot.meta <- subset.meta
toplot.meta <- toplot.meta %>% gather(Signature, Expression, !!signature.names, factor_key = TRUE)
pdf("output/fig4rashesUMAPpred.pdf", width = 6, height = 10)
ggplot(toplot.meta, aes(x=Signature, y=Expression, fill=Signature)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "top")+
  facet_wrap(~Prediction, ncol = 1)+
  scale_fill_manual(values=sign.col)
dev.off()

#### FIG6 - DUPI non-responder

subset.mat <- biopsy.scaled.mat[dupi.evol$Biopsy, der.sign$genes]
subset.meta <- as.data.frame(biopsy.thres.mat[dupi.evol$Biopsy,])
subset.meta$Treatment <- factor(dupi.evol$Treatment, 
                                   levels=c("Pre-treatment", "Post-treatment"))
subset.meta$Patient <- dupi.evol$Patient
selected.pathway <- c("Th1", "Th2", "Th17")

toplot.meta <- subset.meta
toplot.meta <- toplot.meta %>% gather(Signature, Expression, selected.pathway, factor_key = TRUE)
toplot.meta$Expression <- round(100*toplot.meta$Expression)
pdf("output/fig6dupi.pdf", width = 6, height = 4)
ggplot(data=toplot.meta, aes(x=Patient, y=Expression, fill=Signature, label=Expression)) +
  geom_bar(width = 1, colour = "black", stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = sign.col) +
  facet_wrap(.~ Treatment, nrow = 1) + theme_void()
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

repro <- c("PsO_047", "PsO_048",
           "Eryth_019", "PsO_044",
           "AD_052", "AD_053",
           "UR_002", "AD_054")
repro.meta <- biopsy.meta[repro,]

repro.mat <- biopsy.norm.mat[rownames(repro.meta), der.sign$genes]
repro.umap <- predict(sentinels.umap, repro.mat)
repro.meta$UMAP1 <- repro.umap[,1]
repro.meta$UMAP2 <- repro.umap[,2]
plot.umap <- rbind(repro.meta[, c("Phenotype", "UMAP1", "UMAP2")], 
                   sentinels.meta[, c("Phenotype", "UMAP1", "UMAP2")])
plot.umap[rownames(plot.umap) %in% rownames(sentinels.meta), "Biospy"] <- "Sentinel"
plot.umap[is.na(plot.umap$Biospy), "Biospy"] <- "Repro"
pdf("output/suppfig5UMAP.pdf", width = 5, height = 3.5)
plot.umap %>% ggplot(aes(x=UMAP1, y=UMAP2))+
  geom_point(size=2, shape=21, aes(fill=Phenotype, col=Biospy))+
  scale_fill_manual(values=pheno.col)+
  scale_color_manual(values=c("black", NA))+
  theme(legend.position = "top")+
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

pheno.repro <- c("PsO", "AD", "LP")
repro <- biopsy.meta[biopsy.meta$Phenotype %in% pheno.repro, "Biopsy"]
repro.meta <- biopsy.meta[repro,]

pdf("output/suppfig5Th2.pdf", width = 5, height = 3.5)
ggplot(subset(repro.meta, !is.na(Location)), aes(x=Location, y=Th2)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Phenotype, ncol = 3, scales = "free")+
  scale_fill_manual(values=sign.col)+
  ylim(0,1)
dev.off()

#### SUP6 - DEGs BP & DHR

merge.mat <- biopsy.scaled.mat[c(sent1, sent.other), ]
merge.meta <- biopsy.meta[c(sent1, sent.other), ]

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

pdf("output/suppfig6DEGsBP.pdf", width = 7, height = 7)
customVolcano(top.table, levels(merge.meta$Phenotype)[pheno.ind], rownames(pb.dge))
dev.off()

pheno.ind <- 4
top.table <- top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = pheno.ind)
toxmp.dge <- top.table %>% filter(logFC > 1)
toxmp.sign <- der.sign[der.sign$genes %in% rownames(toxmp.dge), ]
toxmp.dge[toxmp.sign$genes, "Signature"] <- toxmp.sign$signature

pdf("output/suppfig6DEGsDHR.pdf", width = 7, height = 7)
customVolcano(top.table, levels(merge.meta$Phenotype)[pheno.ind], rownames(toxmp.dge))
dev.off()

#### SUP7 - Barplot example mismatch

patient <- "NR_004"
patient <- "NR_007"

pdf("output/suppfig7NR_007.pdf", width = 5, height = 3.5)
data.frame(Signature=signature.names, Score=biopsy.thres.mat[patient, signature.names]) %>% 
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
boxplot(biopsy.scaled.mat[sent1, der.sign$genes], las=2)
dev.off()

#### SupTable5 - Basics

biopsy.meta$cat <- biopsy.meta$Phenotype
biopsy.meta[biopsy.meta$erythro, "cat"] <- "erythro"
biopsy.meta[biopsy.meta$rashes, "cat"] <- "rashes"
biopsy.meta[biopsy.meta$post.treatment, "cat"] <- "NR"
table(biopsy.meta[, "cat"])

median(biopsy.meta$Age, na.rm = T)
min(biopsy.meta$Age, na.rm = T)
max(biopsy.meta$Age, na.rm = T)
table(biopsy.meta$Gender)/length(biopsy.meta$Gender)*100
table(biopsy.meta[, "Phenotype"])
table(biopsy.meta[biopsy.meta$rashes, "Phenotype"])
sum(is.na(biopsy.meta[biopsy.meta$rashes, "Phenotype"]))
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
