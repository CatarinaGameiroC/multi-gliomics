setwd("/Users/catarina/Library/Mobile Documents/com~apple~CloudDocs/Tese/code/Wed13Nov")

library("DESeq2")
library("edgeR")


mrna <- as.matrix(read.csv("DataSets/original_assays/assay_rna_coding.csv", row.names = 1))

mrna <- t(mrna)



# normalization as in https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# the difference for what we did is that we do not perform voom stabilization
keep <- rowSums(cpm(mrna) > 1) >= 5
mrna <- mrna[keep,]

mrna <- DGEList(counts = mrna)
mrna <- normLibSizes(mrna, method = "TMM")



classification <- as.data.frame(read.csv("DataSets/SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv", row.names = 1))




#################
# GBM/ASTRO/OLIGO
#################
# groups
classification$classification.2021 <- ifelse(classification $classification.2021 == "glioblastoma", "GBM", 
                                       ifelse(classification $classification.2021  == "astrocytoma", "ASTRO",
                                              ifelse(classification $classification.2021  == "oligodendroglioma", "OLIGO",
                                                     "UNCLASS")))


classification <- classification[,-c(1, 2), drop = FALSE]


missing_samples <- setdiff(colnames(mrna), rownames(classification))
new_rows <- data.frame(classification.2021 = rep("UNCLASS", length(missing_samples)), 
                       row.names = missing_samples)
classification <- rbind(classification, new_rows)


classes <- classification[colnames(mrna), "classification.2021"]


classification[c("TCGA-06-5856", "TCGA-06-0152"), ]

mrna$samples$group <- classification[rownames(mrna$samples), "classification.2021"]


# design matrix
design <- model.matrix(~ 0 + group, data = mrna$samples)
colnames(design) <- levels(factor(mrna$samples$group))
design


# dispersion
mrna <- estimateDisp(mrna, design, robust=TRUE)
mrna$common.dispersion


# fit 
fit <- glmQLFit(mrna, design, robust=TRUE)
head(fit$coefficients)

# the contrasts we want to analyse
con1 <- makeContrasts(GBM - ASTRO, levels=design)
con2 <- makeContrasts(GBM- OLIGO, levels=design)
con3 <- makeContrasts(ASTRO - OLIGO, levels=design)

# do the test
qlf1 <- glmQLFTest(fit, contrast=con1)
qlf2 <- glmQLFTest(fit, contrast=con2)
qlf3 <- glmQLFTest(fit, contrast=con3)



# explore
decideTests(qlf1)
summary(decideTests(qlf1))
# -1*ASTRO 1*GBM
# Down             6634 where is -1
# NotSig           2516 where is 0
# Up               7346 where is 1




# put in csv files the the tables with the p value, FDR, ... and the file about the decision
write.csv(qlf1$table, "DGEresults/mrna_GBMvsASTRO.csv", row.names = TRUE)
write.csv(qlf2$table, "DGEresults/mrna_GBMvsOLIGO.csv", row.names = TRUE)
write.csv(qlf3$table, "DGEresults/mrna_ASTROvsOLIGO.csv", row.names = TRUE)





#################
# GBM/LGG
#################

classification$classification.2021 <- ifelse(classification$classification.2021 == "glioblastoma", "GBM", 
                                             ifelse(classification$classification.2021  == "astrocytoma", "LGG",
                                                    ifelse(classification $classification.2021  == "oligodendroglioma", "LGG",
                                                           "UNCLASS")))


classification <- classification[,-c(1, 2), drop = FALSE]


missing_samples <- setdiff(colnames(mrna), rownames(classification))
new_rows <- data.frame(classification.2021 = rep("UNCLASS", length(missing_samples)), 
                       row.names = missing_samples)
classification <- rbind(classification, new_rows)


classes <- classification[colnames(mrna), "classification.2021"]


mrna$samples$group <- classification[rownames(mrna$samples), "classification.2021"]


# design matrix
design <- model.matrix(~ 0 + group, data = mrna$samples)
colnames(design) <- levels(factor(mrna$samples$group))
design


# dispersion
mrna <- estimateDisp(mrna, design, robust = TRUE)
mrna$common.dispersion


# fit 
fit <- glmQLFit(mrna, design, robust=TRUE)
head(fit$coefficients)

# the contrasts we want to analyse
con1 <- makeContrasts(GBM - LGG, levels=design)

# do the test
qlf1_type <- glmQLFTest(fit, contrast=con1)




# put in csv files the the tables with the p value, FDR, ... and the file about the decision
write.csv(qlf1_type$table, "DGEresults/mrna_GBMvsLGG.csv", row.names = TRUE)







## VOLCANO PLOTS 
qlf1 <- as.data.frame(read_csv("DGEresults/mrna_GBMvsASTRO.csv"))
rownames(qlf1) <- qlf1[,1]
qlf1 <- qlf1[-1]


qlf2 <- as.data.frame(read_csv("DGEresults/mrna_GBMvsOLIGO.csv"))
rownames(qlf2) <- qlf2[,1]
qlf2 <- qlf2[-1]


qlf3 <- as.data.frame(read_csv("DGEresults/mrna_ASTROvsOLIGO.csv"))
rownames(qlf3) <- qlf3[,1]
qlf3 <- qlf3[-1]


qlf1_type <- as.data.frame(read_csv("DGEresults/mrna_GBMvsLGG.csv"))
rownames(qlf1_type) <- qlf1_type[,1]
qlf1_type <- qlf1_type[-1]



library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")


volcano_plot <- function(df, title_string){
  # according to the cutoffs said
  df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > 1 & df$PValue < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -1 & df$PValue < 0.05] <- "DOWN"
  head(df[order(df$PValue) & df$diffexpressed == 'DOWN', ])
  
  df$delabel <- ifelse(rownames(df) %in% rownames(head(df[order(df$PValue), ], 30)), rownames(df), NA)
  df$delabel <- ifelse(!is.na(df$delabel), 
                         bitr(df$delabel, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL, 
                         NA)
  
  ggplot(data = df, aes(x = logFC, y = -log10(PValue), col = diffexpressed,label = delabel)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 1) +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + 
    labs(color = "",x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    scale_x_continuous(breaks = seq(-10, 10, 2)) + 
    ggtitle(title_string) + 
    geom_text_repel(max.overlaps = Inf)
}


volcano_plot(qlf1, "DGE mRNA - GBMvsASTRO")
volcano_plot(qlf2, "DGE mRNA - GBMvsOLIGO")
volcano_plot(qlf3, "DGE mRNA - ASTROvsOLIGO")


volcano_plot(qlf1_type, "DGE mRNA - GBMvsLGG")

