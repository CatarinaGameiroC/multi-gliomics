library("ChAMP")

setwd("/Users/catarina/Library/Mobile Documents/com~apple~CloudDocs/Tese/code/Wed13Nov")

methy <- read.csv("DataSets/original_assays/assay_methylation_for_dge.csv", 
               row.names = 1)


methy_beta <- t(M2B(methy))




classification <- as.data.frame(read.csv("DataSets/SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv", 
                                         row.names = 1))

classification$classification.2021 <- ifelse(classification $classification.2021 == "glioblastoma", "GBM", 
                                             ifelse(classification $classification.2021  == "astrocytoma", "ASTRO",
                                                    ifelse(classification $classification.2021  == "oligodendroglioma", "OLIGO",
                                                           "UNCLASS")))

classification <- classification[,-c(1, 2), drop = FALSE]

missing_samples <- setdiff(colnames(methy_beta), rownames(classification))
new_rows <- data.frame(classification.2021 = rep("UNCLASS", length(missing_samples)), 
                       row.names = missing_samples)
classification <- rbind(classification, new_rows)


classes <- classification[colnames(methy_beta), "classification.2021"]

class <- ifelse(classes == "GBM", "GBM", 
               ifelse(classes  == "ASTRO", "LGG",
                      ifelse(classes  == "OLIGO", "LGG",
                             "UNCLASS")))



gbm_vs_lgg_methy <- champ.DMP(beta = methy_beta,
                   pheno = class,
                   compare.group = c("GBM", "LGG"),
                   adjPVal = 0.05 ,
                   adjust.method = "BH",
                   arraytype = "450K")


astro_vs_oligo_methy <- champ.DMP(beta = methy_beta,
                              pheno = classes,
                              compare.group = c("ASTRO", "OLIGO"),
                              adjPVal = 0.05 ,
                              adjust.method = "BH",
                              arraytype = "450K")

oligo_vs_gbm_methy <- champ.DMP(beta = methy_beta,
                                  pheno = classes,
                                  compare.group = c("GBM", "OLIGO"),
                                  adjPVal = 0.05 ,
                                  adjust.method = "BH",
                                  arraytype = "450K")


library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")


volcano_plot <- function(df){
 
  
  df$diffexpressed <- "NO"
  df$diffexpressed[!is.na(df$logFC)  & df$logFC > 0.3 & df$P.Value < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -0.3 & df$P.Value < 0.05] <- "DOWN"
  
  
  df_out <- df
  df <- df[df$P.Value > 0, ]
  
  df$delabel <- ifelse(rownames(df) %in% rownames(head(df[order(df$P.Value ), ], 30)), rownames(df), NA)
  
  
  p <- ggplot(data = df, aes(x = deltaBeta, y = -log10(P.Value), col = diffexpressed,label = delabel)) +
    geom_vline(xintercept = c(-0.3, 0.3), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 1) +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    coord_cartesian(ylim = c(0, 150), xlim = c(-1, 1)) + 
    labs(color = "",x = expression("DeltaBeta"), y = expression("-log"[10]*"p-value")) + 
    geom_text_repel(max.overlaps = Inf) +
    theme(
      axis.text.y = element_text(size=rel(text_size), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(text_size), vjust=0.5, color='black'),
      axis.title.y= element_text(size=rel(text_size), hjust=1, color='black'),
      axis.title.x=element_text(size=rel(text_size), vjust=1, color='black'),
      legend.text = element_text(size = rel(text_size), color = 'black'),
      legend.title = element_text(size = rel(text_size), color = 'black'),
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      panel.background = element_blank(),
    ) 
  
  
  
  return (list(p, df_out))
}

# promoter: c("TSS1500","TSS200","1stExon")

methy_class <- volcano_plot(gbm_vs_lgg_methy$GBM_to_LGG)

write.csv(methy_class[[2]], "DGEresults/methy_GBMvsLGG.csv", row.names = TRUE)
ggsave(filename = "DGEresults/methy_GBMvsLGG.pdf", plot = methy_class[[1]], width = 10, height = 7, dpi = 300)

methy_class[[2]][methy_sel_f1, "diffexpressed"]
methy_class[[2]][methy_sel_f1, c("deltaBeta", "logFC") ]
methy_class[[2]][methy_sel_f1, ]
methy_class[[1]]


# all methy_sel_f1 are up-regulated in LGG


methy_subclas <- volcano_plot(astro_vs_oligo_methy$ASTRO_to_OLIGO)
write.csv(methy_subclas[[2]], "DGEresults/methy_ASTROvsOLIGO.csv", row.names = TRUE)
ggsave(filename = "DGEresults/methy_ASTROvsOLIGO.pdf", plot = methy_subclas[[1]], width = 10, height = 7, dpi = 300)


methy_subclas[[2]][methy_sel_f3, "diffexpressed"]   # mean(pOLIGO) > mean(pASTRO)
methy_subclas[[2]][methy_sel_f3, c("deltaBeta", "logFC")] 
methy_subclas[[2]][methy_sel_f3, c("gene", "diffexpressed")] 
methy_subclas[[1]]



methy_gbm_oligo_f3 <- volcano_plot(oligo_vs_gbm_methy$GBM_to_OLIGO)

methy_gbm_oligo_f3[[1]]
methy_gbm_oligo_f3[[2]][methy_sel_f3, ]






plot_methylation_vs_expression <- function(data, cpg, gene, groups, save_path = NULL) {
  
  data <- data[!is.na(data$group) & data$group %in% groups, ]

  
  name <- info.mrna[gene,"gene_name"]
  p <- ggplot(data = data, aes_string(x = cpg, y =gene, color = "group")) +
    geom_point() +
    geom_smooth(method = "lm", fullrange = FALSE) +
    scale_color_manual(values = c('#FF7400', '#009999')) +
    labs(
      x =  cpg,
      y =  name,
      color = "Subtype"
    ) +
   # coord_cartesian(xlim = c(min(data$cpg), 1)) +
    theme(
      axis.text.y = element_text(size=rel(text_size), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(text_size), vjust=0.5, color='black'),
      axis.title.y= element_text(size=rel(text_size), hjust=1, color='black'),
      axis.title.x=element_text(size=rel(text_size), vjust=1, color='black'),
      legend.text = element_text(size = rel(text_size), color = 'black'),
      legend.title = element_text(size = rel(text_size), color = 'black'),
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      panel.background = element_blank(),
    )  +
    stat_cor(method = "spearman", aes_string(x = cpg, y = gene, color = "group"), label.x = 0.6)
  
  print(p)
}




df_test <- data.frame(t(M2B(t(omics.list$Methylation[,methy_sel_f3]))), omics.list$mRNA[,mrna_sel_f3_], 
                      group = metadata$Subtype)


cpg_matrix <- as.data.frame(t(M2B(t(omics.list$Methylation[, methy_sel_f3]))))
gene_matrix <- as.data.frame(omics.list$mRNA[, mrna_sel_f3_])


results <- data.frame(CpG = character(), Gene = character(), Correlation = numeric())


for (cpg_name in colnames(cpg_matrix)) {
  for (gene_name in colnames(gene_matrix)) {
    corr_value <- cor(cpg_matrix[[cpg_name]], gene_matrix[[gene_name]], method = "spearman", use = "complete.obs")
    results <- rbind(results, data.frame(CpG = cpg_name, Gene = gene_name, Correlation = corr_value))
  }
}



results[results$Correlation < -0.5,]
p1 <- plot_methylation_vs_expression(data = df_test, cpg = "cg17203063", gene = "ENSG00000162711", groups = c("OLIGO", "GBM"))
p2 <- plot_methylation_vs_expression(data = df_test, cpg = "cg17203063", gene = "ENSG00000162711", groups = c("OLIGO", "ASTRO"))

p3 <- plot_methylation_vs_expression(data = df_test, cpg = "cg17203063", gene = "ENSG00000110876", groups = c("OLIGO", "ASTRO"))
p4 <- plot_methylation_vs_expression(data = df_test, cpg = "cg17203063", gene = "ENSG00000182578", groups = c("OLIGO", "ASTRO"))

ggsave(filename = "DGEresults/cg17203063_ENSG00000162711_oliog_gbm.pdf", plot = p1, width = 10, height = 7, dpi = 300)
ggsave(filename = "DGEresults/cg17203063_ENSG00000162711_oliog_astro.pdf", plot = p2, width = 10, height = 7, dpi = 300)

ggsave(filename = "DGEresults/cg17203063_ENSG00000110876.pdf", plot = p3, width = 10, height = 7, dpi = 300)
ggsave(filename = "DGEresults/cg17203063_ENSG00000182578.pdf", plot = p4, width = 10, height = 7, dpi = 300)


