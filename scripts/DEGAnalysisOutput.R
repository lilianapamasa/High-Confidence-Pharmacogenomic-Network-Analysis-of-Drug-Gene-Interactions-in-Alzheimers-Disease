#set thresholds
logFC_threshold <- 1
pval_threshold <- 0.05

#create significance column
res$Significance <- "Not Significant"
res$Significance[res$logFC > logFC_threshold & res$P.Value < pval_threshold] <- "Upregulated"
res$Significance[res$logFC < -logFC_threshold & res$P.Value < pval_threshold] <- "Downregulated"

#create the volcano plot
library(ggplot2)

volcano <- ggplot(res, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot (Raw P-values)",
       x = "Log2 Fold Change", y = "-Log10(P-value)")

#save as png
ggsave("volcano_plot.png", plot = volcano, width = 10, height = 8, dpi = 300, bg = "white")

#load annotation package
library(hgu133plus2.db)

#subset the top 10 genes by raw P.Value
top_genes <- head(res[order(res$P.Value), ], 10)

probe_ids <- rownames(top_genes) #set rownames to prode ids

#map to gene symbols
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

#add gene symbols to the results
top_genes$GeneSymbol <- gene_symbols

#remove rows where GeneSymbol is missing
top_genes <- na.omit(top_genes)

#show top genes
top_genes[, c("GeneSymbol", "logFC", "P.Value")]

#load ggplot2
library(ggplot2)

#reorder GeneSymbol factor based on logFC
top_genes$GeneSymbol <- factor(top_genes$GeneSymbol, levels = top_genes$GeneSymbol[order(top_genes$logFC)])

#bar plot
bar <- ggplot(top_genes, aes(x = GeneSymbol, y = logFC, fill = logFC > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "steelblue3", "FALSE" = "tomato"), labels = c("Downregulated", "Upregulated")) +
  coord_flip() +
  labs(title = "Top 10 Differentially Expressed Genes",
       x = "Gene Symbol",
       y = "log2 Fold Change",
       fill = "Direction") +
  theme_minimal()

#save as png
ggsave("bar_plot.png", plot = bar, width = 10, height = 8, dpi = 300, bg = "white")
