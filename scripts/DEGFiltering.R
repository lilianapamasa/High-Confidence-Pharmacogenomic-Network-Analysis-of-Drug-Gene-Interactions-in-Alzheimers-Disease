#filter for genes with raw P.Value < 0.05
sig_genes_pval <- subset(res, P.Value < 0.05)

#extract probe IDs
probe_ids <- rownames(sig_genes_pval)  # Assuming rownames contain probe IDs

#convert probe IDs to gene symbols
gene_symbols <- mapIds(hgu133plus2.db, 
                       keys = probe_ids, 
                       column = "SYMBOL", 
                       keytype = "PROBEID", 
                       multiVals = "first")

#add gene symbols to the dataset
sig_genes_pval$GeneSymbol <- gene_symbols

#remove rows where GeneSymbol is missing
sig_genes_pval <- na.omit(sig_genes_pval)

#view the first few rows
head(sig_genes_pval)

#write.csv(sig_genes_pval, "SigGenes_LimmaResults.csv")

#extract unique gene symbols
gene_symbols_list <- unique(sig_genes_pval$GeneSymbol)

#map gene symbols to Entrez IDs
entrez_ids <- bitr(gene_symbols_list, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
#KEGG enrichment analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                          organism = 'hsa', 
                          pvalueCutoff = 0.05)

#view the enrichment results
#summary(kegg_enrich)

#plot the top KEGG pathways
kegg_bar <- barplot(kegg_enrich, showCategory = 10, title = "Top 10 Enriched KEGG Pathways")

#extract KEGG results as a data frame
kegg_results <- as.data.frame(kegg_enrich)

#filter for drug metabolism-related pathways using keywords
drug_metabolism_pathways <- kegg_results[grep('drug|alzheimer', kegg_results$Description, ignore.case = TRUE),]

#view the filtered pathways
drug_metabolism_pathways
