#install essential packages
BiocManager::install(c("GEOquery", "clusterProfiler", "org.Hs.eg.db", "pheatmap", "EnhancedVolcano", "pathview", "biomaRt"))
BiocManager::install(c("affy", "oligo", "Biobase", "limma"), force = TRUE)
library(clusterProfiler)  #pathway enrichment
library(org.Hs.eg.db)     #human gene annotation
library(EnhancedVolcano)  #volcano plots
library(pathview)         #pathway visualization
library(affy)      #For older Affymetrix platforms
library(limma)     # For normalization and transformation
library(Biobase)   # For handling expression data
library(hgu133plus2.db)

setwd("/Users/lilianapamasa/Downloads/GSE28146/GSE28146_RAW") #set the working directory to the raw microarray counts

raw_data <- ReadAffy()

#normalize using RMA (Robust Multi-array Average)
norm_data <- affy::rma(raw_data)

#convert to expression matrix
exprs_matrix <- exprs(norm_data)


#load sample metadata
sample_info <- data.frame(
    row.names = colnames(exprs_matrix),
    condition = c(rep("Control", 8), rep("Alzheimers", 22))
)

#create the design matrix
design <- model.matrix(~ sample_info$condition)

#fit the linear model
fit <- lmFit(exprs_matrix, design)
fit <- eBayes(fit)

#get results
topTable(fit, coef=2, adjust="fdr", number=Inf) -> res

#write.csv(res, "LimmaResults.csv")
