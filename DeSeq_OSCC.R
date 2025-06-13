#Author: Tanvi Nandani
#RNA-seq differential expression and functional enrichment analysis comparing OSCC,premalignant lesions (PML), and control samples.

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Rtsne)
library(ggrepel)

###### Functions #######

# Volcano Plot
make_volcano <- function(res, title){
  res_df <- as.data.frame(res)
  res_df <- res_df %>%
    mutate(Gene = rownames(res_df),
           sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"))
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    theme_minimal() +
    ggtitle(title) +
    xlab("Log2 Fold Change") +
    ylab("-log10 Adjusted p-value") +
    theme(legend.position = "bottom")
}

# GO Enrichment
run_go_enrichment <- function(deg_df, comparison_name) {
  genes <- rownames(deg_df)
  
  message(paste0("Performing GO enrichment for: ", comparison_name))
  message(paste0("Number of input genes: ", length(genes)))
  
  # Convert ENSEMBL IDs to Entrez IDs
  entrez_ids <- bitr(genes,
                     fromType = "ENSEMBL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  entrez_ids <- na.omit(entrez_ids)
  entrez_ids <- entrez_ids[!duplicated(entrez_ids$ENTREZID), ]
  
  message(paste0("Mapped ", nrow(entrez_ids), " genes to Entrez IDs."))
  
  if (nrow(entrez_ids) < 10) {
    warning("Too few genes mapped for GO enrichment. Skipping: ", comparison_name)
    return(NULL)
  }
  
  # Run GO enrichment
  go_result <- enrichGO(gene         = entrez_ids$ENTREZID,
                        OrgDb        = org.Hs.eg.db,
                        ont          = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        readable     = TRUE)
  
  # plot
  out_path <- paste0("~/Desktop/MBI/PATH828/Project/spreadsheets/", comparison_name, "_go_enrichment.csv")
  write.csv(as.data.frame(go_result), out_path)
  print(dotplot(go_result, showCategory = 10, title = paste("GO Enrichment:", comparison_name)))
}

###### Main Analysis Pipeline ######

# Load Data
counts <- read.csv("raw_counts.csv", row.names = 1)
metadata <- read.csv("sample_metadata.csv", row.names = 1)
metadata$Disease_status <- factor(metadata$Disease_status)

# Create DESeqDataSet and run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Disease_status)
dds <- DESeq(dds)

# PCA on all samples
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Disease_status")

# Identify outliers by PCA
pcaData <- plotPCA(vsd, intgroup = "Disease_status", returnData = TRUE)
print(pcaData[order(pcaData$PC1, decreasing = TRUE), ][1:5, ])

# Remove identified outliers
outliers <- c("GSM7110935", "GSM7110936", "GSM7110937")
dds_clean <- dds[, !colnames(dds) %in% outliers]
dds_clean <- DESeq(dds_clean)
vsd_clean <- vst(dds_clean, blind = FALSE)

# PCA on cleaned data
plotPCA(vsd_clean, intgroup = "Disease_status")

# Run t-SNE
tsne_input <- t(assay(vsd_clean))
set.seed(42)
tsne_out <- Rtsne(tsne_input, dims = 2, perplexity = 20, verbose = TRUE)
tsne_df <- data.frame(tsne_out$Y)
tsne_df$Sample <- rownames(tsne_input)
tsne_df$Group <- colData(vsd_clean)$Disease_status

ggplot(tsne_df, aes(X1, X2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE of Samples", x = "tSNE-1", y = "tSNE-2")

# Differential Expression Analysis
dds_clean$Disease_status <- relevel(dds_clean$Disease_status, ref = "Control")
res_oscc_vs_control <- results(dds_clean, contrast = c("Disease_status", "OSCC", "Control"))
res_pml_vs_control <- results(dds_clean, contrast = c("Disease_status", "premalignant lesions", "Control"))
res_pml_vs_oscc <- results(dds_clean, contrast = c("Disease_status", "premalignant lesions", "OSCC"))

# Volcano Plots
make_volcano(res_oscc_vs_control, "OSCC vs Control")
make_volcano(res_pml_vs_control, "PML vs Control")
make_volcano(res_pml_vs_oscc, "PML vs OSCC")

# GO Enrichment
top_oscc_vs_control <- res_oscc_vs_control[which(res_oscc_vs_control$padj < 0.05), ]
top_pml_vs_control <- res_pml_vs_control[which(res_pml_vs_control$padj < 0.05), ]
top_pml_vs_oscc <- res_pml_vs_oscc[which(res_pml_vs_oscc$padj < 0.05), ]

run_go_enrichment(top_oscc_vs_control, "OSCC_vs_Control")
run_go_enrichment(top_pml_vs_control, "PML_vs_Control")
run_go_enrichment(top_pml_vs_oscc, "PML_vs_OSCC")



