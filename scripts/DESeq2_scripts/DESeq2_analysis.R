#R-Script for Step 5-7 of workflow: 
# 5. Exploratory data analysis, 6. Differential expression analysis, 7. Over-representation analysis.


#-----------------------------------------------------------------------------------------------------

#Import libraries
library(DESeq2)
library(ggplot2)
library(corto)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(dplyr)
library(AnnotationDbi)
library(grid)
library(cowplot)
library(viridis)
library(reshape2)
library(enrichplot)

#-----------------------------------------------------------------------------------------------------

#Define the file path for reading the featureCount data
file_path_featureCounts_data = "/Users/mariokummer/Desktop/RNA-Sequencing/data/DESeq_2_data/all_samples_counts.txt"

# Define the file path for saving the plots
base_save_path <- "/Users/mariokummer/Desktop/RNA-Sequencing-Course/results/DESeq2/"

#-----------------------------------------------------------------------------------------------------

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5. Exploratory data analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#=============================================================================================================
 # 5.1 Assessment on how the samples cluster based on their gene expression profiles --> PCA Plot
#=============================================================================================================

# Read FeatureCounts data
counts <- read.table(file_path_featureCounts_data, header=TRUE, row.names=1)

# Formatting: Rename Column header of dataset to only contain sample name, remove first line and the columns containing (Chr, Start, End, Strand and Length)
new_names <- sub(".*(SRR[0-9]+).*", "\\1", colnames(counts))
colnames(counts) <- new_names
counts <- counts[-1, ] 
counts <- counts[, c(0, 6:ncol(counts))]

# Add sample info metadata for one and two factor
sample_info <- data.frame(
  row.names = colnames(counts),
  genotype = c("WT",
               "WT",
               "WT",
               "WT",
               "WT",
               "DKO",
               "DKO",
               "DKO",
               "DKO",
               "WT",
               "WT",
               "WT",
               "DKO",
               "DKO",
               "DKO"),
  infection = c("Case",
                "Case",
                "Case",
                "Case",
                "Case",
                "Case",
                "Case",
                "Case",
                "Case",
                "Control",
                "Control",
                "Control",
                "Control",
                "Control",
                "Control"),
  condition = c("Lung_WT_Case",
                "Lung_WT_Case",
                "Lung_WT_Case",
                "Lung_WT_Case",
                "Lung_WT_Case",
                "Lung_DKO_Case",
                "Lung_DKO_Case",
                "Lung_DKO_Case",
                "Lung_DKO_Case",
                "Lung_WT_Control",
                "Lung_WT_Control",
                "Lung_WT_Control",
                "Lung_DKO_Control",
                "Lung_DKO_Control",
                "Lung_DKO_Control"))

# Create the DESeqDataSet object resp. matrix
dds_one_factor <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

#Run DESeq: Normalization + Dispersion Estimation + Differential Expression Testing
dds_one_factor <- DESeq(dds_one_factor)

#Regularized log transformation with rlog
#Blind = TRUE means the transformation does not use the experimental design
vsd_one_factor <- rlog(dds_one_factor, blind = TRUE)

# Extract PCA plot data
PCA_plot_data <- plotPCA(vsd_one_factor, intgroup = "condition", returnData = TRUE)

# Extract percent variance
percentVar <- attr(PCA_plot_data, "percentVar")
pc1_var <- round(100 * percentVar[1],2)
pc2_var <- round(100 * percentVar[2],2)

#Plot PCA as custom ggplot and store in variable
PCA_plot <- ggplot(PCA_plot_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4.5, shape = 17) +                     
  geom_text_repel(
    aes(label = name),
    direction = "y",
    colour = "black",
    segment.color = NA,
    size = 2) + 
  xlab(paste0("PC1: (", pc1_var, "%)")) +
  ylab(paste0("PC2: (", pc2_var, "%)")) +
  ggtitle("Principal component analysis (PCA) of lung samples") +
  coord_fixed() +                                       
  theme_classic(base_size = 12) +                      
  theme(
    axis.line = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "gray90", linetype = "dotted"),
    panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 8), )  + 
    scale_color_manual(
    values = c(
      "Lung_WT_Case" = "#FF6F61",
      "Lung_DKO_Case" = "#FFA500",
      "Lung_WT_Control" = "#6B5B95",
      "Lung_DKO_Control" = "#8DD3C7"),
    labels = c(
      "Lung_WT_Case" = "Lung Wild-Type Case (T. gondii)",
      "Lung_DKO_Case" = "Lung Double-Knockout (Ifnar−/− × Ifngr−/−) Case (T. gondii)",
      "Lung_WT_Control" = "Lung Wild-Type Control",
      "Lung_DKO_Control" = "Lung Double-Knockout (Ifnar−/− × Ifngr−/−) Control")) +
    labs(color = NULL) +
    guides(color = guide_legend(nrow = 2,byrow = TRUE))

#Show PCA plot
PCA_plot

# Save PCA plot as png
ggsave(paste0(base_save_path,"PCA_plot_Lung_samples.png"), plot = PCA_plot, width = 8, height = 6, dpi = 300)    



#-----------------------------------------------------------------------------------------------------



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6. Differential Expression (DE) analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#=============================================================================================================
  # 6.1.1 Volcano plot for each experimental group contrast (One factor)
#=============================================================================================================

# Extract the DE test results for each comparison between experimental groups
res_wt <- results(dds_one_factor, contrast = c("condition", "Lung_WT_Case", "Lung_WT_Control"))
res_dko <- results(dds_one_factor, contrast = c("condition", "Lung_DKO_Case", "Lung_DKO_Control"))
res_case <- results(dds_one_factor, contrast = c("condition", "Lung_DKO_Case", "Lung_WT_Case"))
res_ctrl <- results(dds_one_factor, contrast = c("condition", "Lung_WT_Control", "Lung_DKO_Control"))

# List of DESeq2 results with descriptive names
res_list <- list(
  WT_vs_Control = res_wt,
  DKO_vs_Control = res_dko,
  WT_vs_DKO_Case = res_case,
  WT_vs_DKO_Control = res_ctrl
)

# Initialize list to store plots
volcano_plots <- list()

# Define custom title for Volcano-Plots
custom_titles <- c(
  "Lung WT Case vs Lung WT Control",
  "Lung DKO Case vs Lung DKO Control",
  "Lung DKO Case vs Lung WT Case",
  "Lung WT Control vs Lung DKO Control"
)

# Define custom subtitle for Volcano-Plots
custom_subtitles <- c(
  "Lung Wild-Type Case (T. gondii) vs Lung Wild-Type Control",
  "Lung Double-Knockout (Ifnar−/− × Ifngr−/−) Case (T. gondii) vs. Lung Double-Knockout (Ifnar−/− × Ifngr−/−) Control",
  "Lung Double-Knockout (Ifnar−/− × Ifngr−/−) Case (T. gondii) vs. Lung Wild-Type Case (T. gondii)",
  "Lung Wild-Type Control vs. Lung Double-Knockout (Ifnar−/− × Ifngr−/−) Control"
)

#padj value threshold for each contrast from res_list, needed for determination which data points with label (Gen name) are shown in plot
padj_treshold = c(100,50,100,10)

# Loop through each DESeq2 result to create individual volcano plots for each contrast between experimental groups
for (i in seq_along(res_list)) {
  
  # Get the right data set and title and subtitle
  res <- res_list[[i]]
  title_text <- custom_titles[i]
  subtitle_text <- custom_subtitles[i]
  
  # Convert res to data frame
  res <- as.data.frame(res_list[[i]])
  
  # Map ENSEMBL IDs to gene symbols: 1. Get Ensemble-ID from rownames of res 2. Convert Ensembl-IDs to gen symbols
  ens_id <- rownames(res)
  symbols <- mapIds(org.Mm.eg.db, keys = ens_id,
                    column = "SYMBOL", keytype = "ENSEMBL")
  
  # Replace NAs with ENSEMBL ID (preserve names!)
  symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
  
  # Assign symbols to res data
  res$symbol <- symbols
  
  # Calculate the -log10 for padj for better representation in plot, assign to new column in res
  res$negLog10padj <- -log10(res$padj)

  # Filter significant genes, not NA & negative Log10(padj) is greater than the corresponding threshold for which show label in plot
  sig_res_label <- res[!is.na(res$negLog10padj) & res$negLog10padj > padj_treshold[i], ]

  # Top up- and down-regulated genes log2FoldChange > |2.5|, again just for showing label names in plot
  up_genes <- sig_res_label[sig_res_label$log2FoldChange > 2.5, ]
  up_genes <- as.data.frame(up_genes)
  top_up <- up_genes[order(up_genes$negLog10padj), ]
  down_genes <- sig_res_label[sig_res_label$log2FoldChange < - 2.5, ]
  down_genes <- as.data.frame(down_genes)
  top_down <- down_genes[order(down_genes$negLog10padj), ]
  
  # All top genes (Up and down reg.), again just for showing label names in plot
  top_genes <- c(top_up$symbol, top_down$symbol)
  
  # Remove NA padj rows
  res_valid <- res[!is.na(res$padj), ]
  
  # Get count of total DE genes, DE genes below padj < 0.005, Up-/down-reg. genes
  total_genes <- nrow(res_valid)
  number_genes_below <- sum(res_valid$padj < 0.005)
  up_reg_genes <- sum(res_valid$padj < 0.005 & res_valid$log2FoldChange > 0)
  down_reg_genes <- sum(res_valid$padj < 0.005 & res_valid$log2FoldChange < 0)
  
  # Create custom volcano plot for each contrast (Lung WT Case vs. Lung DKO Case, ... etc.), for explanation of parameter/arguments used see EnhancedVolcano package
  volcano_plot <- EnhancedVolcano(
    res_valid,                                                       
    pCutoff = 0.005,                                                  
    FCcutoff = 1.0,                                            
    lab = res_valid$symbol,                                        
    x = 'log2FoldChange',                                            
    y = 'padj',                                                     
    subtitle = paste0("Differential expression:", subtitle_text),   
    labSize = 1.5,                                                   
    drawConnectors = TRUE,                                           
    widthConnectors = 0.25,                                          
    typeConnectors = "open",                                        
    endsConnectors = "first",                                        
    colConnectors = "black",                                        
    selectLab = top_genes,                                          
    legendLabels = c('Not significant: p-value > 0.005','High fold change: |Log2(FC)| ≥ 1','Significant: p-value < 0.005', 'High fold change: |Log2(FC)| ≥ 1 & Significant: p-value < 0.005'),
    legendPosition = 'bottom',                                     
    legendLabSize = 5,                                              
    legendIconSize = 4,                                        
    shape = 23,                                                      
    title = paste0("Volcano plot for contrast ", title_text),      
    pointSize = 1.25,                                              
    captionLabSize = 8,                                       
    caption = bquote( #Custom caption/note to plot, Total/Up/Down regulated genes
      "Total DE genes: " * .(total_genes) * 
        "; Significant DE genes (padj < 0.005): " * .(number_genes_below) *
        " [Up-reg.: " * .(up_reg_genes) * ", Down-reg.: " * .(down_reg_genes) * "]")) + 
    theme( 
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    plot.title = element_text(size = 16),      
    plot.subtitle = element_text(size = 8),     
    axis.title = element_text(size = 10),       
    axis.text = element_text(size = 10),        
    axis.title.x = element_text(size = 10, margin = margin(t = 1)),   
    axis.text.y = element_text(size = 10),   
    legend.text = element_text(size = 10),        
    panel.grid.major = element_line(color = "grey80", size = 0.3),  
    panel.grid.minor = element_line(color = "grey90", size = 0.2),   
    legend.spacing.y = unit(0.5, "cm"),  
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  
    legend.box.margin = margin(0, 0, 0, 0),              
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + 
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) 
  
  # Store the plot in the list
  volcano_plots[[i]] <- volcano_plot
}

# Arrange 2x2 grid
combined_plot <- (volcano_plots[[1]] | volcano_plots[[2]]) /
  (volcano_plots[[3]] | volcano_plots[[4]])


# Add separation lines in 2x2 grid --> 4 boxes with separate plot
combined_plot_with_lines <- ggdraw() +
  draw_plot(combined_plot, 0, 0, 1, 1) +
  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "grey80", size = 0.5) +
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey80", size = 0.5) +
  draw_line(x = c(0, 1), y = c(0, 0), color = "grey60", size = 1) +
  draw_line(x = c(0, 1), y = c(1, 1), color = "grey60", size = 1) +
  draw_line(x = c(0, 0), y = c(0, 1), color = "grey60", size = 1) +
  draw_line(x = c(1, 1), y = c(0, 1), color = "grey60", size = 1)

# Save as 2x2 combined Volcano plots as png
ggsave(paste0(base_save_path,"combined_volcano_plots.png"), plot = combined_plot, width = 16, height = 12, dpi = 300)



#=============================================================================================================
  # 6.1.2  Volcano plot two factor: Lung WT Case vs. Lung DKO Case 
#=============================================================================================================

# Create the DESeqDataSet object resp.  for two factor contrast Lung WT Case vs. Lung DKO Case
dds_two_factor <- DESeqDataSetFromMatrix(countData = counts,
                                         colData = sample_info,
                                         design = ~ genotype + infection + genotype:infection)

#Run DESeq: Normalization + Dispersion Estimation + Differential Expression Testing
dds_two_factor <- DESeq(dds_two_factor)

# Extract the DE test results for two factor Lung WT Case vs. Lung DKO Case
res_case_two_factor <- results(dds_two_factor, list(c("genotype_WT_vs_DKO", "genotypeWT.infectionControl")))

# Map ENSEMBL IDs to gene symbols: 1. Get Ensemble-ID from row names of res 2. Convert Ensembl-IDs to gen symbols
ens_id_two_factor <- rownames(res_case_two_factor)
symbols_two_factor <- mapIds(org.Mm.eg.db, keys = ens_id_two_factor,
              column = "SYMBOL", keytype = "ENSEMBL")

# Replace NAs with ENSEMBL ID (preserve names!)
symbols_two_factor[is.na(symbols_two_factor)] <- names(symbols_two_factor)[is.na(symbols_two_factor)]

# Assign symbol to res_two_factor data
res_case_two_factor$symbol <- symbols_two_factor

# Calculate the -log10 for padj for better representation in plot, assign to new column in res
res_case_two_factor$negLog10padj <- -log10(res_case_two_factor$padj)

# Filter significant genes, not NA & negative Log10(padj) is greater than the corresponding threshold 10, only for 
sig_res_two_factor <- res_case_two_factor[!is.na(res_case_two_factor$negLog10padj) & res_case_two_factor$negLog10padj > 10, ]

# Top up- and down-regulated genes log2FoldChange > |1|
up_genes_two_factor <- sig_res_two_factor[sig_res_two_factor$log2FoldChange > 1, ]
up_genes_two_factor <- as.data.frame(up_genes_two_factor)
top_up_two_factor <- up_genes_two_factor[order(up_genes_two_factor$negLog10padj), ]
down_genes_two_factor <- sig_res_two_factor[sig_res_two_factor$log2FoldChange < - 1, ]
down_genes_two_factor <- as.data.frame(down_genes_two_factor)
top_down_two_factor <- down_genes_two_factor[order(down_genes_two_factor$negLog10padj), ]

# All top genes (Up and down reg.)
top_genes_two_factor <- c(top_up_two_factor$symbol, top_down_two_factor$symbol)

# Remove NA padj rows
res_valid_two_factor <- res_case_two_factor[!is.na(res_case_two_factor$padj), ]

# Get count of total DE genes, DE genes below padj < 0.05, Up-/down-reg. genes
total_genes_two_factor <- nrow(res_valid_two_factor)
number_genes_below_two_factor <- sum(res_valid_two_factor$padj < 0.005)
up_reg_genes_two_factor <- sum(res_valid_two_factor$padj < 0.005 & res_valid_two_factor$log2FoldChange > 0)
down_reg_genes_two_factor <- sum(res_valid_two_factor$padj < 0.005 & res_valid_two_factor$log2FoldChange < 0)

# Create custom volcano plot for two factor contrast (Lung WT Case vs. Lung DKO Case), for explanation of parameter/arguments used see EnhancedVolcano package
volcano_plot_two_factor <- EnhancedVolcano(
res_valid_two_factor,                                                  
pCutoff = 0.005,                                               
FCcutoff = 1.0,                                              
drawConnectors = TRUE,                                       
widthConnectors = 0.25,                                          
typeConnectors = "open",                                        
endsConnectors = "first",                                      
colConnectors = "black",                                      
lab = res_valid_two_factor$symbol,                                         
x = 'log2FoldChange',                                     
y = 'padj',                                                    
ylab= '-log10(padj)',                                        
subtitle = "Differential expression question: Does genotype (DKO vs WT) change the response to infection (T. gondii) in the lung?",
pointSize = 1.75,                                               
labSize = 1.75,                                                   
selectLab = top_genes_two_factor,                                           
legendLabels = c('Not significant: p-value > 0.005','High fold change: |Log2(FC)| ≥ 1','Significant: p-value < 0.005', 'High fold change: |Log2(FC)| ≥ 1 & Significant: p-value < 0.005'),
legendPosition = 'bottom',                                       
legendLabSize = 5,                                           
legendIconSize = 4,                                              
shape = 23,                                                     
title = "Volcano plot for contrast Lung Wild-Type Case (T. gondii) vs. Lung Double-Knockout (Ifnar−/− × Ifngr−/−) Case (T. gondii)",     # subtitle for each plot 
captionLabSize = 8,                                              
caption = bquote( #Custom caption/note to plot, Total/Up/Down regulated genes
"Total DE genes: " * .(total_genes_two_factor) * 
  "; Significant DE genes (padj < 0.005): " * .(number_genes_below_two_factor) *
  " [Up-reg.: " * .(up_reg_genes_two_factor) * ", Down-reg.: " * .(down_reg_genes_two_factor) * "]")) + 
theme( 
panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
plot.title = element_text(size = 16),      
plot.subtitle = element_text(size = 8),     
axis.title = element_text(size = 10),       
axis.text = element_text(size = 10),        
axis.title.x = element_text(size = 10, margin = margin(t = 1)),   
axis.text.y = element_text(size = 10),   
legend.text = element_text(size = 10),        
panel.grid.major = element_line(color = "grey80", size = 0.3),  
panel.grid.minor = element_line(color = "grey90", size = 0.2),   
legend.spacing.y = unit(0.5, "cm"),  
legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  
legend.box.margin = margin(0, 0, 0, 0),              
plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + 
guides(color = guide_legend(nrow = 2, byrow = TRUE)) 

# Save as PNG
ggsave(paste0(base_save_path,"volcano_plot_two_factor_plots.png"), plot = volcano_plot_two_factor, width = 16, height = 12, dpi = 300)



#=============================================================================================================
 # 6.2.1 Heatmap with genes of interest for one factor: Lung WT Case vs. Lung DKO Case 
#=============================================================================================================
 
# Get Ensemble-Id of genes of interest ( padj > 0.005, |log2FoldChange| > 1)   
genes_of_interest_one_factor <- rownames(res_case[!is.na(res_case$padj) & res_case$padj < 0.005 & abs(res_case$log2FoldChange) > 1 , ])

# Extract normalized counts
norm_counts_genes_one_factor <- counts(dds_one_factor, normalized = TRUE)

# Sorted by padj (smallest first) and then by absolute value of log2FoldChange
res_case_filtered_log2FoldChange <- res_case[order(res_case$padj, -abs(res_case$log2FoldChange)), ]

# Select top 25 genes
res_case_top25_genes <- head(rownames(res_case_filtered_log2FoldChange), 25)

# Get subset --> only genes of interest --> top 25 genes
norm_counts_genes_interest_one_factor_top_25 <- norm_counts_genes_one_factor[res_case_top25_genes, , drop = FALSE]

# Map ENSEMBL IDs to gene symbols
ens_id_gene_interest_one_factor_top_25 <- rownames(norm_counts_genes_interest_one_factor_top_25)
symbols_gene_interest_one_factor_top_25 <- mapIds(org.Mm.eg.db, keys = ens_id_gene_interest_one_factor_top_25,
                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace NAs with ENSEMBL ID (preserve names!)
symbols_gene_interest_one_factor_top_25[is.na(symbols_gene_interest_one_factor_top_25)] <- names(symbols_gene_interest_one_factor_top_25)[is.na(symbols_gene_interest_one_factor_top_25)]

# Replace rownames in matrix with gene symbols instead of Ensembl-ID
rownames(norm_counts_genes_interest_one_factor_top_25) <- symbols_gene_interest_one_factor_top_25

# Prepare sample annotation labels
sample_info$condition <- factor(
  sample_info$condition,
  levels = c("Lung_WT_Case", "Lung_DKO_Case", "Lung_WT_Control", "Lung_DKO_Control"),
  labels = c(
    "Lung Wild-Type Case (T. gondii)",
    "Lung Double-Knockout (Ifnar−/− x Ifngr−/−) Case (T. gondii)",
    "Lung Wild-Type Control",
    "Lung Double-Knockout (Ifnar−/− x Ifngr−/−) Control"
  )
)

# Define annotation colors
ann_colors <- list(
  condition = c(
    "Lung Wild-Type Case (T. gondii)" = "#FF6F61",
    "Lung Double-Knockout (Ifnar−/− x Ifngr−/−) Case (T. gondii)" = "#FFA500",
    "Lung Wild-Type Control" = "#6B5B95",
    "Lung Double-Knockout (Ifnar−/− x Ifngr−/−) Control" = "#8DD3C7"
  )
)

# Create annotation data frame
ann_df <- data.frame(
  condition = sample_info$condition)

# Ensure Row names of ann_df are names of norm_counts_genes_interest
rownames(ann_df) <- colnames(norm_counts_genes_interest_one_factor_top_25)

#Define colors for the heatmap gradient for the norm counts
col_fun <- colorRamp2(c(min(norm_counts_genes_interest_one_factor_top_25), max(norm_counts_genes_interest_one_factor_top_25)), c("#2A7BFF", "#FF4C4C"))

# Create top annotation 
top_anno <- HeatmapAnnotation(
  condition = ann_df$condition,
  col = ann_colors)

#Plot custom heatmap and save in variable heatmap_one_factor
heatmap_one_factor <- Heatmap(
  norm_counts_genes_interest_one_factor_top_25,
  name = "Expression:\nNormalized gene count",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = ann_df$Condition,
  show_column_names = TRUE,
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = top_anno,
  column_title = "Top 25 genes: Normalized gene count for genes of interest one factor (padj < 0.005 & log2FoldChange > |1|)", 
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.0f", norm_counts_genes_interest_one_factor_top_25[i, j]), 
      x, y,
      gp = gpar(fontsize = 10, col = "black")
    )
  }
)

# Save as PNG
png(paste0(base_save_path, "heatmap_gene_interest_one_factor.png"), width = 7000, height = 3500, res = 300)

# Print heatmap for top 25 genes
print(heatmap_one_factor)

# Close current graphics device
dev.off()



#=============================================================================================================
 # 6.2.2 Heatmap with genes of interest for one factor: Lung WT Case vs. Lung DKO Case 
#=============================================================================================================

# Get Ensemble-Id of genes of interest ( padj > 0.005, |log2FoldChange| > 1)   
genes_of_interest_two_factor <- rownames(res_case_two_factor[!is.na(res_case_two_factor$padj) & res_case_two_factor$padj < 0.005 & abs(res_case_two_factor$log2FoldChange) > 1 , ])

# Extract normalized counts
norm_counts_genes_two_factor <- counts(dds_two_factor, normalized = TRUE)

# Order by absolute log2FoldChange (largest first)
res_case_two_factor_filtered_log2FoldChange <- res_case_two_factor[order(res_case_two_factor$padj, -abs(res_case_two_factor$log2FoldChange)), ]

# Select top 25 genes
res_case_two_factor_top25_genes <- head(rownames(res_case_two_factor_filtered_log2FoldChange), 25)

# Get subset --> only genes of interest --> top 25 genes
norm_counts_genes_interest_two_factor_top_25 <- norm_counts_genes_two_factor[res_case_two_factor_top25_genes, , drop = FALSE]

# Map ENSEMBL IDs to gene symbols
ens_id_gene_interest_two_factor_top_25 <- rownames(norm_counts_genes_interest_two_factor_top_25)
symbols_gene_interest_two_factor_top_25 <- mapIds(org.Mm.eg.db, keys = ens_id_gene_interest_two_factor_top_25,
                                                   column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace NAs with ENSEMBL ID (preserve names!)
symbols_gene_interest_two_factor_top_25[is.na(symbols_gene_interest_two_factor_top_25)] <- names(symbols_gene_interest_two_factor_top_25)[is.na(symbols_gene_interest_two_factor_top_25)]

# Replace rownames in matrix with gene symbols instead Ensemble-ID
rownames(norm_counts_genes_interest_two_factor_top_25) <- symbols_gene_interest_two_factor_top_25

#Plot custom heatmap and save in variable heat_plot
heatmap_two_factor <- Heatmap(
  norm_counts_genes_interest_two_factor_top_25,
  name = "Expression:\nNormalized gene count",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = ann_df$Condition,
  show_column_names = TRUE,
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = top_anno,
  column_title = "Top 25 genes: Normalized gene count for genes of interest two factor (padj < 0.005 & log2FoldChange > |1|", 
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.0f", norm_counts_genes_interest_two_factor_top_25[i, j]), 
      x, y,
      gp = gpar(fontsize = 10, col = "black")
    )
  }
)

# Save as PNG
png(paste0(base_save_path, "heatmap_gene_interest_two_factor.png"), width = 7000, height = 3500, res = 300)

# Print heatmap for top 25 genes
print(heatmap_two_factor)

# Close current graphics device
dev.off()


#=============================================================================================================
 # 6.3.1  Barplot with genes of interest for one factor: Lung WT Case vs. Lung DKO Case 
#=============================================================================================================

# Reshape matrix into long format
df_melted_one_factor <- data.frame()
df_melted_one_factor <- melt(norm_counts_genes_interest_one_factor_top_25)

# Rename colnames of df_melted_one_factor
colnames(df_melted_one_factor) <- c("Gene", "Sample", "Expression")

# Get sample annotation for df from sample_info
df_melted_one_factor$Condition <- sample_info$condition[df_melted_one_factor$Sample]

# Plot custom barplot for genes of interest one factor
barplot_one_factor <- ggplot(df_melted_one_factor, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_col(position = "dodge") +
  facet_grid(~ Condition, scales = "free_x", space = "fixed") +
  theme_bw() +
  ggtitle("Top 25 genes: Normalized gene count for genes of interest one factor (padj < 0.005 & log2FoldChange > |1|)") +
  xlab("Sample") +
  ylab("Normalized gene count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, max(df_melted_one_factor$Expression), by = 25000)) +
  scale_fill_viridis_d(option = "turbo") +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.5, "lines"),
    strip.text = element_text(size = 5, colour = "black")  
  )

# Save as PNG
ggsave(paste0(base_save_path, "barplot_one_factor.png"), barplot_one_factor, width = 10, height = 6, dpi = 300)



#=============================================================================================================
 # 6.3.2 Barplot with genes of interest for one factor: Lung WT Case vs. Lung DKO Case
#=============================================================================================================

# reshape matrix into long format
df_melted_two_factor <- data.frame()
df_melted_two_factor <- melt(norm_counts_genes_interest_two_factor_top_25)

# Rename colnames of df_melted_two_factor
colnames(df_melted_two_factor) <- c("Gene", "Sample", "Expression")

# Get sample annotation for df from sample_info
df_melted_two_factor$Condition <- sample_info$condition[df_melted_two_factor$Sample]  

# Plot custom barplot for three selected genes two factor
barplot_two_factor <- ggplot(df_melted_two_factor, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_col(position = "dodge") +
  facet_grid(~ Condition, scales = "free_x", space = "fixed") +
  theme_bw() +
  ggtitle("Top 25 genes: Normalized gene count for genes of interest two factor (padj < 0.005 & log2FoldChange > |1|)") +
  xlab("Sample") +
  ylab("Normalized gene count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, max(df_melted_two_factor$Expression), by = 25000)) +
  scale_fill_viridis_d(option = "turbo") +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.5, "lines"),
    strip.text = element_text(size = 5, colour = "black")   # ← fixed
  )

# Save as PNG
ggsave(paste0(base_save_path, "barplot_two_factor.png"), barplot_two_factor, width = 10, height = 6, dpi = 300)



#=============================================================================================================
 # 6.4.1 Barplot & Heatmap with three selected genes for two factor: Lung WT Case vs. Lung DKO Case
#=============================================================================================================

# Three selected (Two factor) genes based on publication and volcano plot: 1. Tap1, 2. Ifit1, 3. Oas2
three_genes_selected = c("ENSMUSG00000032690", "ENSMUSG00000037321", "ENSMUSG00000034459")

# Get subset --> only genes of interest
norm_counts_genes_interest_two_factor_three_selected <- norm_counts_genes_two_factor[three_genes_selected, , drop = FALSE]

# Map ENSEMBL IDs to gene symbols
ens_id_gene_interest_two_factor_three_selected <- rownames(norm_counts_genes_interest_two_factor_three_selected)

symbols_gene_interest_two_factor_three_selected <- mapIds(org.Mm.eg.db, keys = ens_id_gene_interest_two_factor_three_selected,
                                                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace rownames in matrix with gene symbols instead of Ensembl-ID
rownames(norm_counts_genes_interest_two_factor_three_selected) <- symbols_gene_interest_two_factor_three_selected

#Plot custom heatmap and save in variable heatmap_two_factor_three_selected
heatmap_two_factor_three_selected <- Heatmap(
  norm_counts_genes_interest_two_factor_three_selected,
  name = "Expression:\nNormalized gene count",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = ann_df$Condition,
  show_column_names = TRUE,
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = top_anno,
  column_title = "Normalized gene count for three selected genes two factor (Tap1, Ifit1, Oas2)", 
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.0f", norm_counts_genes_interest_two_factor_three_selected[i, j]), 
      x, y,
      gp = gpar(fontsize = 10, col = "black")
    )
  }
)

# Save as PNG
png(paste0(base_save_path, "heatmap_two_factor_three_selected.png"), width = 7000, height = 3500, res = 300)

# Print heatmap for top 25 genes
print(heatmap_two_factor_three_selected)

# Close current graphics device
dev.off()



# Reshape matrix into long format
df_melted_two_factor_three_selected <- data.frame()
df_melted_two_factor_three_selected <- melt(norm_counts_genes_interest_two_factor_three_selected)

# Get sample annotation for df from sample_info
colnames(df_melted_two_factor_three_selected) <- c("Gene", "Sample", "Expression")
df_melted_two_factor_three_selected$Condition <- sample_info$condition[df_melted_two_factor_three_selected$Sample]  

# Plot custom barplot for three selected genes
barplot_three_selected_two_factor <- ggplot(df_melted_two_factor_three_selected, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_col(position = "dodge") +
  facet_grid(~ Condition, scales = "free_x", space = "fixed") +
  theme_bw() +
  ggtitle("Normalized gene count for three selected genes (Oas2, Tap1, Ifit1 ) two factor (padj < 0.005 & log2FoldChange > |1|)") +
  xlab("Sample") +
  ylab("Normalized gene count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, max(df_melted_two_factor_three_selected$Expression)+ 2000, by = 2000)) +
  scale_fill_viridis_d(option = "turbo") +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.5, "lines"),
    strip.text = element_text(size = 5, colour = "black")   # ← fixed
  )

# Save as PNG
ggsave(paste0(base_save_path, "barplot_three_selected_two_factor.png"), barplot_three_selected_two_factor, width = 10, height = 6, dpi = 300)



#-----------------------------------------------------------------------------------------------------



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 7. Over-representation analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#=============================================================================================================
#7.1 Over-representation analysis for gene of interest for two factor: Lung WT Case vs. Lung DKO Case
#=============================================================================================================

# Create new variable of res_cas (Contrast: Lung WT Case vs. Lung DKO Case)
res_case_one_factor <- res_case

# Get the Ensemble-Id from res_case_one_factor with filter padj < 0.005 & log2FoldChange > |1|
de_genes_one_factor <- rownames(res_case_one_factor[!is.na(res_case_one_factor$padj) & res_case_one_factor$padj < 0.005 & abs(res_case_one_factor$log2FoldChange) > 1,])

#Get all DE genes Ensemble ID's from res_case_one_factor --> For argument universe in enrichGO
all_genes <- rownames(res_case_one_factor)

# Get count of genes of interest
count_interest_genes <- length(de_genes_one_factor)

# Get count of total DE genes of this contrast 
count_total_universe_genes <- length(all_genes)

# Remove NA value in column log2FoldChange
res_case_one_factor <- res_case_one_factor[!is.na(res_case_one_factor$log2FoldChange), ]

# Get enrichGo terms for one factor 
ego_gene_one_factor<- enrichGO(
  gene          = de_genes_one_factor,        
  universe = all_genes,
  OrgDb         = org.Mm.eg.db,     
  ont           = "ALL",             
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.005,
  qvalueCutoff  = 1,
  pAdjustMethod = "BH",         
  readable      = TRUE    
  )

# Dotplot of top 25 GO terms
go_dotplot_one_factor <- dotplot(ego_gene_one_factor, showCategory=25)  +
              theme(axis.text.x = element_text(size = 5),  
                    axis.text.y = element_text(size = 5),  
                    axis.title.x = element_text(size = 10), 
                    axis.title.y = element_text(size = 10),
                    plot.title = element_text(size = 12, face = "bold")
                    ) + ggtitle("GO Enrichment Dotplot (Gene of interest one factor): Lung DKO Case vs. Lung WT Case") +
              labs(caption = paste("Total DE genes:", count_total_universe_genes, "| Genes of interest:", count_interest_genes))

# Save as PNG
ggsave(paste0(base_save_path, "enrichGO_dotplot_one_factor.png"), go_dotplot_one_factor, width = 10, height = 6, dpi = 300)

# Enrichment Map plot of top 25 GO terms
go_enrich_map_plot_one_factor <- emapplot(pairwise_termsim(ego_gene_one_factor)) + 
                      ggtitle("GO Enrichment Map (Gene of interest one factor): Lung DKO Case vs. Lung WT Case") +
                      theme(
                        plot.title = element_text(size = 20, face = "bold"))

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_enrichmentMap_one_factor.png"), go_enrich_map_plot_one_factor, width = 12, height = 12, dpi = 300, bg = "white")

# Gene concept Network plot of top 25 GO terms
go_gene_network_plot_one_factor <- cnetplot(ego_gene_one_factor, showCategory = 25, color_category='firebrick', 
                                 color_item="#3B78E7", howCategoryLabel = TRUE, showCategorySize = TRUE, 
                                 color_connector = 'black') + 
                    theme_void() +
                      theme(
                        plot.background = element_rect(fill = "white"),
                        plot.title = element_text(
                          size = 24,
                          face = "bold",
                          hjust = 0.5,       
                          margin = margin(b = 5, t = 5) 
                        ),
                        legend.position = "right",
                        legend.margin = margin(r = 15),  
                        legend.spacing.x = unit(0.2, "cm")  
                      ) + ggtitle("GO Enrichment Gen Network (Gene of interest one factor): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_gene_network_one_factor.png"), go_gene_network_plot_one_factor, width = 25, height = 15, dpi = 300, bg = "white")

# Order by smallest padj value and then by decreasing order of log2FoldChange
res_case_one_factor_filtered_log2FoldChange <- res_case_one_factor[order(res_case_one_factor$padj, -abs(res_case_one_factor$log2FoldChange)), ]

# Select top 25 genes of the filtered (from above)
res_case_one_factor_top25_genes <- head(rownames(res_case_one_factor_filtered_log2FoldChange), 25)

# Get top 25 genes of the filtered (1. by padj 2. log2FoldChange), get log2FoldChange value and its corresponding Ensembl-ID
geneList <- res_case_one_factor_filtered_log2FoldChange[res_case_one_factor_top25_genes,]$log2FoldChange
names(geneList) <- res_case_one_factor_top25_genes

# Get enrichGo terms for only the top 25 filtered genes, this is onyl for plotting Heatmap of most relevant ones
ego_gene_one_factor_top_25<- enrichGO(
  gene          = res_case_one_factor_top25_genes,        
  universe = all_genes,
  OrgDb         = org.Mm.eg.db,    
  ont           = "ALL",             
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.005,
  qvalueCutoff  = 1,
  pAdjustMethod = "BH",             
  readable      = TRUE             
)


# Heatmap GO terms of top 25 filtered genes 
go_heatmap_one_factor <- heatplot(
  ego_gene_one_factor_top_25,
  showCategory = 25,
  foldChange = geneList
) +
  scale_fill_gradient2(
    low = "blue",       
    mid = "grey90",      
    high = "red",
    midpoint = 0,
    limits = c(
      -max(abs(geneList), na.rm = TRUE),
      max(abs(geneList), na.rm = TRUE)
    ),
  ) +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  ggtitle("GO Enrichment Heatmap (Gene of interest one factor): Lung DKO Case vs. Lung WT Case") + labs(fill = "log2(FoldChange)")

  
# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_heatmap_one_factor.png"), go_heatmap_one_factor, width = 20, height = 6, dpi = 900)



#=============================================================================================================
#7.2  Over-representation analysis for gene of interest for two factor: Lung WT Case vs. Lung DKO Case
#=============================================================================================================

#Get the Ensemble-Id from res_case_one_factor with filter padj < 0.005 & log2FoldChange > |1|
de_genes_two_factor  <- rownames(res_case_two_factor[!is.na(res_case_two_factor$padj) & res_case_two_factor$padj < 0.005 & abs(res_case_two_factor$log2FoldChange) > 1,])

#Get all DE genes Ensemble ID's from res_case_one_factor --> For argument universe in enrichGO
all_genes_two_factor <- rownames(res_case_two_factor)

# Get count of genes of interest
count_interest_genes <- length(de_genes_two_factor)

# Get count of total DE genes of this contrast 
count_total_universe_genes <- length(all_genes_two_factor)

# Remove NA value in column log2FoldChange
res_case_two_factor <- res_case_two_factor[!is.na(res_case_two_factor$log2FoldChange), ]

# Get enrichGo terms for two factor 
ego_gene_two_factor <- enrichGO(
  gene          = de_genes_two_factor, 
  universe = all_genes,
  OrgDb         = org.Mm.eg.db,     
  ont           = "ALL",             
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",             
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE            
)

# Dotplot of top 25 GO terms
go_dotplot <- dotplot(ego_gene_two_factor, showCategory=25)  +
  theme(axis.text.x = element_text(size = 5),  
        axis.text.y = element_text(size = 5),  
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10)) + 
  ggtitle("GO Enrichment Dotplot (Gene of interest two factor): Lung DKO Case vs. Lung WT Case") +
  labs(caption = paste("Total DE genes:", count_total_universe_genes, "| Genes of interest:", count_interest_genes))


# Save as PNG
ggsave(paste0(base_save_path, "enrichGO_dotplot_two_factor.png"), go_dotplot, width = 10, height = 6, dpi = 300)

# Enrichment Map plot of top 25 GO terms
go_enrich_map_plot <- emapplot(pairwise_termsim(ego_gene_two_factor)) + 
  ggtitle("GO Enrichment Map (Gene of interest two factor): Lung DKO Case vs. Lung WT Case") +
  theme(plot.title = element_text(size = 20, face = "bold"))

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_enrichmentMap_two_factor.png"), go_enrich_map_plot, width = 20, height = 12, dpi = 300, bg = "white")

# Gene concept Network plot of top 25 GO terms
go_gene_network_plot <- cnetplot(ego_gene_two_factor, , showCategory = 25, color_category='firebrick', 
                                 color_item="#3B78E7", howCategoryLabel = TRUE, showCategorySize = TRUE, 
                                 color_connector = 'black') + 
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(
      size = 24,
      face = "bold",
      hjust = 0.5,       
      margin = margin(b = 5, t = 5) 
    ),
    legend.position = "right",
    legend.margin = margin(r = 15),  
    legend.spacing.x = unit(0.2, "cm")  
  ) + ggtitle("GO Enrichment Gen Network (Gene of interest two factor): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_gene_network_two_factor.png"), go_gene_network_plot, width = 20, height = 12, dpi = 300, bg = "white")

# Order by smallest padj value and then by decreasing order of log2FoldChange
res_case_two_factor_filtered_log2FoldChange <- res_case_two_factor[order(res_case_two_factor$padj, -abs(res_case_two_factor$log2FoldChange)), ]

# Select top 25 genes of the filtered (from above)
res_case_two_factor_top25_genes <- head(rownames(res_case_two_factor_filtered_log2FoldChange), 25)

# Get top 25 genes of the filtered (1. by padj 2. log2FoldChange), get log2FoldChange value and its corresponding Ensembl-ID
geneList_two_factor <- res_case_two_factor_filtered_log2FoldChange[res_case_two_factor_top25_genes,]$log2FoldChange
names(geneList_two_factor) <- res_case_two_factor_top25_genes


# Get enrichGo terms for only the top 25 filtered genes, this is onyl for plotting Heatmap of most relevant ones
ego_gene_two_factor_top_25<- enrichGO(
  gene          = res_case_two_factor_top25_genes,        
  universe = all_genes_two_factor,
  OrgDb         = org.Mm.eg.db,     
  ont           = "ALL",             
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.005,
  qvalueCutoff  = 1,
  pAdjustMethod = "BH",             
  readable      = TRUE              
)


# Heatmap of top 25 GO terms
go_heatmap_two_factor <- heatplot(
  ego_gene_two_factor_top_25,
  showCategory = 25,
  foldChange = geneList_two_factor
) +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey90",
    high = "red",
    midpoint = 0,
    limits = c(
      -max(abs(geneList_two_factor), na.rm = TRUE),
      max(abs(geneList_two_factor), na.rm = TRUE)
    ),
  ) +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  ggtitle("GO Enrichment Heatmap (Gene of interest two factor): Lung DKO Case vs. Lung WT Case") + labs(fill = "log2(FoldChange)")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_heatmap_two_factor.png"), go_heatmap_two_factor, width = 10, height = 6, dpi = 300)

