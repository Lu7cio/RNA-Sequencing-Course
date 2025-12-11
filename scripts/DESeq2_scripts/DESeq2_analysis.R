#R-Script for Step 5-7 of workflow: 
# 5. Exploratory data analysis, 6. Differential expression analysis, 7. Overrepresentation analysis.


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


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#5. Exploratory data analysis
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# Read FeatureCounts data
counts <- read.table(file_path_featureCounts_data, header=TRUE, row.names=1)

# Formatting: Rename Column header of dataset to only contain sample name, remove first line and the columns containing (Chr, Start, End, Strand and Length)
new_names <- sub(".*(SRR[0-9]+).*", "\\1", colnames(counts))
colnames(counts) <- new_names
counts <- counts[-1, ] 
counts <- counts[, c(0, 6:ncol(counts))]

# Add sample info metadata
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
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

#Run DESeq: Normalization + Dispersion Estimation + Differential Expression Testing
dds <- DESeq(dds)

#Regularized log transformation with rlog
#Blind = TRUE means the transformation does not use the experimental design
vsd <- rlog(dds, blind = TRUE)

# Extract PCA plot data
PCA_plot_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

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
ggsave(paste0(base_save_path,"PCA_plot_Lung_samples.png"),
       plot = PCA_plot,      # your ggplot object
       width = 8,          # in inches
       height = 6,         # in inches
       dpi = 300)    


#-----------------------------------------------------------------------------------------------------


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#6. Differential Expression (DE) analysis
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



#================================================
  #1. Volcano plot for each experimental group
#================================================

# Extract the DE test results for each comparison between experimental groups
res_wt <- results(dds, contrast = c("condition", "Lung_WT_Case", "Lung_WT_Control"))
res_dko <- results(dds, contrast = c("condition", "Lung_DKO_Case", "Lung_DKO_Control"))
res_case <- results(dds, contrast = c("condition", "Lung_DKO_Case", "Lung_WT_Case"))
res_ctrl <- results(dds, contrast = c("condition", "Lung_WT_Control", "Lung_DKO_Control"))

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

#padj value threshold for each contrast  from res_list, needed for determination which data points with label (Gen name) are shown in plot
padj_treshold = c(100,50,100,10)

# Loop through each DESeq2 result to create individual volcano plots
for (i in seq_along(res_list)) {
  
  # Get the right data set and title and subtitle
  res <- res_list[[i]]
  title_text <- custom_titles[i]
  subtitle_text <- custom_subtitles[i]
  
  # Suppose res_list[[i]] is your data frame
  res <- as.data.frame(res_list[[i]])
  
  # Map ENSEMBL IDs to gene symbols
  ens_id <- rownames(res)
  symbols <- mapIds(org.Mm.eg.db, keys = ens_id,
                    column = "SYMBOL", keytype = "ENSEMBL")
  
  # Replace NAs with ENSEMBL ID (preserve names!)
  symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
  
  # Assigne symboly to res data
  res$symbol <- symbols
  
  # Calculate the -log10 for padj for better representation in plot, assigen to new column in res
  res$negLog10padj <- -log10(res$padj)

  # Filter significant genes, not NA & negative Log10(padj) is greater than the corresponding threshold
  sig_res <- res[!is.na(res$negLog10padj) & res$negLog10padj > padj_treshold[i], ]

  # Top up- and down-regulated genes log2FoldChange > |2.5|
  up_genes <- sig_res[sig_res$log2FoldChange > 2.5, ]
  up_genes <- as.data.frame(up_genes)
  top_up <- up_genes[order(up_genes$negLog10padj), ]
  down_genes <- sig_res[sig_res$log2FoldChange < - 2.5, ]
  down_genes <- as.data.frame(down_genes)
  top_down <- down_genes[order(down_genes$negLog10padj), ]
  
  # All top genes (Up and down reg.)
  top_genes <- c(top_up$symbol, top_down$symbol)
  
  # Remove NA padj rows
  res_valid <- res[!is.na(res$padj), ]
  
  # Get count of total DE genes, DE genes below padj < 0.05, Up-/down-reg. genes
  total_genes <- nrow(res_valid)
  number_genes_below <- sum(res_valid$padj < 0.05)
  up_reg_genes <- sum(res_valid$padj < 0.05 & res_valid$log2FoldChange > 0)
  down_reg_genes <- sum(res_valid$padj < 0.05 & res_valid$log2FoldChange < 0)
  
  # Create custom volcano plot for each contrast (Lung WT Case vs. Lung DKO Case, ... etc.)
  volcano_plot <- EnhancedVolcano(
    res_valid,                                                       # data for the plot
    pCutoff = 0.05,                                                  # p-value Cutoff (pointed line in plot)      
    FCcutoff = 1.0,                                                  # log2FoldChange Cutoff (pointed line in plot)
    lab = res_valid$symbol,                                          # labels for data points (only for interesting ones)
    x = 'log2FoldChange',                                            # x-label
    y = 'padj',                                                      # y-label
    subtitle = paste("Differential expression:", subtitle_text),     # subtitle for each plot 
    pointSize = 1.25,                                                # point size of data points
    labSize = 1.5,                                                   # label size of data points labels
    drawConnectors = TRUE,                                           # enable connector from data point to label name of data point
    widthConnectors = 0.25,                                          # width of connector
    typeConnectors = "open",                                         # type of connector end
    endsConnectors = "first",                                        # where to apply end connector type
    colConnectors = "black",                                         # connector color
    selectLab = top_genes,                                           # set selection for which data points labels are shown
    legendLabels = c('Not significant: p-value > 0.05','High fold change: |Log2(FC)| ≥ 1','Significant: p-value < 0.05', 'High fold change: |Log2(FC)| ≥ 1 & Significant: p-value < 0.05'),
    legendPosition = 'bottom',                                       #Legend position
    legendLabSize = 5,                                               # Legend label size
    legendIconSize = 4,                                              #Legend icon size
    shape = 23,                                                      # data points shape
    title = title_text,                                              # title
    captionLabSize = 8,                                               # Label size for caption/note
    caption = bquote( #Custom caption/note to plot, Total/Up/Down regulated genes
      "Total DE genes: " * .(total_genes) * 
        "; Significant DE genes (padj < 0.05): " * .(number_genes_below) *
        " [Up-reg.: " * .(up_reg_genes) * ", Down-reg.: " * .(down_reg_genes) * "]")) + 
    theme( #Customization for nice plot theme/layout
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

# Save as PNG
ggsave(paste0(base_save_path,"combined_volcano_plots.png"), plot = combined_plot, width = 16, height = 12, dpi = 300)


#================================================
#2.1 Heatmap with genes of interest
#================================================

#Define some different expressed genes of interest
res_case <- results(dds, contrast = c("condition", "Lung_DKO_Case", "Lung_WT_Case"))
res_case$negLog10padj <- -log10(res_case$padj)
genes_of_interest <- rownames(res_case[!is.na(res_case$negLog10padj) & res_case$negLog10padj > 200 & res_case$log2FoldChange < -2.5, ])

genes_of_interest

# Extract normalized counts
norm_counts_genes <- counts(dds, normalized = TRUE)

# Get subset --> only genes of interest
norm_counts_genes_interest <- norm_counts_genes[genes_of_interest, , drop = FALSE]

# Map ENSEMBL IDs to gene symbols
ens_id_gene_interest <- rownames(norm_counts_genes_interest)
symbols_gene_interest <- mapIds(org.Mm.eg.db, keys = ens_id_gene_interest,
                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace rownames in matrix with gene symbols instead Ensemble-ID
rownames(norm_counts_genes_interest) <- symbols_gene_interest

# Prepare sample annotation
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

# Define annottion colours
ann_colors <- list(
  condition = c(
    "Lung Wild-Type Case (T. gondii)" = "#FF6F61",
    "Lung Double-Knockout (Ifnar−/− x Ifngr−/−) Case (T. gondii)" = "#FFA500",
    "Lung Wild-Type Control" = "#6B5B95",
    "Lung Double-Knockout (Ifnar−/− x Ifngr−/−) Control" = "#8DD3C7"
  )
)

# Create annotation dataframe
ann_df <- data.frame(
  condition = sample_info$condition)

# Rownames of ann_df should be colnames of norm_counts_genes_interest
rownames(ann_df) <- colnames(norm_counts_genes_interest)

#Define colors for the heatmap
col_fun <- colorRamp2(c(min(norm_counts_genes_interest), max(norm_counts_genes_interest)), c("#3B78E7", "#FF4C4C"))

# Create top annotation correctly
top_anno <- HeatmapAnnotation(
  condition = ann_df$condition,
  col = ann_colors)

#Plot heatmap and save in variable heat_plot
heat_plot <- Heatmap(
  norm_counts_genes_interest,
  name = "Expression:\nNormalized gene count",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = ann_df$Condition,
  show_column_names = TRUE,
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = top_anno,
  column_title = "Normalized gene count for genes of interest (padj > 200 & log2FoldChange < -2.5)", 
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.0f", norm_counts_genes_interest[i, j]), # show as integer
      x, y,
      gp = gpar(fontsize = 10, col = "black")
    )
  }
)

# Save as PNG
png(paste0(base_save_path, "heatmap_gene_interest.png"), width = 7000, height = 3500, res = 300)

# Draw heatmap
draw(heat_plot)

dev.off()


#================================================
#2.2 Heatmap 3 selected genes
#================================================

#Define the 3 selected genes 
genes_three_selected <- c("ENSMUSG00000028270", "ENSMUSG00000078853", "ENSMUSG00000046879")

# Get subset --> only 3 selected genes
norm_counts_genes_three_selected <- norm_counts_genes[genes_three_selected, , drop = FALSE]

# Map ENSEMBL IDs to gene symbols
ens_id_three_selected <- rownames(norm_counts_genes_three_selected)
symbols_three_selected <- mapIds(org.Mm.eg.db, keys = ens_id_three_selected,
                                column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace rownames in matrix with gene symbols instead Ensemble-ID
rownames(norm_counts_genes_three_selected) <- symbols_three_selected

# Rownames of ann_df should be colnames of norm_counts_genes_three_selected
rownames(ann_df) <- colnames(norm_counts_genes_three_selected)

#Plot heatmap and save in variable heat_plot
heat_plot <- Heatmap(
  norm_counts_genes_three_selected,
  name = "Expression:\nNormalized gene count",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = ann_df$Condition,
  show_column_names = TRUE,
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = top_anno,
  column_title = "Normalized gene count for 3 selected genes (Gbp2, Igtp, Irgm1)", 
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.0f", norm_counts_genes_interest[i, j]), # show as integer
      x, y,
      gp = gpar(fontsize = 10, col = "black")
    )
  }
)

# Save as PNG
png(paste0(base_save_path, "heatmap_gene_three_selected.png"), width = 7000, height = 3500, res = 300)

# Draw heatmap
draw(heat_plot)

dev.off()




#================================================
#3. Barplot with genes of interest
#================================================

# reshape matrix into long format
df <- data.frame()
df <- melt(norm_counts_genes_interest)

colnames(df) <- c("Gene", "Sample", "Expression")

df$Condition <- sample_info$condition[df$Sample]  

barplot <- ggplot(df, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_col(position = "dodge") +
  facet_grid(~ Condition, scales = "free_x", space = "fixed") +
  theme_bw() +
  ggtitle("Normalized gene count for genes of interest (padj > 200 & log2FoldChange < -2.5)") +
  xlab("Sample") +
  ylab("Normalized gene count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, max(df$Expression), by = 25000)) +
  scale_fill_viridis_d(option = "turbo") +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.5, "lines"),
    strip.text = element_text(size = 5, colour = "black")   # ← fixed
  )

# Draw barplot
print(barplot)

# Save as PNG
ggsave(paste0(base_save_path, "barplot_gene_interest.png"), barplot, width = 10, height = 6, dpi = 300)


#================================================
#3. Barplot with  3 selected genes
#================================================

# reshape matrix into long format
df <- data.frame()
df <- melt(norm_counts_genes_three_selected)

colnames(df) <- c("Gene", "Sample", "Expression")

df$Condition <- sample_info$condition[df$Sample]  

barplot <- ggplot(df, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_col(position = "dodge") +
  facet_grid(~ Condition, scales = "free_x", space = "fixed") +
  theme_bw() +
  ggtitle("Normalized gene count for 3 selected genes (Gbp2, Igtp, Irgm1)") +
  xlab("Sample") +
  ylab("Normalized gene count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, max(df$Expression), by = 25000)) +
  scale_fill_viridis_d(option = "virids") +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.5, "lines"),
    strip.text = element_text(size = 5, colour = "black")   # ← fixed
  )

# Draw barplot
print(barplot)

# Save as PNG
ggsave(paste0(base_save_path, "barplot_gene_three_selected.png"), barplot, width = 10, height = 6, dpi = 300)


#-----------------------------------------------------------------------------------------------------


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#7.1 Overrepresentation analysis for gene of interest
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# Get DE genes1
de_genes <- ens_id_gene_interest
all_genes <- rownames(res)

# Get enrichGo terms
ego_gene_interest<- enrichGO(
  gene          = de_genes,        # DE genes
  universe = all_genes,
  OrgDb         = org.Mm.eg.db,    # mouse annotation
  ont           = "BP",             # Biological Process
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",             # Benjamini-Hochberg correction
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert IDs to gene symbols in output
)

# View top terms
head(ego_gene_interest)

# Number of enriched terms
nrow(ego_gene_interest)

# Dotplot of top 25 GO terms
go_dotplot <- dotplot(ego_gene_interest, showCategory=25)  +
              theme(axis.text.x = element_text(size = 5),  
                    axis.text.y = element_text(size = 5),  
                    axis.title.x = element_text(size = 10), 
                    axis.title.y = element_text(size = 10)  
                    ) + ggtitle("GO Enrichment Dotplot (Gene of interest): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path, "enrichGO_dotplot_gene_interest.png"), go_dotplot, width = 10, height = 6, dpi = 300)

# Enrichment Map plot of top 25 GO terms
go_enrich_map_plot <- emapplot(pairwise_termsim(ego_gene_interest)) + 
                      ggtitle("GO Enrichment Map (Gene of interest): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_enrichmentMap_gene_interest.png"), go_enrich_map_plot, width = 20, height = 12, dpi = 300, bg = "white")

# Gene concept Network plot of top 25 GO terms
go_gene_network_plot <- cnetplot(ego_gene_interest, showCategory = 25, color_category='firebrick', 
                                 color_item="#3B78E7", howCategoryLabel = TRUE, showCategorySize = TRUE, 
                                 color_connector = 'black') + 
                    theme_void() +
                      theme(
                        plot.background = element_rect(fill = "white"),
                        plot.title = element_text(
                          size = 16,
                          face = "bold",
                          hjust = 0.5,       
                          margin = margin(b = 5, t = 5) 
                        ),
                        legend.position = "right",
                        legend.margin = margin(r = 15),  
                        legend.spacing.x = unit(0.2, "cm")  
                      ) + ggtitle("GO Enrichment Gen Network (Gene of interest): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_gene_network_gene_interest.png"), go_gene_network_plot, width = 20, height = 12, dpi = 300, bg = "white")

# Heatmap of top 25 GO terms
go_heatmap <- heatplot(ego_gene_interest, showCategory = 25) +
                theme(axis.text.x = element_text(size = 5),  
                      axis.text.y = element_text(size = 5),  
                      axis.title.x = element_text(size = 10), 
                      axis.title.y = element_text(size = 10))  + 
                        ggtitle("GO Enrichment Heatmap (Gene of interest): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_heatmap_gene_interest.png"), go_heatmap, width = 10, height = 6, dpi = 300)


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#7.2 Overrepresentation analysis for three selected genes
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

de_genes <- ens_id_three_selected
all_genes <- rownames(res)
ego_gene_three_selected <- enrichGO(
  gene          = de_genes, # DE genes
  universe = all_genes,
  OrgDb         = org.Mm.eg.db,    # mouse annotation
  ont           = "BP",             # Biological Process
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",             # Benjamini-Hochberg correction
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert IDs to gene symbols in output
)

# View top terms
head(ego_gene_three_selected)

# Number of enriched terms
nrow(ego_gene_three_selected)

# Dotplot of top 25 GO terms
go_dotplot <- dotplot(ego_gene_three_selected, showCategory=25)  +
  theme(axis.text.x = element_text(size = 5),  
        axis.text.y = element_text(size = 5),  
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10)  
  ) + ggtitle("GO Enrichment Dotplot (Gbp2, Igtp, Irgm1): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path, "enrichGO_dotplot_gene_three_selected.png"), go_dotplot, width = 10, height = 6, dpi = 300)

# Enrichment Map plot of top 25 GO terms
go_enrich_map_plot <- emapplot(pairwise_termsim(ego_gene_three_selected)) + 
  ggtitle("GO Enrichment Map (Gbp2, Igtp, Irgm1): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_enrichmentMap_gene_three_selected.png"), go_enrich_map_plot, width = 20, height = 12, dpi = 300, bg = "white")

# Gene concept Network plot of top 25 GO terms
go_gene_network_plot <- cnetplot(ego_gene_three_selected, showCategory = 25, color_category='firebrick', 
                                 color_item="#3B78E7", howCategoryLabel = TRUE, showCategorySize = TRUE, 
                                 color_connector = 'black') + 
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,       
      margin = margin(b = 5, t = 5) 
    ),
    legend.position = "right",
    legend.margin = margin(r = 15),  
    legend.spacing.x = unit(0.2, "cm")  
  ) + ggtitle("GO Enrichment Gen Network (Gbp2, Igtp, Irgm1): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_gene_network_gene_three_selected.png"), go_gene_network_plot, width = 20, height = 12, dpi = 300, bg = "white")

# Heatmap of top 25 GO terms
go_heatmap <- heatplot(ego_gene_three_selected, showCategory = 25) +
  theme(axis.text.x = element_text(size = 5),  
        axis.text.y = element_text(size = 5),  
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10))  + 
  ggtitle("GO Enrichment Heatmap (Gbp2, Igtp, Irgm1): Lung DKO Case vs. Lung WT Case")

# Save as PNG
ggsave(paste0(base_save_path,"enrichGO_heatmap_gene_three_selected.png"), go_heatmap, width = 10, height = 6, dpi = 300)





