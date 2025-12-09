#R-Script for Step 5-7 of Workflow: 
# Exploratory data analysis, Differential expression analysis, Overrepresentation analysis.


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

#-----------------------------------------------------------------------------------------------------



#5. Exploratory data analysis

#File path to featureCount data file
file_path_featureCounts_data = "/Users/mariokummer/Desktop/RNA-Sequencing/data/DESeq_2_data/all_samples_counts.txt"

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

# Create the diffrent pastle colours for the data points in the plot
intense_pastel <- function(n) {
  base_cols <- c(
    "#FF6F61", # coral
    "#FFA500", # orange
    "#6B5B95", # purple
    "#8DD3C7"  # teal
  )
  cols <- base_cols[1:n]
  adjustcolor(cols, alpha.f = 1)
}

#Plot PCA as custome ggplot and store in variable
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
ggsave("/Users/mariokummer/Desktop/RNA-Sequencing/results/DESeq2/PCA_plot_Lung_samples.png",
       plot = PCA_plot,      # your ggplot object
       width = 8,          # in inches
       height = 6,         # in inches
       dpi = 300)    


#-----------------------------------------------------------------------------------------------------


#6. Differential Expression (DE) analysis

#================================================
  #1. Volcano plot for each experimental group
#================================================


# Extract the DE test results for each comaparison between experimental groups
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
  
  
  res$symbol <- symbols
  
  # Filter significant genes
  res$negLog10padj <- -log10(res$padj)
  sig_res <- res[!is.na(res$negLog10padj) & res$negLog10padj > padj_treshold[i], ]

  
  # Top 50 up- and down-regulated genes
  up_genes <- sig_res[sig_res$log2FoldChange > 2.5, ]
  up_genes <- as.data.frame(up_genes)
  top_up <- up_genes[order(up_genes$negLog10padj), ]
  
  down_genes <- sig_res[sig_res$log2FoldChange < - 2.5, ]
  down_genes <- as.data.frame(down_genes)
  top_down <- down_genes[order(down_genes$negLog10padj), ]
  
  # All top genes (Up and down reg.)
  top_genes <- c(top_up$symbol, top_down$symbol)
  
  top_genes
  
  # Remove NA padj rows
  res_valid <- res[!is.na(res$padj), ]
  
  # Get count of total DE genes, DE genes belo padj < 0.05, Up-/down-reg. genes
  total_genes <- nrow(res_valid)
  number_genes_below <- sum(res_valid$padj < 0.05)
  up_reg_genes <- sum(res_valid$padj < 0.05 & res_valid$log2FoldChange > 0)
  down_reg_genes <- sum(res_valid$padj < 0.05 & res_valid$log2FoldChange < 0)
  
  # Create custome volcano plot
  volcano_plot <- EnhancedVolcano(
    res_valid,
    pCutoff = 0.05,      
    FCcutoff = 1.0, 
    lab = res_valid$symbol,
    x = 'log2FoldChange',
    y = 'padj',
    subtitle = paste("Differential expression:", subtitle_text),
    pointSize = 1.25,      
    labSize = 1.5,  
    drawConnectors = TRUE,
    widthConnectors = 0.25,
    typeConnectors = "open",
    endsConnectors = "first", 
    colConnectors = "black",
    selectLab = top_genes,
    legendLabels = c('Not significant: p-value > 0.05','High fold change: |Log2(FC)| ≥ 1','Significant: p-value < 0.05', 'High fold change: |Log2(FC)| ≥ 1 & Significant: p-value < 0.05'),
    legendPosition = 'bottom',
    legendLabSize = 5,   
    legendIconSize = 4,
    shape = 23,
    title = title_text,
    caption = bquote(
      "Total DE genes: " * .(total_genes) * 
        "; Significant DE genes (padj < 0.05): " * .(number_genes_below) *
        " [Up-reg.: " * .(up_reg_genes) * ", Down-reg.: " * .(down_reg_genes) * "]"),
    captionLabSize = 8
  ) + theme(
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
    legend.spacing.y = unit(0.0005, "cm"),  
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


# Add seperation lines in 2x2 grid --> 4 boxes with seperate plot
combined_plot <- ggdraw(combined_plot) +
  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "grey80", size = 0.5) +  # vertical line
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey80", size = 0.5) +    # horizontal line
  draw_line(x = c(0, 1), y = c(0, 0),
            color = "grey60", size = 1) +
  # top border
  draw_line(x = c(0, 1), y = c(1, 1),
            color = "grey60", size = 1) +
  # left border
  draw_line(x = c(0, 0), y = c(0, 1),
            color = "grey60", size = 1) +
  # right border
  draw_line(x = c(1, 1), y = c(0, 1),
            color = "grey60", size = 1)

# Display the combined plot
print(combined_plot)

  
# Save as PNG
png_filename <- paste0("/Users/mariokummer/Desktop/RNA-Sequencing/results/DESeq2/","combined_volcano_plots.png")
ggsave(filename = png_filename, plot = combined_plot, width = 16, height = 12, dpi = 300)





#================================================
#2. Heatmap with genes of interest
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
png("/Users/mariokummer/Desktop/RNA-Sequencing/results/DESeq2/heatmap_gene_interest.png", width = 7000, height = 3500, res = 300)

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
ggsave("/Users/mariokummer/Desktop/RNA-Sequencing/results/DESeq2/barplot_gene_interest.png", barplot, width = 10, height = 6, dpi = 300)





















#7. Overrepresentation analysis
de_genes <- ens_id_gene_interest
all_genes <- rownames(res)
ego <- enrichGO(
  gene          = de_genes,        # DE genes
  OrgDb         = org.Mm.eg.db,    # mouse annotation
  ont           = "BP",             # Biological Process
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",             # Benjamini-Hochberg correction
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert IDs to gene symbols in output
)

# View top terms
head(ego)

# Number of enriched terms
nrow(ego)


# Barplot of top 10 GO terms
barplot(ego, showCategory=25, title="Top GO terms (BP)")

# Dotplot (alternative)
dotplot(ego, showCategory=25)








