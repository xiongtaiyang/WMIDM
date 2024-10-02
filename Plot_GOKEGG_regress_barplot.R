# Load necessary libraries
library(zlibbioc)           # BioConductor library for compressed data
library(GO.db)              # Gene Ontology database
library(HDO.db)             # Human Disease Ontology database
library(org.Hs.eg.db)       # Human gene annotation package (OrgDb)
library(clusterProfiler)    # Library for enrichment analysis
library(tidyverse)          # Collection of packages for data manipulation and visualization
library(ggplot2)            # Plotting library
library(forcats)            # For factor manipulation
library(dplyr)              # Data manipulation library

# Read the gene list from the CSV file
gene_files <- (".../PLS3_geneWeights.csv")
data <- read.csv(gene_files)

# Rename the columns of the data for clarity
names(data) <- c("gene", "fdr", "zscore")

#### Extract significant genes
# Extract genes that are significantly downregulated (zscore < -2) and fdr < 0.05
sig_gene <- data %>%
  dplyr::filter(zscore < -2 & fdr < 0.05)

#### GO enrichment analysis
# Convert gene symbols to Entrez IDs for GO analysis
ids <- bitr(sig_gene$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db', drop = TRUE)

# Save the mapping of SYMBOL to ENTREZID to a CSV file
write.csv(ids, file = 'PLS3_Gene_ENTREZID_Neg.csv', row.names = FALSE)

# Print the first 5 rows of the mapping to check the output
print(ids[1:5,])

# Merge the significant genes with their ENTREZID
sig_gene <- merge(sig_gene, ids, by.x = 'gene', by.y = 'SYMBOL')

# Sort the genes by z-score in descending order
sig_gene <- arrange(sig_gene, desc(zscore))

# Extract the ENTREZID column for GO analysis
gene_diff <- sig_gene$ENTREZID

# Perform GO enrichment analysis
ego <- enrichGO(gene_diff,
                OrgDb = "org.Hs.eg.db",
                qvalueCutoff = 0.05,
                pvalueCutoff = 0.05,
                ont = "all")

#### KEGG enrichment analysis
# Perform KEGG enrichment analysis using ENTREZID
KEGG <- enrichKEGG(gene_diff,
                   organism = "hsa",  # 'hsa' is the organism code for Homo sapiens
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05)

# Convert enrichResult objects to data frames for further processing
ego_df <- as.data.frame(ego)
kegg_df <- as.data.frame(KEGG)

# Add a new column 'ONTOLOGY' to the KEGG data to differentiate between GO and KEGG results
kegg_df$subcategory <- NULL
kegg_df$category <- NULL
kegg_df$ONTOLOGY <- 'KEGG'

# Combine GO and KEGG results into one data frame
combined_df <- rbind(ego_df, kegg_df)

# Capitalize the first letter of each word in the 'Description' column
combined_df <- combined_df %>%
  mutate(Description = str_to_title(Description))

# Custom function to capitalize the first letter of each word
capitalize <- function(s) {
  s <- tolower(s)
  s <- strsplit(s, " ")[[1]]
  s <- paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
  return(paste(s, collapse = " "))
}

# Apply the capitalization function to the 'Description' column
combined_df$Description <- sapply(combined_df$Description, capitalize)

# Update the GO enrichment result with the combined data frame
ego@result <- combined_df

# Save the combined enrichment result to a CSV file
write.csv(ego, file = 'PLS3_EnrichResult_Neg.csv', row.names = FALSE)

# Get the number of rows in the combined data frame
nums <- nrow(combined_df)

# Create a bar plot of the enrichment results, split by ONTOLOGY
barplot(ego,
        x = "Count",          # X-axis represents the count of genes in each category
        color = "p.adjust",   # Bar color based on the adjusted p-value
        showCategory = 5,     # Show top 5 categories
        split = 'ONTOLOGY',   # Split the plot by ONTOLOGY (GO vs KEGG)
        label_format = Inf) + 
  facet_grid(ONTOLOGY ~ .,   # Separate plots for each ONTOLOGY category
             space = 'free_y', # Adjust panel size based on the y-axis
             scale = 'free_y'  # Adjust y-axis scale independently for each plot
  ) +
  theme(text = element_text(size = 24),         # General text size
        axis.text.x = element_text(size = 23),  # X-axis label size
        axis.text.y = element_text(size = 23),  # Y-axis label size
        axis.title = element_text(size = 28),   # Axis title size
        legend.text = element_text(size = 20),  # Legend text size
        legend.title = element_text(size = 20), # Legend title size
        plot.title = element_text(size = 36, hjust = 1))  # Plot title size

# Save the plot as a PNG file
ggsave(filename = "Regress_PLS_zuhe_barplot_Neg.png", width = 10, height = 10, dpi = 300)
