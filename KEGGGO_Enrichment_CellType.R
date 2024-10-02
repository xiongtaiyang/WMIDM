# Load necessary libraries for GO and KEGG enrichment analysis, data manipulation, and visualization
library(zlibbioc)
library(GO.db)
library(HDO.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(forcats)
library(dplyr)

# Set the folder path where the input CSV files are located
folder_path <- ""

# Get a list of all CSV files in the folder
file_names <- list.files(path = folder_path, full.names = TRUE)

# Initialize an empty data frame to store all combined enrichment analysis results
all_combined_data <- data.frame()

# Loop over each file in the folder
for (file in file_names) {
  
  # Read the CSV file into a data frame
  data <- read.csv(file)
  
  # Rename the columns for clarity
  names(data) <- c("gene", "fdr", "zscore")
  
  # Store the data for significant genes
  sig_gene <- data
  
  #### GO enrichment analysis (Gene Ontology)
  # Convert gene symbols to ENTREZIDs for enrichment analysis
  ids <- bitr(sig_gene$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db', drop = TRUE)
  print(ids[1:5,])  # Print the first 5 rows of the converted IDs for verification
  sig_gene <- merge(sig_gene, ids, by.x = 'gene', by.y = 'SYMBOL')

  # Extract the ENTREZIDs for the significant genes
  gene_diff <- sig_gene$ENTREZID
  
  # Check if there are no ENTREZIDs; if so, exit the loop for this file
  if (is.null(gene_diff)) {
    break
  }
  
  # Perform Gene Ontology (GO) enrichment analysis
  ego <- enrichGO(gene_diff,
                  OrgDb = "org.Hs.eg.db",
                  qvalueCutoff = 0.05,
                  pvalueCutoff = 0.05,
                  ont = "all")  # 'ont' specifies which GO category (BP, CC, MF, or all)

  #### KEGG enrichment analysis (Kyoto Encyclopedia of Genes and Genomes)
  KEGG <- enrichKEGG(gene_diff,
                     organism = "hsa",  # 'hsa' for human
                     qvalueCutoff = 0.05,
                     pvalueCutoff = 0.05)
  
  # Convert enrichResult objects to data frames for GO and KEGG
  ego_df <- as.data.frame(ego)
  kegg_df <- as.data.frame(KEGG)
  
  # Combine the results from GO and KEGG into one data frame
  if (nrow(kegg_df) != 0 & nrow(ego_df) != 0) {
    # If both KEGG and GO results are available, merge them
    kegg_df$subcategory <- NULL  # Remove unnecessary columns from KEGG
    kegg_df$category <- NULL
    kegg_df$ONTOLOGY <- 'KEGG'  # Add ontology type to KEGG results
    
    # Combine the two data frames (GO and KEGG)
    combined_df <- rbind(ego_df, kegg_df)
  } else if (nrow(kegg_df) == 0 & nrow(ego_df) != 0) {
    # If only GO results are available, use GO data
    combined_df <- rbind(ego_df)
  } else if (nrow(kegg_df) != 0 & nrow(ego_df) == 0) {
    # If only KEGG results are available, use KEGG data
    combined_df <- rbind(kegg_df)
  } else {
    # If no results are available, set combined_df to NULL
    combined_df <- NULL
  }
  
  # If there are combined results, format and process them
  if (!is.null(combined_df)) {
    # Capitalize the first letter of each word in the Description column
    combined_df <- combined_df %>%
      mutate(Description = str_to_title(Description))

    # Define a custom function to capitalize the first letter of each word
    capitalize <- function(s) {
      s <- tolower(s)
      s <- strsplit(s, " ")[[1]]
      s <- paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
      return(paste(s, collapse = " "))
    }

    # Apply the capitalize function to the Description column
    combined_df$Description <- sapply(combined_df$Description, capitalize)
    
    # Add task name based on the file name
    task <- basename(file)
    task <- sub("\\.csv$", "", task)  # Remove '.csv' extension from the file name
    combined_df$task = task  # Assign task name to the new column
  }
  
  # Append the current file's results to the overall combined data frame
  all_combined_data <- rbind(all_combined_data, combined_df)

}

# Write the combined results from all files to a CSV file
write.csv(all_combined_data, file = 'PLS3_EachCellType_Pos.csv', row.names = F)

# Store the number of rows in the combined data for reference
nums <- nrow(all_combined_data)

# Update the 'ego' object with the combined results (note: this may not be necessary)
ego@result <- all_combined_data

# Code for plotting (currently not implemented)
# plot““
