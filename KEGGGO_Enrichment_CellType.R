library(zlibbioc)
library(GO.db)
library(HDO.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(forcats)
library(dplyr)

# 设置文件夹路径
folder_path <- "E:/HCP/WM-Getm-over/Results/Gene_expression/Regress/2bk_0bk_map_GeneResult/CellEnrichment/Pos_eachCell_Type"

# 获取文件夹中的所有文件
file_names <- list.files(path = folder_path, full.names = TRUE)

# 存储所有的富集分析结果
all_combined_data <- data.frame()

for (file in file_names) {
  
  data <- read.csv(file)
  
  # 修改变量名
  names(data) <- c("gene", "fdr","zscore")
  
  sig_gene <- data
  #### GO富集分析
  ids <- bitr(sig_gene$gene,'SYMBOL','ENTREZID','org.Hs.eg.db',drop = TRUE)
  print(ids[1:5,])
  sig_gene <- merge(sig_gene,ids,by.x='gene',by.y='SYMBOL')

  gene_diff <- sig_gene$ENTREZID
  
  if (is.null(gene_diff)){
    break
  }
  
  
  ego <- enrichGO(gene_diff,
                  OrgDb = "org.Hs.eg.db",
                  qvalueCutoff = 0.05,
                  pvalueCutoff = 0.05,
                  ont="all")
  
  #### KEGG富集分析
  KEGG <- enrichKEGG(gene_diff,
                     organism = "hsa",
                     qvalueCutoff = 0.05,
                     pvalueCutoff = 0.05)
  
  
  # 画出组合barplot富集图
  # Convert enrichResult objects to data frames
  ego_df <- as.data.frame(ego)
  kegg_df <- as.data.frame(KEGG)
  
  if (nrow(kegg_df) != 0 & nrow(ego_df) != 0){
    # Add ONTOLOGY column to KEGG data frame
    kegg_df$subcategory <- NULL
    kegg_df$category <- NULL
    kegg_df$ONTOLOGY <- 'KEGG'
  
    # Combine the two data frames
    combined_df <- rbind(ego_df, kegg_df)
  }else if ((nrow(kegg_df) == 0 & nrow(ego_df) != 0)){
    
    combined_df <- rbind(ego_df)
    
  }else if ((nrow(kegg_df) != 0 & nrow(ego_df) == 0)){
    
    combined_df <- rbind(kegg_df)
    
  }else{
    combined_df <- NULL
    
  }
  
  
  if (!is.null(combined_df)){
    # Assuming combined_df is your data frame
    combined_df <- combined_df %>%
      mutate(Description = str_to_title(Description))
  
    # Custom function to capitalize first letter of each word
    capitalize <- function(s) {
      s <- tolower(s)
      s <- strsplit(s, " ")[[1]]
      s <- paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
      return(paste(s, collapse = " "))
    }
  
    # Apply the function to the Description column
    combined_df$Description <- sapply(combined_df$Description, capitalize)
  
    # Add task
    task <- basename(file)
    task <- sub("\\.csv$", "", task)
    combined_df$task = task
  
  }
  # Append the combined_df to all_combined_data
  all_combined_data <- rbind(all_combined_data, combined_df)

}

write.csv(all_combined_data, file = 'PLS3_EachCellType_Pos.csv',row.names = F)

nums <- nrow(all_combined_data)
ego@result <- all_combined_data

# plot




