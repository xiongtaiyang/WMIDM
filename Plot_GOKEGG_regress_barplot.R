
library(zlibbioc)
library(GO.db)
library(HDO.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(forcats)
library(dplyr)


# 读取CSV文件中的基因列表
gene_files <- ("E:/HCP/WM-Getm-over/Results/Gene_expression/Regress/2bk_0bk_map_r_zscore_new_gene_PLS/PLS3_geneWeights.csv")


data <- read.csv(gene_files)

# 修改变量名
names(data) <- c("gene", "fdr","zscore")

#### 差异基因提取
#提取case组较control组表达的基因，提取zscore>0.2的数据
sig_gene <- data %>%
  dplyr::filter(zscore < -2 & fdr<0.05)

#### GO富集分析
ids <- bitr(sig_gene$gene,'SYMBOL','ENTREZID','org.Hs.eg.db',drop = TRUE)
write.csv(ids, file = 'PLS3_Gene_ENTREZID_Neg.csv',row.names = F)
print(ids[1:5,])
sig_gene <- merge(sig_gene,ids,by.x='gene',by.y='SYMBOL')
sig_gene <- arrange(sig_gene, desc(zscore))
gene_diff <- sig_gene$ENTREZID
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

# Add ONTOLOGY column to KEGG data frame
kegg_df$subcategory <- NULL
kegg_df$category <- NULL
kegg_df$ONTOLOGY <- 'KEGG'
# Combine the two data frames
combined_df <- rbind(ego_df, kegg_df)

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


ego@result <- combined_df

write.csv(ego, file = 'PLS3_EnrichResult_Neg.csv',row.names = F)
nums <- nrow(combined_df)
barplot(ego,
        x = "Count",
        color = "p.adjust",
        showCategory=5,
        split='ONTOLOGY',
        label_format = Inf)+
  facet_grid(ONTOLOGY~.,
             space = 'free_y',#面板大小根据y轴自行调整
             scale='free_y'#子图坐标轴根据y轴自行调整
  )+
  #scale_y_discrete(labels=function(x) str_wrap(x, width = 100))+
  theme(text = element_text(size=24),
        axis.text.x = element_text(size = 23),  # x-axis labels
        axis.text.y = element_text(size = 23),  # y-axis labels
        axis.title = element_text(size=28),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=36,hjust = 1))

ggsave(filename = "Regress_PLS_zuhe_barplot_Neg.png",width = 10, height = 10,dpi=300)

