# load merged data:
#install.packages("data.table")
library(data.table)
merged_counts <- readRDS('FCA_expression.dt.rds')
merged_counts
#RPMK (embryo data)
file_path <- "gene_rpkm_report_fb_2024_01.tsv"
RPKM_data <- read.table(file_path, header = TRUE, sep = "\t", col.names = c("##Release_ID", "FBgn#", "GeneSymbol", "Parent_library_FBlc#", 
                                                                            "Parent_library_name", "RNASource_FBlc#", "RNASource_name", 
                                                                            "RPKM_value", "Bin_value", "Unique_exon_base_count", 
                                                                            "Total_exon_base_count", "Count_used"))
RPKM_development <- RPKM_data[RPKM_data$Parent_library_name == "modENCODE_mRNA-Seq_development" ,]
#load promoter data
promoter_data <- read.csv("drive_performance.csv")
target_gene_FLYBASE <- promoter_data$FLYBASE
target_gene_data <- merged_counts[merged_counts$FLYBASE %in% target_gene_FLYBASE, ]
target_gene_data
target_gene_info <- target_gene_data[,1:4]
target_gene_info
target_gene <- target_gene_info$SYMBOL
target_gene
target_gene_RPLM_development <- subset(RPKM_development, FBgn. %in% target_gene_FLYBASE)
columns_to_remove <- c(1, 4, 5, 6, 9, 10, 11, 12)
target_gene_RPLM_development <- target_gene_RPLM_development[, -columns_to_remove]

library(dplyr)
library(tidyr)
grouped_data <- target_gene_RPLM_development %>%
  group_by(FBgn., GeneSymbol, .add = TRUE)
processed_data <- grouped_data %>%
  spread(RNASource_name, RPKM_value)
final_table <- processed_data %>%
  select(-FBgn.) %>%
  group_by(GeneSymbol) %>%
  summarise_all(list(~first(.))) %>%
  ungroup()
final_table <- final_table %>%
  select(FBgn., GeneSymbol,
         `mE_mRNA_em0-2hr`,       
         `mE_mRNA_em2-4hr`,       
         `mE_mRNA_em4-6hr`,       
         `mE_mRNA_em6-8hr`,      
         `mE_mRNA_em8-10hr`,     
         `mE_mRNA_em10-12hr`,     
         `mE_mRNA_em12-14hr`,     
         `mE_mRNA_em14-16hr`,     
         `mE_mRNA_em16-18hr`,     
         `mE_mRNA_em18-20hr`,     
         `mE_mRNA_em20-22hr`,     
         `mE_mRNA_em22-24hr`,     
         `mE_mRNA_L1`,            
         `mE_mRNA_L2`,            
         `mE_mRNA_L3_12hr`,       
         `mE_mRNA_L3_PS1-2`,      
         `mE_mRNA_L3_PS3-6`,      
         `mE_mRNA_L3_PS7-9`,      
         `mE_mRNA_WPP`,           
         `mE_mRNA_P5`,            
         `mE_mRNA_P6`,            
         `mE_mRNA_P8`,            
         `mE_mRNA_P9-10`,         
         `mE_mRNA_P15`,           
         `mE_mRNA_AdF_Ecl_1days`,
         `mE_mRNA_AdF_Ecl_5days`,
         `mE_mRNA_AdF_Ecl_30days`,
         `mE_mRNA_AdM_Ecl_1days`,
         `mE_mRNA_AdM_Ecl_5days`,
         `mE_mRNA_AdM_Ecl_30days`)
final_table <- as.data.frame(final_table)
merge_FT <- left_join(final_table, target_gene_info, by = c("FBgn." = "FLYBASE"))
merge_FT$GeneSymbol <- merge_FT$SYMBOL
merge_FT <- merge_FT[, -( (ncol(merge_FT) - 2):ncol(merge_FT) )]
merge_FT[, 3:ncol(merge_FT)] <- log2(merge_FT[, 3:ncol(merge_FT)] + 1)
merge_FT

cell_types <- colnames(merge_FT)
cell_types <- cell_types[-(1:2)]
cell_types
#embryo resistance Embryo>65% embryo<5% yellow/cinnabar f/cinnabar m>40
good <- c('FBgn0262598','FBgn0034724','FBgn0038204','FBgn0031074',
          'FBgn0052814','FBgn0030647','FBgn0035861','FBgn0030667',
          'FBgn0034998','FBgn0262101','FBgn0031296')
bad <- c('FBgn0037538','FBgn0005390','FBgn0004400','FBgn0033979'
         ,'FBgn0000042','FBgn0283442','FBgn0037549','FBgn0002962'
         ,'FBgn0051477','FBgn0034837','FBgn0033921','FBgn0004650'
         ,'FBgn0024177','FBgn0003943','FBgn0036809','FBgn0039585')
good
bad

good <- merge_FT[merge_FT$FBgn. %in% good,]
bad <- merge_FT[merge_FT$FBgn. %in% bad,]
# function to find p-value of a cell type:
cell.type.test <- function(query_cell_type) {
  good_values <- good[[query_cell_type]]
  bad_values <- bad[[query_cell_type]]
  utest <- wilcox.test(good_values, bad_values)#Wilcoxon rank-sum test
  p <- utest$p.value
  return(p)
}
#install.packages("pbapply")
library(pbapply)
p_values <- pblapply(cell_types, cell.type.test)
# form result:
p_values <- unlist(p_values)
result <- data.table(cell.type = cell_types,
                     p = p_values)
result <- result[order(p, decreasing = F)]
result[1:15]


library(ggplot2)
library(ggpubr)
sig_cell_types <- result[p < 0.05, cell.type]
sig_cell_types <- head(result[p < 0.05, cell.type], 9)
sig_cell_types
for (query_cell_type in sig_cell_types) {
  
  good_values <- good[[query_cell_type]]
  bad_values <- bad[[query_cell_type]]
  
  dt <- data.table(
    value = c(good_values, bad_values),
    group = factor(c(rep("good", length(good_values)), rep("bad", length(bad_values))),
                   levels = c("bad", "good"))
  )
  
  p <- ggplot(dt, aes(x = group, y = log2(value + 1), color = group)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
    scale_color_manual(values = c("good" = "#4477AA", "bad" = "#EE6677")) +
    scale_x_discrete(labels = c("good" = "Low Resistance", "bad" = "High Resistance")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 24, face = "plain", hjust = 0.5),
      axis.title.y = element_text(size = 22, face = "plain", color = "black"),
      axis.text = element_text(size = 22, color = "black"),
      axis.text.x = element_text(size = 22, face = "plain", color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = query_cell_type,
      x = NULL,
      y = expression(log[2]~"(TPM + 1)")
    )
  
  print(p)
}

dev.off()