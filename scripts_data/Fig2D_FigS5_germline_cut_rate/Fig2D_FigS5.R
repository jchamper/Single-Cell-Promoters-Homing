# load data
library(data.table)
fca <- readRDS("FCA_expression.dt.rds")
cell_types <- colnames(fca)
cell_types <- cell_types[-(1:5)]

#cut rate>90 / cut rate < 25 yellow####
good <- c('FBgn0037549','FBgn0024177','FBgn0033921','FBgn0037538',
          'FBgn0034837','FBgn0031620','FBgn0051477','FBgn0028852',
          'FBgn0004650','FBgn0283442','FBgn0036809','FBgn0002962',
          'FBgn0039585','FBgn0032089','FBgn0000042','FBgn0005390',
          'FBgn0004400','FBgn0031296','FBgn0003943','FBgn0033979')
bad <- c('FBgn0262598','FBgn0036654','FBgn0037296','FBgn0262566','FBgn0052086')

#pvalue analysis ----
good <- fca[FLYBASE %in% good]
bad <- fca[FLYBASE %in% bad]

# function to find p-value of a cell type:
cell.type.test <- function(query_cell_type) {
  good_values <- good[[query_cell_type]]
  bad_values <- bad[[query_cell_type]]
  utest <- wilcox.test(good_values, bad_values)#Wilcoxon rank-sum test
  p <- utest$p.value
  return(p)
}
# run:
#install.packages("pbapply")
library(pbapply)
p_values <- pblapply(cell_types, cell.type.test)
# form result:
p_values <- unlist(p_values)
result <- data.table(cell.type = cell_types,
                     p = p_values)
result <- result[order(p, decreasing = F)]
result[1:15]


#FigS5 Generate candidate expression plots within a single cell type.####
library(ggplot2)
library(ggpubr)
library(stringr)
sig_cell_types <- result[p < 0.05, cell.type]
sig_cell_types <- head(result[p < 0.05, cell.type], 10)
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
    scale_x_discrete(labels = c("good" = "High cut", "bad" = "Low cut")) +
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
      title = str_wrap(query_cell_type, width = 40),
      x = NULL,
      y = expression(log[2]~"(TPM + 1)")
    )
  
  print(p)
}

dev.off()

#Fig2D Generate volcano plots for five cell types####
library(data.table)
library(pbapply)
library(ggplot2)
library(ggrepel)
library(scales)
abbrev_map <- fread("celltype_abbreviations.csv")
head(abbrev_map)

#difference of log2£¨TPM+1£©
median_safe <- function(x) median(x[is.finite(x)], na.rm = TRUE)
get_log2fc <- function(ct) {
  mg <- median_safe(good[[ct]])
  mb <- median_safe(bad[[ct]])
  log2((mg + 1) / (mb + 1))
}
log2fc_values <- sapply(cell_types, get_log2fc)
p_values <- pblapply(cell_types, cell.type.test)
p_values <- unlist(p_values)
result <- data.table(cell_type = cell_types,
                     log2FC = log2fc_values,
                     p = p_values)
result[, neg_log10_p := -log10(p)]
result[, neg_log10_p := -log10(p + 1e-300)]
result <- merge(result, abbrev_map, by = "cell_type", all.x = TRUE)
result[, sig := ifelse(p < 0.05 & abs(log2FC) > 1, "Significant", "Not Significant")]
result[, status := fifelse(p < 0.05 & log2FC > 1, "up",
                           fifelse(p < 0.05 & log2FC < -1, "down", "stable"))]
label_points <- result[sig == "Significant"]

celltype_info <- fread("celltype_annotation.csv")
colnames(celltype_info) <- c("cell_name", "cell_class")
result <- merge(result, celltype_info, by.x = "cell_type", by.y = "cell_name", all.x = TRUE)
result[, cell_class := factor(cell_class,
                              levels = c("Fgermline","Fsomatic","Mgermline","Msomatic","somatic"))]
result[, color_group := cell_class]

color_values <- c(
  "Fgermline" = "#CC79A7",  
  "Fsomatic"  = "#F7C6DE",  
  "Mgermline" = "#B36A00",  
  "Msomatic"  = "#F6C141",  
  "somatic"   = "#56B4E9"   
)
top10_labels <- result[sig == "Significant"][order(p)][1:10]
ggplot(result, aes(x = log2FC, y = neg_log10_p, color = color_group)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_text_repel(data = top10_labels,
                  aes(label = Abbreviation),
                  size = 2.5,
                  max.overlaps = Inf) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey20") +
  scale_color_manual(values = color_values, drop = FALSE) +
  labs(title = "High - Low Germline Cut Rate Groups",
       x = expression(Delta~log[2]~"(TPM+1)"),
       y = expression(-log[10]~"(p-value)"),
       color = "Cell class") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5) ) 









