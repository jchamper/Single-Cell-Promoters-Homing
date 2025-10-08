# load:
library(data.table)
fca <- readRDS("C:/Users/MOCHIKUKU/Desktop/promoter_selection/FCA_expression.dt.rds")
fca[1:2, 1:10]
# cell types:
cell_types <- colnames(fca)
cell_types <- cell_types[-(1:5)]
#cut rate>90 / cut rate < 25 yellow####
good <- c('FBgn0037549','FBgn0024177','FBgn0033921','FBgn0037538',
          'FBgn0034837','FBgn0031620','FBgn0051477','FBgn0028852',
          'FBgn0004650','FBgn0283442','FBgn0036809','FBgn0002962',
          'FBgn0039585','FBgn0032089','FBgn0000042','FBgn0005390',
          'FBgn0004400','FBgn0031296','FBgn0003943','FBgn0033979')
bad <- c('FBgn0262598','FBgn0036654','FBgn0037296','FBgn0262566','FBgn0052086')

#female drive conversion FDC>40 FDC<15 yellow####
good <- c('FBgn0000042','FBgn0034837','FBgn0033921','FBgn0036809'
          ,'FBgn0037538','FBgn0283442','FBgn0024177','FBgn0039585'
          ,'FBgn0003943','FBgn0002962','FBgn0004400','FBgn0031296'
          ,'FBgn0031620','FBgn0028852','FBgn0051477','FBgn0004650'
          ,'FBgn0262101','FBgn0031074','FBgn0037549','FBgn0032089'
          ,'FBgn0038204','FBgn0034998','FBgn0030647','FBgn0005390'
          ,'FBgn0003401')
bad <- c('FBgn0036654', 'FBgn0262566','FBgn0037296',
         'FBgn0052086')
#female drive conversion  FDC>40 FDC<15 cnU1####
good <- c('FBgn0031296','FBgn0037549','FBgn0037538',
          'FBgn0002962','FBgn0283442','FBgn0004650','FBgn0033921',
          'FBgn0262101','FBgn0036809','FBgn0003943',
          'FBgn0000042','FBgn0024177','FBgn0038204','FBgn0031074',
          'FBgn0031620','FBgn0030647','FBgn0039585',
          'FBgn0028852','FBgn0034837','FBgn0034998','FBgn0051477',
          'FBgn0034724','FBgn0004400','FBgn0003401','FBgn0005390',
          'FBgn0035861','FBgn0052814','FBgn0032089','FBgn0030667',
          'FBgn0033979')
bad <- c('FBgn0036654', 'FBgn0262566','FBgn0037296',
         'FBgn0052086')

#male drive conversion MDC>40,FDC<10####
good <- c('FBgn0030667','FBgn0034998','FBgn0038204','FBgn0003401','FBgn0262101',
          'FBgn0035861','FBgn0037549','FBgn0031074','FBgn0028852','FBgn0262598',
          'FBgn0031620','FBgn0039585','FBgn0024177','FBgn0032089','FBgn0031296')
bad <- c('FBgn0033979','FBgn0036654','FBgn0262566','FBgn0037296','FBgn0052086')

#embryo resistance Embryo>65% embryo<5% yellow/cinnabar f/cinnabar m>40####
good <- c('FBgn0262598','FBgn0034724','FBgn0038204','FBgn0031074',
          'FBgn0052814','FBgn0030647','FBgn0035861','FBgn0030667',
          'FBgn0034998','FBgn0262101','FBgn0031296')
bad <- c('FBgn0037538','FBgn0005390','FBgn0004400','FBgn0033979'
         ,'FBgn0000042','FBgn0283442','FBgn0037549','FBgn0002962'
         ,'FBgn0051477','FBgn0034837','FBgn0033921','FBgn0004650'
         ,'FBgn0024177','FBgn0003943','FBgn0036809','FBgn0039585')

#somatic expression  S=1/S=0 yellowFD/CMDC/CFDC>40####
good <- c('FBgn0002962','FBgn0004650','FBgn0033979','FBgn0052814',
          'FBgn0030647','FBgn0035861','FBgn0030667','FBgn0034998')
bad <- c('FBgn0037549','FBgn0037538','FBgn0283442','FBgn0004400',
         'FBgn0036809','FBgn0003943','FBgn0000042','FBgn0024177',
         'FBgn0039585','FBgn0028852','FBgn0034837','FBgn0051477')
#FDC/MDC no sig / F>M exclude low DC####
good <- c('FBgn0031296','FBgn0037549','FBgn0032089','FBgn0003401','FBgn0005390',
          'FBgn0004400','FBgn0262101','FBgn0038204','FBgn0031074','FBgn0031620',
          'FBgn0030647','FBgn0035861','FBgn0036654','FBgn0028852','FBgn0034998',
          'FBgn0052086')
bad <- c('FBgn0037538','FBgn0002962','FBgn0283442','FBgn0004650','FBgn0033979',
         'FBgn0033921','FBgn0036809','FBgn0034724','FBgn0003943','FBgn0000042',
         'FBgn0024177','FBgn0052814','FBgn0039585','FBgn0034837','FBgn0051477')

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
#Generate candidate expression plots within a single cell type.####
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
    scale_x_discrete(labels = c("good" = "High", "bad" = "Low")) +
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

#Generate volcano plots for five cell types####
library(data.table)
library(pbapply)
library(ggplot2)
library(ggrepel)
library(scales)
abbrev_map <- fread("C:/Users/MOCHIKUKU/Desktop/cell type pvalue/celltype_abbreviations.csv")
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
celltype_info <- fread("C:/Users/MOCHIKUKU/Desktop/figure/celltype4.csv")
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
  labs(title = "High - Low Female Drive Conversion Groups",
       x = expression(Delta~log[2]~"(TPM+1)"),
       y = expression(-log[10]~"(p-value)"),
       color = "Cell class") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5) ) 







