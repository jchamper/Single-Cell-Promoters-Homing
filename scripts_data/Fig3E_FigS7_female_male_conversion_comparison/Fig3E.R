# load data
data <- read.csv("cnU1compare.csv")
# calculate drive inheritance
data$drive_inheritance <- data$Group.1 / (data$Group.1 + data$Group.2)
data$drive_conversion <- (data$drive_inheritance - 0.5) * 2
#install.packages("dplyr")
library(dplyr)
# generate data frame
results <- data.frame(
  Experimentname = character(),
  Test = character(),
  F_mean = numeric(),
  M_mean = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
# ttest for each experiment
for (exp in unique(data$Experimentname)) {
  subset_data <- subset(data, Experimentname == exp)
  f_vals <- subset_data$drive_conversion[subset_data$gender == "F"]
  m_vals <- subset_data$drive_conversion[subset_data$gender == "M"]
  if (length(f_vals) > 1 && length(m_vals) > 1) {
    test <- t.test(f_vals, m_vals, var.equal = FALSE)
    
    results <- rbind(results, data.frame(
      Experimentname = exp,
      Test = "Welch t-test",
      F_mean = mean(f_vals),
      M_mean = mean(m_vals),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

print(results)
FLYBASE20 <- subset(results, p_value >= 0.05)$Experimentname
FLYBASE20
#[1] "32086"  "CG4415" "CG7878" "fsM3"   "ms1"    "pp1"    "pp4"    "ps2"    "ps3"    "rcd"    "RO3"   
#[12] "shu"    "sp0"    "sp2"    "sp3_2"  "sp3"    "sp4"    "sp5"  
#cinnabar female vs cinnabar male
FLYBASE20 <- c('FBgn0031296','FBgn0037549','FBgn0032089','FBgn0003401',
               'FBgn0005390','FBgn0004400','FBgn0262101','FBgn0038204',
               'FBgn0031074','FBgn0031620','FBgn0030647','FBgn0035861',
               'FBgn0036654','FBgn0262566','FBgn0037296','FBgn0028852',
               'FBgn0034998','FBgn0052086')
# load merged data:
#install.packages("data.table")
library(data.table)
merged_counts <- readRDS('FCA_expression.dt.rds')
merged_counts
#load promoter data
promoter_data <- read.csv("drive_performance.csv")
#install.packages("dplyr")
library(dplyr)
# add SYMBOL to promoter_data
promoter_data <- promoter_data %>%
  left_join(merged_counts %>% select(FLYBASE, SYMBOL), by = "FLYBASE")
#Convert to numeric
# Column groups
cols_percent <- c(2:7, 9:12)        # Columns with percentage values 
col_numeric_only <- 8              # Column to convert directly without cleaning
# Process percentage columns (strip %, convert to numeric, divide by 100)
promoter_data[cols_percent] <- lapply(promoter_data[cols_percent], function(col) {
  clean_col <- gsub("%", "", trimws(col))           # Remove percentage sign and surrounding spaces
  clean_col <- gsub("[^0-9.-]", "", clean_col)      # Remove non-numeric characters (safe handling)
  as.numeric(clean_col) / 100
})
# Convert column 6 directly to numeric
promoter_data[[col_numeric_only]] <- as.numeric(promoter_data[[col_numeric_only]])

#rho
rho_test <- cor.test(promoter_data$cnU1_MDC, promoter_data$cnU1_FDC, method = "spearman")
print(rho_test)

matched_rows <- merged_counts[merged_counts$FLYBASE %in% FLYBASE20, ]
matched_symbols <- matched_rows$SYMBOL
matched_symbols
promoter_data$highlight_group <- ifelse(promoter_data$SYMBOL %in% matched_symbols, "matched", "other")

#figure
library(ggplot2)
library(ggrepel)
ggplot(promoter_data, aes(x = cnU1_MDC, y = cnU1_FDC, label = SYMBOL)) +
  geom_point(aes(color = highlight_group), alpha = 0.7, size = 2) +
  geom_text_repel(
    size = 3, 
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3
  ) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  annotate(
    "text",
    x = -0.05, y = 0.9,
    label = "Spearman's ¦Ñ = 0.19\np = 0.27",
    hjust = 0, size = 4
  )+
  scale_color_manual(values = c("matched" = "orange", "other" = "black")) +
  scale_x_continuous(limits = c(-0.1, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, 1), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(
    x = "cinnabar Male Drive Conversion Rate",
    y = "cinnabar Female Drive Conversion Rate",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave("cinnabar_FMplot.png", width = 6, height = 6, units = "in", dpi = 300)

