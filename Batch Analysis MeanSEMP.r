#install.packages("lme4")
#install.packages("emmeans")
library(lme4)
library(emmeans)

setwd("/Users/weizhechen/Documents/Champerlab/Multiplexing/Figures/Batch\ effect-Figure1/")
# Configuration parameters
data_file = "combined-resistance.csv"
f="combined-resistance.csv"
expectation_value = 0.5
output_conf_intervals_and_coefs = FALSE  # Takes much longer to run.


# Read the combined CSV file
combined_data <- read.csv(data_file, as.is = TRUE, header = TRUE, check.names = FALSE,
                          na.strings = "", blank.lines.skip = TRUE)

# Extract unique experiment names
experiments <- unique(combined_data$Experiment)

# Loop through each experiment 循环遍历所有实验，对每个实验的数据进行分析
for (experiment in experiments) {
  # Subset the data for the current experiment
  cur_data <- subset(combined_data, Experiment == experiment)
  group = NULL
  drive = NULL
  model = NULL
  dataframe = NULL
  for (i in 1:nrow(cur_data)) {
    dr = cur_data[i, 1]
    res = cur_data[i, 2]
    dr_vial = rep(i - 1, dr) #生成组编号
    res_vial = rep(i - 1, res)
    dr_inds = rep(1, dr) #驱动记录1
    res_inds = rep(0, res) #非驱动记录0
    group = c(group, dr_vial, res_vial) 
    drive =  c(drive, dr_inds, res_inds) #把drive/non-drive转换成1/0
  }
  dataframe <- data.frame(group, drive)
  #glmer 拟合广义线性混合模型，二项回归
  model <- glmer(drive ~ 1 + (1 | group), data = dataframe, family = binomial, nAGQ = 25)
  print(summary(model))
  outfile <- file(paste(experiment, "_analysis.txt", sep = ""), "w")#结果写入文件
  sink(outfile)
  print("Model summary:")
  print(summary(model))
  writeLines("\n\nCalculate expected value:")
  print(emmeans(model, ~1, type="response")) #计算模型预期值
  writeLines("\n\nCompare the experiment to a null expectation value.")
  print(test(emmeans(model, ~1), null = qlogis(expectation_value))) #统计检验比较结果是否和expectation——value一致
  
  if (output_conf_intervals_and_coefs) { #输出置信区间和系数 目前是false
    writeLines("\n\nModel confidence intervals:")
    print(confint(model))
    writeLines("\n\nModel coefficients:")
    print(coef(model))
  }
  sink()
  close(outfile)
}


#在多个实验的情况下，建立 广义线性混合模型 (GLMM)，并对不同实验进行比较。下面是逐步解析：
if (length(experiments) > 1) {
  group = NULL
  drive = NULL
  Experiment = NULL
  model = NULL
  dataframe = NULL
  for (i in 1:nrow(combined_data)) {
    dr = combined_data[i, 1]
    res = combined_data[i, 2]
    exp = combined_data[i, 3]
    dr_vial = rep(i - 1, dr)
    res_vial = rep(i - 1, res)
    dr_inds = rep(1, dr)
    res_inds = rep(0, res)
    exp_dr = rep(exp, dr)
    exp_res = rep(exp, res)
    group = c(group, dr_vial, res_vial) #储存实验组编号
    drive = c(drive, dr_inds, res_inds) #储存而分类变量
    Experiment = c(Experiment, exp_dr, exp_res) #储存实验名称
  }
  dataframe <- data.frame(group, drive, as.factor(Experiment))
  print(dataframe)
  # 计算不同实验的驱动基因效应，不同实验数据可共享部分参数，提高统计效率
  model <- glmer(drive ~Experiment + (1 | group), data = dataframe, family = binomial, nAGQ = 25)
  print(summary(model))
  outfile <- file(paste(f, "_analysis.txt", sep=""), "w")
  sink(outfile)
  print("Model summary:")
  print(summary(model))
  writeLines("\n\nJoint test to see if there's any difference between any pairs of experiments.")
  print(joint_tests(model)) #检查实验组之间是否存在统计学显著差异。
  writeLines("\n\nCalculate expected value for each experiment:")
  print(emmeans(model, ~Experiment, type="response"))
  writeLines("\n\nCompare expected values for each experiment to one another as odds ratios.")
  print(emmeans(model, pairwise~Experiment, type="response")) #计算不同实验之间的比值比
  writeLines("\n\nShow diferences as a difference in proportions")
  print(emmeans(model, pairwise~Experiment, regrid="response"))   #显示实验间的差异
  writeLines("\n\nCompare each experiment to a null expectation value.")
  print(test(emmeans(model, ~Experiment), null = qlogis(expectation_value))) #比较每个实验的预期值和 理论值（expectation_value = 0.5）。
  
  if (output_conf_intervals_and_coefs) { #False here
    writeLines("\n\nModel confidence intervals:")
    print(confint(model))
    writeLines("\n\nModel coefficients:")
    print(coef(model))
  }
  sink()
  close(outfile)
}
print("Finished!")
