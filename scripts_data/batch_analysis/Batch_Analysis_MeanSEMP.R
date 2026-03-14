
install.packages("emmeans")  
install.packages("Matrix", type = "binary")  
library(Matrix) 
install.packages("lme4", type = "binary")  
library(lme4)  
library(emmeans) 

data_file = "cinnabar_phenotype.csv"  # CSV files should be placed in a folder with the suffix "_drive_performance"
f = "paternal"  
expectation_value = 0.5 
output_conf_intervals_and_coefs = FALSE 


combined_data <- read.csv(data_file, as.is = TRUE, header = TRUE, check.names = FALSE, 
                          na.strings = "", blank.lines.skip = TRUE)


experiments <- unique(combined_data$Experiment)

for (experiment in experiments) {
  cur_data <- subset(combined_data, Experiment == experiment)
  group = NULL
  drive = NULL
  model = NULL
  dataframe = NULL
  

  for (i in 1:nrow(cur_data)) {
    dr = cur_data[i, 1]  
    res = cur_data[i, 2]  
    dr_vial = rep(i - 1, dr)  
    res_vial = rep(i - 1, res) 
    dr_inds = rep(1, dr)  
    res_inds = rep(0, res)  
    group = c(group, dr_vial, res_vial) 
    drive =  c(drive, dr_inds, res_inds)  
  }
  
  dataframe <- data.frame(group, drive)
  model <- glmer(drive ~ 1 + (1 | group), data = dataframe, family = binomial, nAGQ = 25)
  print(summary(model))  
  
  
  outfile <- file(paste(experiment, "_analysis.txt", sep = ""), "w")
  sink(outfile)
  print("Model summary:")
  print(summary(model))
  writeLines("\n\nCalculate expected value:")
  print(emmeans(model, ~1, type="response"))  
  writeLines("\n\nCompare the experiment to a null expectation value.")
  print(test(emmeans(model, ~1), null = qlogis(expectation_value)))  
  

  if (output_conf_intervals_and_coefs) { 
    writeLines("\n\nModel confidence intervals:")
    print(confint(model))
    writeLines("\n\nModel coefficients:")
    print(coef(model))
  }
  sink()
  close(outfile)
}

