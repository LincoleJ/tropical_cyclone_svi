# compute covariate average and put in table
# pre-balancing: simple average; post-balancing: weighted mean of synthetic 
# control units
rm(list=ls())
library(dplyr)
library(tidyverse)
library(readr)

#####
####
## prepare the datasets
for (year in 2005:2018) {
  covariates <- c("exposure_", "svi_", "tmean_", "ppt_mean_", "pct_hsg_", "pct_black_", "pct_hispanic_","pci_",
                  "pct_pov_", "mf_ratio_", "pct_youth_", "pct_elderly_", "pop_density_",
                  "pop_size_")
  ncov = length(covariates)
  file1 = paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds")
  df = readRDS(file1) %>% 
    drop_na() %>%
    filter(fips != 42033,
           fips != 35039)
  d = df[, "exposed"][[1]]
  treated_units = df[d == 1, ]
  control_units = df[d == 0, ]
  X = df[, -1] %>% select(-exposed)
  X = as.matrix(X)
  X0 = t(X[d==0,]); X1 = t(X[d==1,])
  weights = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
                               pattern = paste(year),
                               full.names = TRUE))
  df_weight = df
  df_weight$weight = df_weight$exposed *weights$weights.1 + 
    (1-df_weight$exposed)*weights$weights.0
  df_weight = df_weight %>% mutate(weight = case_when(is.nan(weight) ~ 1,
                                                      TRUE ~ weight))
  treated = NULL
  untreated = NULL
  pre_bal = NULL
  post_bal = NULL
  for (i in 1:(ncov)) {
    covariate = covariates[i]
    treated[i] = mean(rowMeans(X1)[grep(covariate, names(rowMeans(X1)))])
    untreated[i] = mean(rowMeans(X0)[grep(covariate, names(rowMeans(X0)))])
    pre_bal[i] = mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df_weight))]))
    post_bal[i] = sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df_weight))]) * 
                        subset(df_weight, exposed == 0)$weight) / sum(subset(df_weight, exposed == 0)$weight)
    
  }
  
  write.csv(cbind(covariates, treated, untreated, post_bal), 
            paste0("./latest_missing_data_fixed/figures_tables/conglomerate_tables/",
                   year, "_cov_summary.csv"))
}

# for the average amoung unexposed regions

###########*\\\\\\\\
###########*
###########*
#### read 2005 and 2018
cov_2005 = read_csv("./latest_missing_data_fixed/figures_tables/conglomerate_tables/2005_cov_summary.csv")[, -1]
cov_2005

cov_2018 = read_csv("./latest_missing_data_fixed/figures_tables/conglomerate_tables/2018_cov_summary.csv")[, -1]
cov_2018
