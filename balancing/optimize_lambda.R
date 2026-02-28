rm(list=ls())
library("dplyr")
source("./balancing/cbps_ATT.R")
library(ggplot2)
library(Hmisc)
library(gridExtra)

# optimizes lambda to improve balancing conditions 
# based on optimal average |SMD|
# smd_1 = mean(X_1) - mean(X_0) / sd(X_1)
# smd_2 = mean(X_1) - mean(X_0) / wtd_sd(X_0)
# smd_3 = mean(X_1) - mean(X_0) / sqrt((var(X_1) + var(X_0)) / 2)

covariates <- c("exposure_", "svi_", "tmean_", "ppt_mean_", "pct_hsg_", 
                "pct_black_", "pct_hispanic_","pci_",
                "pct_pov_", "mf_ratio_", "pct_youth_", "pct_elderly_", 
                "pop_density_", "pop_size_")

################################################################################
#################--------------all exposure years--------------#################
################################################################################
# smd_1 = mean(X_1) - mean(X_0) / sd(X_1)
# smd_2 = mean(X_1) - mean(X_0) / wtd_sd(X_0)
# smd_3 = mean(X_1) - mean(X_0) / sqrt((var(X_1) + var(X_0)) / 2)
# asmd = mean(|smd|)

# print exposed units
for (year in 2005:2018) {
  df = readRDS(paste0("./processed-data/df_", year, ".rds")) %>%
    filter(fips != 42033, fips != 35039) %>%
    drop_na()
  
  n_exposed = sum(df$exposed == 1)
  n_control = sum(df$exposed == 0)
  
  print(paste0(year, ": exposed = ", n_exposed, ", control = ", n_control))
}

lambda_selection_upd = data.frame(year = numeric(),
                                  lambda = numeric(), 
                                  asmd_1 = numeric(),
                                  asmd_2 = numeric(),
                                  asmd_3 = numeric(),
                                  msmd = numeric(), 
                                  t = numeric())
for (year in c(2005, 2008, 2017, 2018)) {
  file_directory = paste0("./processed-data/df_", year, ".rds")
  df = readRDS(file_directory) %>%
    filter(fips != 42033,
           fips != 35039)
  X = df %>% drop_na()
  W = X$exposed
  X$exposed = NULL
  X$fips = NULL
  X.mean <- colMeans(X)
  X.sd <- apply(X, 2, sd)
  X.sd[X.sd == 0] <- 1 # in case Xj is constant.
  X.scl <- scale(X, center = X.mean, scale = X.sd)
  for (t in 1:8) {
    res_regu.list <- lapply(1:8, function(n) {
      system.time(res <- cbps_att(as.matrix(X.scl),
                                  W,
                                  theta.init = rep(0, ncol(X)+1),
                                  #method = "Nelder-Mead",
                                  control = list(trace=10, maxit = 5000),
                                  lambda = rep(10^{n-t}, ncol(X)))) 
      return(res)
    })
    converge_set = (sapply(res_regu.list, function(res) res$convergence))
    res = res_regu.list[[min(which(converge_set == 0))]]
    lambda = 10^{min(which(converge_set == 0) ) - t}
    df_weight = df %>% drop_na()
    df_weight$weight = df_weight$exposed *res$weights.1 + (1-df_weight$exposed)*res$weights.0
    post_bal_1 = NULL
    post_bal_2 = NULL
    post_bal_3 = NULL
    ncov = length(covariates)
    for (i in 1:(ncov)) {
      covariate = covariates[i]
      post_bal_1[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
                          sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
                                subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
        sd(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))])) 
      post_bal_2[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
                          sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
                                subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
        sqrt(wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
                     subset(df_weight, exposed == 0)$weight))
      post_bal_3[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
                          sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
                                subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
        sqrt((wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
                      subset(df_weight, exposed == 0)$weight) + 
                var(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))) / 2)
    }
    lambda_selection_upd = rbind(lambda_selection_upd, 
                                 data.frame(year = year,
                                            lambda = lambda,
                                            asmd_1 = mean(abs(post_bal_1)),
                                            asmd_2 = mean(abs(post_bal_2)),
                                            asmd_3 = mean(abs(post_bal_3)),
                                            msmd = max(abs(post_bal_3)),
                                            t = t))
    print(paste0(year, ", ", lambda))
  }
}
write.csv(lambda_selection_upd,
          "./balancing/sensitivity/lambda_sens_fin.csv")

# choose minimal ASMD_2
min_asmd = lambda_selection_upd %>% 
  drop_na() %>%
  group_by(year) %>%
  slice(which.min(asmd_2)) %>%
  select(year, lambda, asmd_2, t)

write.csv(min_asmd, "./balancing/sensitivity/opt_asmd.csv")
