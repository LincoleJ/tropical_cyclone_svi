### conduct in-space placebo test for each exposure year
## n1 = # of exposed units per exposure year
print(Sys.time())
rm(list = ls())
library(dplyr)
library(tidyverse)
library(Hmisc)
library(splines)
source("./balancing/cbps_ATT.R")

# function to exclude rows whose fips codes are not from contiguous states
exclude_non_contiguous_states <- function(df) {
  contiguous_states <- c("01", "04", "05", "06", "08", "09", "10", "11", "12", 
                         "13", "16", "17", "18", "19", "20", "21", "22", "23", 
                         "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", 
                         "34", "35", "36", "37", "38", "39", "40", "41", "42", "44", 
                         "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", 
                         "56")
  
  # filter out rows with FIPS codes not from contiguous states
  contiguous_rows <- df[substr(df$fips, 1, 2) %in% contiguous_states, ]
  
  return(contiguous_rows)
}

# load tropical cyclone data, 
# define placebo pool to regions not exposed to any tropical cyclone 1995-2018
# exclude non-contiguous counties
hurr_dat <- read_csv("./latest_missing_data_fixed/raw-data/tc-data/hurr_dat.csv")[,-1] %>%
  exclude_non_contiguous_states() %>%
  filter(fips != 42033,
         fips != 35039)
col_names_hurr_dat = names(hurr_dat)
col_names_hurr_dat[-1] = paste0("tc_", col_names_hurr_dat[-1])
names(hurr_dat) = col_names_hurr_dat
# View(hurr_dat %>% filter(across(-1, ~ . == 0)))
placebo_pool = hurr_dat %>% filter(if_all(-1, ~. == 0))
# nrow(hurr_dat) - nrow(placebo_pool) 

covariates <- c("exposure_", "svi_", "tmean_", "ppt_mean_", "pct_hsg_", 
                "pct_black_", "pct_hispanic_","pci_",
                "pct_pov_", "mf_ratio_", "pct_youth_", "pct_elderly_", 
                "pop_density_", "pop_size_")

# parallel processing
library(parallel)
library(foreach)
library(doParallel)
num_cores <- detectCores()
print(num_cores)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# conduct analysis on all exposure years
# using parallel computing
weights_data = data.frame(exposure_year = numeric(),
                          iter = numeric(),
                          fips = numeric(),
                          placebo = numeric(),
                          weight = numeric())

for (T_0 in 2005:2018) {
  df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", T_0, ".rds")) %>% 
    drop_na() %>%
    filter(fips != 42033,
           fips != 35039)
  exposed = df$fips[df$exposed == 1]
  n1 = length(exposed)
  print(paste0("year: ", T_0))
  for (iter in 1:50) {
    set.seed(200149-iter^3)
    t1 = Sys.time()
    print(paste0("iteration: ", iter))
    print(paste0("year: ", T_0))
    placebo_pool_upd = placebo_pool[placebo_pool$fips %in% df$fips, ]
    placebo = placebo_pool_upd[sample(1:nrow(placebo_pool_upd), n1, replace = FALSE), ]$fips
    df_p = df
    df_p$placebo = 0
    df_p[df_p$fips %in% placebo, ]$placebo = 1
    df_p$exposed = NULL 
    
    # implement cbps
    X = df_p
    # X$exposed = NULL
    W = X$placebo
    X$fips = NULL
    X$placebo = NULL
    X.mean <- colMeans(X)
    X.sd <- apply(X, 2, sd)
    X.sd[X.sd == 0] <- 1 # in case Xj is constant.
    X.scl <- scale(X, center = X.mean, scale = X.sd)
    optimize_lambda = foreach(t = 1:8, 
                              .combine = rbind,
                              .packages = c("dplyr", "Hmisc")) %dopar% {
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
      df_weight = df_p
      df_weight$weight = df_weight$placebo *res$weights.1 + 
        (1-df_weight$placebo)*res$weights.0
      df_weight = df_weight %>% 
        mutate(weight = case_when(is.nan(weight) ~ 1,
                                  TRUE ~ weight))# fix infinite weights to be 1
      post_bal = NULL
      ncov = length(covariates)
      for (i in 1:(ncov)) {
        covariate = covariates[i]
        post_bal[i] <- (mean(rowMeans(subset(df_weight, placebo == 1)[,grep(covariate, colnames(df_p))]))  - 
                          sum(rowMeans(subset(df_weight, placebo == 0)[,grep(covariate, colnames(df_p))])* 
                                subset(df_weight, placebo == 0)$weight)/sum(subset(df_weight, placebo == 0)$weight))/
          sqrt(wtd.var(rowMeans(subset(df_weight, placebo == 0)[,grep(covariate, colnames(df))]), 
                       subset(df_weight, placebo == 0)$weight))
      }
      data.frame(t = t, asmd = mean(abs(post_bal)))
    }
    optimal_t = optimize_lambda$t[which.min(optimize_lambda$asmd)]
    res_regu.list <- lapply(1:8, function(n) {
      system.time(res <- cbps_att(as.matrix(X.scl),
                                  W,
                                  theta.init = rep(0, ncol(X)+1),
                                  #method = "Nelder-Mead",
                                  control = list(trace=10, maxit = 5000),
                                  lambda = rep(10^{n - optimal_t}, ncol(X))))
      return(res)
    })
    converge_set = (sapply(res_regu.list, function(res) res$convergence))
    res = res_regu.list[[min(which(converge_set == 0))]]
    lambda = 10^{min(which(converge_set == 0)) - optimal_t}
    df_weight = df_p
    df_weight$weight = df_weight$placebo *res$weights.1 + 
      (1-df_weight$placebo)*res$weights.0
    df_weight = df_weight %>% 
      mutate(weight = case_when(is.nan(weight) ~ 1,
                                TRUE ~ weight))# fix infinite weights to be 1
    print(paste0("for each iteration: ", Sys.time() - t1, " seconds"))
    weights_data = rbind(weights_data, data.frame(exposure_year = T_0,
                                                  iter = iter,
                                                  fips = df_weight$fips,
                                                  placebo = df_weight$placebo,
                                                  weight = df_weight$weight))
    }
}
stopImplicitCluster()
write.csv(weights_data,
          "./latest_missing_data_fixed/analysis/in_space_placebo/weights_data_(set_seed).csv")

###################################################################################
## test code for 2005
# 不用parallel computing
time1 = Sys.time()
weights_data = data.frame(exposure_year = numeric(),
                          iter = numeric(),
                          fips = numeric(),
                          placebo = numeric(),
                          weight = numeric())
T_0 = 2005
df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", T_0, ".rds")) %>% 
  drop_na() %>%
  filter(fips != 42033,
         fips != 35039)
exposed = df$fips[df$exposed == 1]
n1 = length(exposed)
print(paste0("year: ", T_0))
for (iter in 1:2) {
  print(paste0("iteration: ", iter))
  print(paste0("year: ", T_0))
  placebo_pool_upd = placebo_pool[placebo_pool$fips %in% df$fips, ]
  placebo = placebo_pool_upd[sample(1:nrow(placebo_pool_upd), n1, replace = FALSE), ]$fips
  df_p = df
  df_p$placebo = 0
  df_p[df_p$fips %in% placebo, ]$placebo = 1
  df_p$exposed = NULL 
  
  # implement cbps
  X = df_p
 # X$exposed = NULL
  W = X$placebo
  X$fips = NULL
  X$placebo = NULL
  X.mean <- colMeans(X)
  X.sd <- apply(X, 2, sd)
  X.sd[X.sd == 0] <- 1 # in case Xj is constant.
  X.scl <- scale(X, center = X.mean, scale = X.sd)
  optimize_lambda = data.frame(t = numeric(),
                              asmd = numeric())
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
    df_weight = df_p
    df_weight$weight = df_weight$placebo *res$weights.1 + 
      (1-df_weight$placebo)*res$weights.0
    df_weight = df_weight %>% 
      mutate(weight = case_when(is.nan(weight) ~ 1,
                                TRUE ~ weight))# fix infinite weights to be 1
    post_bal = NULL
    ncov = length(covariates)
    for (i in 1:(ncov)) {
      covariate = covariates[i]
      post_bal[i] <- (mean(rowMeans(subset(df_weight, placebo == 1)[,grep(covariate, colnames(df_p))]))  - 
                        sum(rowMeans(subset(df_weight, placebo == 0)[,grep(covariate, colnames(df_p))])* 
                              subset(df_weight, placebo == 0)$weight)/sum(subset(df_weight, placebo == 0)$weight))/
        sqrt(wtd.var(rowMeans(subset(df_weight, placebo == 0)[,grep(covariate, colnames(df))]), 
                     subset(df_weight, placebo == 0)$weight))
    }
    optimize_lambda = rbind(optimize_lambda, data.frame(t = t,
                                                        asmd = mean(abs(post_bal))))
  }
  optimal_t = optimize_lambda$t[which.min(optimize_lambda$asmd)]
  res_regu.list <- lapply(1:8, function(n) {
    system.time(res <- cbps_att(as.matrix(X.scl),
                                W,
                                theta.init = rep(0, ncol(X)+1),
                                #method = "Nelder-Mead",
                                control = list(trace=10, maxit = 5000),
                                lambda = rep(10^{n - optimal_t}, ncol(X))))
    return(res)
  })
  converge_set = (sapply(res_regu.list, function(res) res$convergence))
  res = res_regu.list[[min(which(converge_set == 0))]]
  lambda = 10^{min(which(converge_set == 0)) - optimal_t}
  df_weight = df_p
  df_weight$weight = df_weight$placebo *res$weights.1 + 
    (1-df_weight$placebo)*res$weights.0
  df_weight = df_weight %>% 
    mutate(weight = case_when(is.nan(weight) ~ 1,
                              TRUE ~ weight))# fix infinite weights to be 1
  weights_data = rbind(weights_data, data.frame(exposure_year = T_0,
                                                iter = iter,
                                                fips = df_weight$fips,
                                                placebo = df_weight$placebo,
                                                weight = df_weight$weight))
}
time2 = Sys.time()
time2 - time1

# full data frame
for (T_0 in 2005:2018) {
  df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", T_0, ".rds")) %>% 
    drop_na() %>%
    filter(fips != 42033,
           fips != 35039)
  exposed = df$fips[df$exposed == 1]
  n1 = length(exposed)
  print(paste0("year: ", T_0))
  for (iter in 1:100) {
    print(paste0("iteration: ", iter))
    print(paste0("year: ", T_0))
    placebo_pool_upd = placebo_pool[placebo_pool$fips %in% df$fips, ]
    placebo = placebo_pool_upd[sample(1:nrow(placebo_pool_upd), n1, replace = FALSE), ]$fips
    
    df$placebo = 0
    df[df$fips %in% placebo, ]$placebo = 1
    df$exposed = NULL
    
    # implement cbps
    X = df
    W = X$placebo
    X$fips = NULL
    X$placebo = NULL
    X.mean <- colMeans(X)
    X.sd <- apply(X, 2, sd)
    X.sd[X.sd == 0] <- 1 # in case Xj is constant.
    X.scl <- scale(X, center = X.mean, scale = X.sd)
    res_regu.list <- lapply(1:8, function(n) {
      system.time(res <- cbps_att(as.matrix(X.scl),
                                  W,
                                  theta.init = rep(0, ncol(X)+1),
                                  #method = "Nelder-Mead",
                                  control = list(trace=10, maxit = 5000),
                                  lambda = rep(10^{n-7}, ncol(X))))
      return(res)
    })
    converge_set = (sapply(res_regu.list, function(res) res$convergence))
    res = res_regu.list[[min(which(converge_set == 0))]]
    lambda = 10^{min(which(converge_set == 0)) - 7}
    df_weight = df
    df_weight$weight = df_weight$placebo *res$weights.1 + 
      (1-df_weight$placebo)*res$weights.0
    df_weight = df_weight %>% mutate(weight = case_when(is.nan(weight) ~ 1,
                                                        TRUE ~ weight)) # fix infinite weights to be 1
    control = df_weight$fips[df_weight$weight != 0 & df_weight$weight != 1]
    df_final = df_weight[df_weight$fips %in% c(control, placebo), ] %>% 
      select(fips, weight)
    
    svi = svi_grouped[svi_grouped$fips %in% c(control, placebo), ]
    svi_final = merge(df_final, svi, by = "fips")
    y = numeric()
    for (i in 1:(ncol(svi) - 1)) {
      y[i] = mean(subset(svi_final, weight == 1)[, i + 2]) - 
        (sum(subset(svi_final, weight != 1)$weight * 
               subset(svi_final, weight != 1)[, i + 2]) /
           sum(subset(svi_final, weight != 1)$weight))
    }
    
    pop_size = pop_grouped[pop_grouped$fips %in% c(control, placebo), ]
    pop_final = merge(df_final, pop_size, by = "fips")
    
    wtd_pop_size = sum(subset(pop_final, weight != 1)[,grep(paste0("pop_", T_0), names(pop_final))] *
                         subset(pop_final, weight != 1)$weight) / 
      sum(subset(pop_final, weight != 1)$weight)
    
    result = data.frame(exposure_year = T_0,
                        acee_2010 = y[1],
                        acee_2014 = y[2],
                        acee_2016 = y[3],
                        acee_2018 = y[4],
                        acee_2020 = y[5],
                        wtd_pop_size = wtd_pop_size,
                        iter = iter)
    
    results = rbind(results, result)
  }
}
