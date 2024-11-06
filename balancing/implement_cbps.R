rm(list = ls())
#library("sf")
#library("tigris")
# options(tigris_class = "sf")
library("tidyverse")
# library("magrittr")
library("dplyr")
source("./balancing/cbps_ATT.R")

##############---------with selected penalization parameter lambda---------##############
opt_asmd = read_csv("./latest_missing_data_fixed/balancing/sensitivity/opt_asmd.csv")[, -1]

for (i in 1:nrow(opt_asmd)) {
  year = opt_asmd$year[i]
  t = opt_asmd$t[i]
  file_directory = paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds")
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
  res_regu.list <- lapply(1:8, function(n) {
    system.time(res <- cbps_att(as.matrix(X.scl),
                                W,
                                theta.init = rep(0, ncol(X)+1),
                                #method = "Nelder-Mead",
                                control = list(trace=10, maxit = 5000),
                                lambda = rep(10^{n-t}, ncol(X)))) # reset grid search
    return(res)
  })
  converge_set = (sapply(res_regu.list, function(res) res$convergence))
  res = res_regu.list[[min(which(converge_set == 0))]]
  lambda = 10^{min(which(converge_set == 0)) - t}
  saveRDS(res, paste0("./latest_missing_data_fixed/balancing/weights-data-upd/weights-data_", year, "_tc_", lambda ,".RDS"))
}

for (year in 2005:2018) {
  file_directory = paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds")
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
  res_regu.list <- lapply(1:8, function(n) {
    system.time(res <- cbps_att(as.matrix(X.scl),
                                W,
                                theta.init = rep(0, ncol(X)+1),
                                #method = "Nelder-Mead",
                                control = list(trace=10, maxit = 5000),
                                lambda = rep(10^{n-7}, ncol(X)))) # reset grid search
    return(res)
  })
  converge_set = (sapply(res_regu.list, function(res) res$convergence))
  res = res_regu.list[[min(which(converge_set == 0))]]
  lambda = 10^{min(which(converge_set == 0)) - 7}
  saveRDS(res, paste0("./latest_missing_data_fixed/balancing/weights-data/weights-data_", year, "_tc_", lambda ,".RDS"))
}
