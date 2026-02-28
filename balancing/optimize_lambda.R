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
# lambda_selection = data.frame(year = numeric(),
#                               lambda = numeric(), 
#                               asmd = numeric(),
#                               msmd = numeric())
# for (year in c(2005, 2009, 2010, 2011, 2013, 2014, 2016, 2017, 2018)) {
#   file_directory = paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds")
#   df = readRDS(file_directory) %>%
#     filter(fips != 42033,
#            fips != 35039)
#   X = df %>% drop_na()
#   W = X$exposed
#   X$exposed = NULL
#   X$fips = NULL
#   X.mean <- colMeans(X)
#   X.sd <- apply(X, 2, sd)
#   X.sd[X.sd == 0] <- 1 # in case Xj is constant.
#   X.scl <- scale(X, center = X.mean, scale = X.sd)
#   for (t in 1:8) {
#     res_regu.list <- lapply(1:8, function(n) {
#       system.time(res <- cbps_att(as.matrix(X.scl),
#                                   W,
#                                   theta.init = rep(0, ncol(X)+1),
#                                   #method = "Nelder-Mead",
#                                   control = list(trace=10, maxit = 5000),
#                                   lambda = rep(10^{n-t}, ncol(X)))) 
#       return(res)
#     })
#     converge_set = (sapply(res_regu.list, function(res) res$convergence))
#     res = res_regu.list[[min(which(converge_set == 0))]]
#     lambda = 10^{min(which(converge_set == 0) ) - t}
#     df_weight = df %>% drop_na()
#     df_weight$weight = df_weight$exposed *res$weights.1 + (1-df_weight$exposed)*res$weights.0
#     post_bal = NULL
#     
#     ncov = length(covariates)
#     for (i in 1:(ncov)) {
#       covariate = covariates[i]
#       post_bal[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
#                         sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
#                               subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
#         sd(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))])) 
#     }
#     lambda_selection = rbind(lambda_selection, data.frame(year = year,
#                                                           lambda = lambda,
#                                                           asmd = mean(abs(post_bal)),
#                                                           msmd = max(abs(post_bal))))
#   }
# }
# 
# # create side-by-sde 3x3 plots to showcase
# # plot each lambda against ASMD & MSMD
# plot_lambda = function(y) {
#   dat = lambda_selection %>% filter(year == y,
#                                     asmd <= 1.5)
#   ggplot(dat, aes(x = log10(lambda))) + 
#     geom_point(aes(y = asmd), size = 1.5, color = "darkblue") +
#     geom_line(aes(y = asmd), linetype = "solid", color = "darkblue") + 
#     labs(title = paste0(y), x = "log(lambda)", y = "|SMD|") + 
#     scale_y_continuous(breaks = seq(0, 1.5, by = 0.05)) +
#     scale_x_continuous(breaks = seq(-8, 0, by = 1)) + 
#     theme_minimal() +
#     geom_hline(yintercept = .1, color = "red", linetype = "dashed", size = 1)
# }
# 
# years_need_balancing = c(2005, 2009, 2010, 2011, 2013, 2014, 2016, 2017, 2018)
# plot_ls <- list()
# for (i in 1:length(years_need_balancing)) {
#   y = years_need_balancing[i]
#   plot_ls[[i]] <- plot_lambda(y)
# }
# grid.arrange(grobs = plot_ls, ncol = 3)

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

# lambda_vals = lambda_selection_upd %>% 
#   drop_na() %>%
#   pivot_longer(cols = asmd_1:msmd,
#                names_to = "measure",
#                values_to = "asmd") 
# 
# plot_asmds = function(y) {
#   dat = lambda_vals %>% filter(year == y,
#                                measure != "msmd") %>%
#     filter(asmd <= 1.5)
#   ggplot(dat, aes(x = log10(lambda), y = asmd, group = measure)) + 
#     geom_point(aes(color = measure)) +
#     geom_line(aes(color = measure)) + 
#     labs(title = paste0(y), x = "log(lambda)", y = "|SMD|") + 
#     scale_y_continuous(breaks = seq(0, 1.5, by = 0.05)) +
#     scale_x_continuous(breaks = seq(-8, 0, by = 1)) + 
#     theme_minimal() +
#     geom_hline(yintercept = .1, color = "red", linetype = "dashed", size = 1)
# }
# 
# # for the exposure years with previously no satisfying lambda values
# poorly_balanced_years = c(2009, 2010, 2014, 2016)
# plot_ls_1 <- list()
# for (i in 1:length(poorly_balanced_years)) {
#   y = poorly_balanced_years[i]
#   plot_ls_1[[i]] <- plot_asmds(y)
# }
# grid.arrange(grobs = plot_ls_1, ncol = 2)
# 
# # other exposure years
# other_years = setdiff(2005:2018, poorly_balanced_years)
# plot_ls_2 = list()
# for (i in 1:length(other_years)) {
#   y = other_years[i]
#   plot_ls_2[[i]] <- plot_asmds(y)
# }
# grid.arrange(grobs = plot_ls_2, ncol = 2)

# choose minimal ASMD_2
min_asmd = lambda_selection_upd %>% 
  drop_na() %>%
  group_by(year) %>%
  slice(which.min(asmd_2)) %>%
  select(year, lambda, asmd_2, t)

write.csv(min_asmd, "./balancing/sensitivity/opt_asmd.csv")
# 
# weights_df = data.frame(fips = numeric(),
#                         year = numeric(),
#                         weight = numeric())
# for (i in 1:nrow(min_asmd)) {
#   year = min_asmd$year[i]
#   t = min_asmd$t[i]
#   file_directory = paste0("./processed-data/df_", year, ".rds")
#   df = readRDS(file_directory) %>%
#     filter(fips != 42033,
#            fips != 35039)
#   X = df %>% drop_na()
#   W = X$exposed
#   X$exposed = NULL
#   X$fips = NULL
#   X.mean <- colMeans(X)
#   X.sd <- apply(X, 2, sd)
#   X.sd[X.sd == 0] <- 1 # in case Xj is constant.
#   X.scl <- scale(X, center = X.mean, scale = X.sd)
#   res_regu.list <- lapply(1:8, function(n) {
#     system.time(res <- cbps_att(as.matrix(X.scl),
#                                 W,
#                                 theta.init = rep(0, ncol(X)+1),
#                                 #method = "Nelder-Mead",
#                                 control = list(trace=10, maxit = 5000),
#                                 lambda = rep(10^{n-t}, ncol(X))))
#     return(res)
#   })
#   converge_set = (sapply(res_regu.list, function(res) res$convergence))
#   res = res_regu.list[[min(which(converge_set == 0))]]
#   df_weight = df %>% drop_na()
#   df_weight$weight = df_weight$exposed *res$weights.1 + (1-df_weight$exposed)*res$weights.0
#   weights_df = rbind(weights_df, data.frame(fips = df_weight$fips, 
#                                             year = year,
#                                             weight = df_weight$weight))
# }
# 
# # maximum weights by year
# max_weight = weights_df %>% 
#   group_by(year) %>%
#   slice(which.max(weight)) %>% 
#   select(year, weight)
# 
# # number of active units
# active_controls = weights_df %>%
#   filter(weight != 1) %>%
#   group_by(year) %>%
#   summarise(active_controls = sum(weight != 0))
# 
# summary_cov = merge(merge(min_asmd %>% select(year, lambda, asmd_2), 
#                           max_weight, 
#                           by = "year"), 
#                     active_controls,
#                     by = "year")
# library(knitr)
# colnames(summary_cov) = c("year", "lambda", "ASMD", "max weight", "# of active controls")
# kable(summary_cov)
# ##################----------test code for 2018------------#########################
# year = 2015
# file_directory = paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds")
# df = readRDS(file_directory) %>%
#   filter(fips != 42033,
#          fips != 35039) %>% drop_na()
# X = df %>% drop_na()
# W = X$exposed
# X$exposed = NULL
# X$fips = NULL
# X.mean <- colMeans(X)
# X.sd <- apply(X, 2, sd)
# X.sd[X.sd == 0] <- 1 # in case Xj is constant.
# X.scl <- scale(X, center = X.mean, scale = X.sd)
# lambda_selection = data.frame(lambda = numeric(), 
#                               asmd_1 = numeric(),
#                               asmd_2 = numeric(),
#                               asmd_3 = numeric(),
#                               msmd = numeric(),
#                               t = numeric())
# for (t in 1:8) {
#   res_regu.list <- lapply(1:8, function(n) {
#     system.time(res <- cbps_att(as.matrix(X.scl),
#                                 W,
#                                 theta.init = rep(0, ncol(X)+1),
#                                 #method = "Nelder-Mead",
#                                 control = list(trace=10, maxit = 5000),
#                                 lambda = rep(10^{n-t}, ncol(X)))) 
#     return(res)
#   })
#   converge_set = (sapply(res_regu.list, function(res) res$convergence))
#   res = res_regu.list[[min(which(converge_set == 0))]]
#   lambda = 10^{min(which(converge_set == 0) ) - t}
#   df_weight = df %>% drop_na()
#   df_weight$weight = df_weight$exposed *res$weights.1 + (1-df_weight$exposed)*res$weights.0
#   post_bal_1 = NULL
#   post_bal_2 = NULL
#   post_bal_3 = NULL
#   ncov = length(covariates)
#   for (i in 1:(ncov)) {
#     covariate = covariates[i]
#     post_bal_1[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
#                       sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
#                             subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
#       sd(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))])) 
#     post_bal_2[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
#                         sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
#                               subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
#       sqrt(wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
#                    subset(df_weight, exposed == 0)$weight))
#     post_bal_3[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
#                         sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
#                               subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
#       sqrt((wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
#                    subset(df_weight, exposed == 0)$weight) + 
#               var(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))) / 2)
#   }
#   lambda_selection = rbind(lambda_selection, 
#                            data.frame(lambda = lambda,
#                                       asmd_1 = mean(abs(post_bal_1)),
#                                       asmd_2 = mean(abs(post_bal_2)),
#                                       asmd_3 = mean(abs(post_bal_3)),
#                                       msmd = max(abs(post_bal_3)),
#                                       t = t))
# }
