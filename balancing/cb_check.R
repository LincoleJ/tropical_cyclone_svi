print(Sys.time())
rm(list = ls())
library("tidyverse")
library("dplyr")
library("ggplot2")
library("Hmisc")

covariates <- c("exposure_", "svi_", "tmean_", "ppt_mean_", "pct_hsg_", 
                "pct_black_", "pct_hispanic_","pci_",
                "pct_pov_", "mf_ratio_", "pct_youth_", "pct_elderly_", 
                "pop_density_", "pop_size_")
cb_condition <- data.frame(covariates)
for (year in 2005:2018) {
  df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds")) %>% 
    drop_na() %>%
    filter(fips != 35039,
           fips != 42033)
  weights = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
                               pattern = paste(year),
                               full.names = TRUE))
  df_weight = df 
  df_weight$weight = df_weight$exposed *weights$weights.1 + (1-df_weight$exposed)*weights$weights.0
  post_bal <- NULL
  pre_bal <- NULL
  ncov = length(covariates)
  for (i in 1:(ncov)) {
    covariate = covariates[i]
    post_bal[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
                      sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
                            subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
      sqrt(wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
                   subset(df_weight, exposed == 0)$weight))
    pre_bal[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
                     mean(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])))/
      sqrt(wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
                   subset(df_weight, exposed == 0)$weight))
  }
  
  post_bal <- replace(post_bal, is.infinite(post_bal), 0)
  pre_bal <- replace(pre_bal, is.infinite(pre_bal), 0)
  
  order = order((pre_bal))
  balance <- data.frame(matrix(NA, nrow = ncov * 2, ncol = 0))
  balance$Covariates <- rep(1:ncov, 2)
  balance$SMD <- (c(sort((pre_bal)), post_bal[order]))
  balance$covariates_name <- rep(covariates[order], 2)
  balance$Scenarios <- c(rep("Unweighted: Orginal Data", ncov), rep("Weighted: Synthetic Control", ncov))
  balance$Scenarios <- factor(balance$Scenarios, levels = c("Unweighted: Orginal Data", "Weighted: Synthetic Control"))
  balance_p <- ggplot(balance, aes(x=SMD, y=Covariates, colour=Scenarios)) + 
    scale_y_discrete(limit = 1:ncov,
                     labels = c("Tropical Cyclone Exposure", "SVI",
                                "Summer Temperature Mean", "Summer Precipitation Mean", 
                                "Percent High School Grad",
                                "Percent Black", "Percent Hispanic", "Per Capita Income",
                                "Percent Poverty", "Male-Female Ratio", "Percent Below 15",
                                "Percent Above 65", "Population Density", "Population Size")[order]) + 
    scale_color_manual(breaks = c("Unweighted: Orginal Data", "Weighted: Synthetic Control"),
                       values=c("red", "blue")) + 
    geom_point() +
    geom_path() +
    theme_bw() +
    theme(plot.margin = unit(c(0.2, 1, 0.2, 0.2), "lines"),
          plot.title = element_text(hjust = 0.5, size = 18),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          text = element_text(size = 16),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position="bottom",
          legend.box="vertical",
          axis.text.y = element_text(angle = 30, hjust = 1)) +
    guides(color=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(paste(year)) + 
    xlab("Standardized Mean Differences") +
    xlim(min(balance$SMD) - 0.01, max(balance$SMD) + 0.01)
  ggsave(file.path("./latest_missing_data_fixed/balancing/cb_figures", paste0("Covariate_balance", year, ".jpeg")), 
         balance_p, 
         width = 8.5 / 1.6,
         height = 11 / 1.6,
         units = "in")
}

# figure that shows average |SMD|
balance_upd <- data.frame(year = NA,
                          smd = NA,
                          covariates_name = NA, 
                          scenarios = NA)
covariates_upd = covariates
ncov_upd = length(covariates)
for (year in 2005:2018) {
  df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds")) %>% 
    drop_na() %>%
    filter(fips != 35039,
           fips != 42033)
  weights = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
                               pattern = paste(year),
                               full.names = TRUE))
  df_weight = df 
  df_weight$weight = df_weight$exposed *weights$weights.1 + (1-df_weight$exposed)*weights$weights.0
  post_bal <- NULL
  pre_bal <- NULL
  ncov = length(covariates)
  for (i in 1:(ncov)) {
    covariate = covariates[i]
    post_bal[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
                      sum(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])* 
                            subset(df_weight, exposed == 0)$weight)/sum(subset(df_weight, exposed == 0)$weight))/
      sqrt(wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
                   subset(df_weight, exposed == 0)$weight))
    
    pre_bal[i] <- (mean(rowMeans(subset(df_weight, exposed == 1)[,grep(covariate, colnames(df))]))  - 
                     mean(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))])))/
      sqrt(wtd.var(rowMeans(subset(df_weight, exposed == 0)[,grep(covariate, colnames(df))]), 
                   subset(df_weight, exposed == 0)$weight))
  }
  
  post_bal <- replace(post_bal, is.infinite(post_bal), 0)
  pre_bal <- replace(pre_bal, is.infinite(pre_bal), 0)
  
  order = order((pre_bal))
  balance <- data.frame(matrix(NA, nrow = ncov_upd * 2, ncol = 0))
  balance$year = year
  balance$smd <- (c(sort((pre_bal)), post_bal[order]))
  balance$covariates_name <- rep(covariates_upd[order], 2)
  balance$scenarios <- c(rep("unweighted", ncov_upd), rep("weighted", ncov_upd))
  balance$scenarios <- factor(balance$scenarios, levels = c("unweighted", "weighted"))
  balance_upd = rbind(balance_upd, balance)
}
balance_upd = balance_upd[-1, ]
x1 = balance_upd %>% 
  filter(scenarios == "weighted") %>%
  group_by(year) %>%
  summarise(abs_smd = mean(abs(smd))) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = abs_smd)) + 
  geom_point(size = 3, shape = 15) +  
#  ggtitle("Average ASMD by year") +
  scale_x_continuous(breaks = seq(2005, 2018, by = 1)) +
  xlab("year") + ylab("average ASMD") +
  geom_hline(yintercept = .1, color = "red", linetype = "dashed", size = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  theme_minimal() #+
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank())

# plot boxplots 
x2 = balance_upd %>% filter(scenarios == "weighted") %>%
  select(year, smd) %>%
  mutate(year = as.character(year)) %>%
  ggplot() +
  geom_boxplot(aes(y = year, x = abs(smd))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = .1, color = "red", linetype = "dashed", size = 1) +
 # scale_x_continuous(limits = c(0, 0.3)) +
  labs(xlab = "ASMD", title = "Distribution of Post-Balalancing Covariates")

ggsave("./latest_missing_data_fixed/balancing/asmd_by_year_boxplot.png", x2)

# plot only dot plots
x3 = balance_upd %>% 
  filter(scenarios == "weighted") %>%
  group_by(year) %>%
  summarise(abs_smd = mean(abs(smd))) %>%
  ungroup() %>%
  ggplot(aes(y = year, x = abs_smd)) +
  geom_point(size = 5, shape = 15) + 
  scale_y_continuous(breaks = seq(2005, 2018, by = 1)) +
  geom_vline(xintercept = .1, color = "coral", linetype = "dotdash", size = 1) +
  geom_vline(xintercept = .2, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed",
                                          linewidth = 0.25),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(x = "mean ASMD",
       y = "year of exposure") 
ggsave("./latest_missing_data_fixed/figures_tables/overall_balancing.jpg", x3)

