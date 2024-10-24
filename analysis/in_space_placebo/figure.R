# calculate ACEE, weighted population size at time of exposure
# plot acee v.s. for each iteration
rm(list=ls())
library(dplyr)
library(tidyverse)
library(Hmisc)
library(splines)

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

# import weights data
weights_data = read_csv("./latest_missing_data_fixed/analysis/in_space_placebo/weights_data_(set_seed).csv")[, -1]

################################################################################ 
#########################---------main figure---------##########################
################################################################################ 
# import outcome variables (i.e., SVI from 2010, 2014, 2016, 2018, 2020)
svi_grouped = read_csv("./latest_missing_data_fixed/raw-data/svi-data/svi_dat.csv")[-c(1, 3)] %>%
  filter(fips != 42033,
         fips != 35039)
svi_colnames = c("fips", paste0("svi_", colnames(svi_grouped)[-1]))
colnames(svi_grouped) = svi_colnames

# import population size variable
pop_grouped = list.files(path = "./raw-data/census-data/cdc_population_monthly_infer/pct_hispanic",
                         pattern = "*.csv", full.name = TRUE) %>% 
  read_csv() %>% 
  filter(month == 6) %>% 
  select(-month) %>% 
  group_by(fips, year) %>%
  summarise(pop = sum(pop), .groups = "keep") %>%
  pivot_wider(names_from = year,
              values_from = pop,
              names_prefix = "pop_") %>%
  exclude_non_contiguous_states() %>%
  filter(fips != 42033,
         fips != 35039)

# calculate ACEE
results = data.frame(exposure_year = numeric(),
                     acee_2010 = numeric(),
                     acee_2014 = numeric(),
                     acee_2016 = numeric(),
                     acee_2018 = numeric(),
                     acee_2020 = numeric(),
                     wtd_pop_size = numeric(),
                     iter = numeric())
for (T_0 in 2005:2018) {
  for(iteration in 1:50) {
    df_weight = weights_data %>% filter(exposure_year == T_0,
                                        iter == iteration)
    placebo = subset(df_weight, df_weight$placebo == 1)$fips
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
    
    # weighted population size at time of exposure
    pop_size = pop_grouped[pop_grouped$fips %in% c(control, placebo), ]
    pop_final = merge(df_final, pop_size, by = "fips")
    wtd_pop_size = sum(subset(pop_final, weight != 1)[,grep(paste0("pop_", T_0), names(pop_final))] *
                         subset(pop_final, weight != 1)$weight) / 
      sum(subset(pop_final, weight != 1)$weight)
    
    results = rbind(results, data.frame(exposure_year = T_0,
                                        acee_2010 = y[1],
                                        acee_2014 = y[2],
                                        acee_2016 = y[3],
                                        acee_2018 = y[4],
                                        acee_2020 = y[5],
                                        wtd_pop_size = wtd_pop_size,
                                        iter = iteration))
  }
}

# outcome model for each placebo run.
all.lags = 1:15
all.outcome.years = c(2010, 2014, 2016, 2018, 2020)
all.lags.ns = ns(all.lags, df = 1)
cmp.ns = ggplot()
for (i in 1:49) {
  results_summary = results %>% 
    filter(iter == i) %>%
    select(-iter)
  acee_pe_placebo = results_summary %>%
    pivot_longer(cols = starts_with("acee_"), 
                 names_to = "outcome_year",
                 names_prefix = "acee_",
                 values_to = "acee") %>%
    mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
    filter(lag > 0)
  jackfun.placebo = function(end.year) {
    jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
    reg.jack.placebo = lm(acee ~ ns(lag, df = 1),
                          weights = wtd_pop_size,
                          data = jack.data.placebo)
    coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
  }
  full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
  mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
    mean(sapply(full.reg.placebo, function(v) v[i]))
  })
  cmp.ns = cmp.ns + 
    geom_line(data = data.frame(Lag = all.lags, Diff = mean.reg.placebo),
              aes(x=Lag, y = Diff), 
              size = 0.4,
              color = "grey")
}

i = 50
results_summary = results %>% 
  filter(iter == i) %>%
  select(-iter)
acee_pe_placebo = results_summary %>%
  pivot_longer(cols = starts_with("acee_"), 
               names_to = "outcome_year",
               names_prefix = "acee_",
               values_to = "acee") %>%
  mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
  filter(lag > 0)
jackfun.placebo = function(end.year) {
  jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
  reg.jack.placebo = lm(acee ~ ns(lag, df = 1),
                        weights = wtd_pop_size,
                        data = jack.data.placebo)
  coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
}
full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
  mean(sapply(full.reg.placebo, function(v) v[i]))
})
plot_placebo_df = data.frame(Lag = all.lags, Diff = mean.reg.placebo, type = "placebo controls")
plot_placebo_df = rbind(plot_placebo_df, data.frame(Lag = all.lags, Diff = mean.reg, type = "exposed"))
cmp.ns.placebo = cmp.ns + 
  geom_line(data = plot_placebo_df,
            aes(x=Lag, y = Diff, group = type, color = type, size = type)) +
  scale_color_manual(values = c("placebo controls" = "grey", "exposed" = "red"), 
                     name = NULL) +
  scale_size_manual(values = c("placebo controls" = 0.4, "exposed" = 1.25), 
                    name = NULL) +
  labs(x = "Lag", y = "Weighted Difference in SVI") +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13, 15)) +
  theme_bw() +
  theme(legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.background = element_rect(fill = "white", colour = "black"), # box around legend
        legend.key = element_rect(fill = "transparent"))

ggsave("./latest_missing_data_fixed/figures_tables/placebo_main.jpg", cmp.ns.placebo)

#########---------for subcategories of SVI-------
for (svi_subcat in c("housing", "ses", "minority", "hhd")) {
  svi_file_name = paste0("./latest_missing_data_fixed/raw-data/svi-data/svi_dat_", 
                         svi_subcat, ".csv")
  svi_grouped = read_csv(svi_file_name)[-c(1, 3)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames = c("fips", paste0("svi_", colnames(svi_grouped)[-1]))
  colnames(svi_grouped) = svi_colnames
  
  svi_subcat_results = data.frame(exposure_year = numeric(),
                                  acee_2010 = numeric(),
                                  acee_2014 = numeric(),
                                  acee_2016 = numeric(),
                                  acee_2018 = numeric(),
                                  acee_2020 = numeric(),
                                  wtd_pop_size = numeric(),
                                  iter = numeric())
  for (T_0 in 2005:2018) {
    for(iteration in 1:50) {
      df_weight = weights_data %>% filter(exposure_year == T_0,
                                          iter == iteration)
      placebo = subset(df_weight, df_weight$placebo == 1)$fips
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
      
      # weighted population size at time of exposure
      pop_size = pop_grouped[pop_grouped$fips %in% c(control, placebo), ]
      pop_final = merge(df_final, pop_size, by = "fips")
      wtd_pop_size = sum(subset(pop_final, weight != 1)[,grep(paste0("pop_", T_0), names(pop_final))] *
                           subset(pop_final, weight != 1)$weight) / 
        sum(subset(pop_final, weight != 1)$weight)
      
      svi_subcat_results = rbind(svi_subcat_results, 
                                 data.frame(exposure_year = T_0,
                                            acee_2010 = y[1],
                                            acee_2014 = y[2],
                                            acee_2016 = y[3],
                                            acee_2018 = y[4],
                                            acee_2020 = y[5],
                                            wtd_pop_size = wtd_pop_size,
                                            iter = iteration))
    }
  }
  write.csv(svi_subcat_results, 
            paste0("./latest_missing_data_fixed/analysis/in_space_placebo/subcat_svi/",
                   svi_subcat, "_results.csv"))
}
# 
# for (svi_subcat in c("housing", "ses", "minority", "hhd")) {
#   cmp.placebo.subcat = ggplot()
#   file_name = paste0("latest_missing_data_fixed/analysis/in_space_placebo/subcat/svi", svi_subcat, "_results.csv")
#   svi_subcat_results = read_csv(file_name)
#   for (i in 1:50) {
#     results_summary = svi_subcat_results %>% filter(iter == i) %>%
#       select(-iter)
#     acee_pe_placebo = results_summary %>% 
#       pivot_longer(cols = starts_with("acee_"),
#                    names_to = "outcome_year",
#                    names_prefix = "acee_",
#                    values_to = "acee") %>%
#       mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
#       filter(lag > 0)
#     jackfun.placebo = function(end.year) {
#       jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
#       reg.jack.placebo = lm(acee ~ ns(lag, df = 1),
#                             weights = wtd_pop_size,
#                             data = jack.data.placebo)
#       coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
#     }
#     full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
#     mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
#       mean(sapply(full.reg.placebo, function(v) v[i]))
#     })
#     cmp.ns.placebo.subcat = cmp.ns.placebo.subcat + 
#       geom_line(data = data.frame(Lag = all.lags, Diff = mean.reg.placebo),
#                 aes(x=Lag, y = Diff), 
#                 size = 0.4, 
#                 color = "grey")
#   }
# }

cmp.ns.placebo.housing = ggplot()
svi_housing_results = read_csv("./latest_missing_data_fixed/analysis/in_space_placebo/subcat_svi/housing_results.csv")
for (i in 1:50) {
  results_summary = svi_housing_results %>% 
    filter(iter == i) %>%
    select(-iter)
  acee_pe_placebo = results_summary %>%
    pivot_longer(cols = starts_with("acee_"), 
                 names_to = "outcome_year",
                 names_prefix = "acee_",
                 values_to = "acee") %>%
    mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
    filter(lag > 0)
  jackfun.placebo = function(end.year) {
    jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
    reg.jack.placebo = lm(acee ~ ns(lag, df = 1),
                          weights = wtd_pop_size,
                          data = jack.data.placebo)
    coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
  }
  full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
  mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
    mean(sapply(full.reg.placebo, function(v) v[i]))
  })
  cmp.ns.placebo.housing = cmp.ns.placebo.housing + 
    geom_line(data = data.frame(Lag = all.lags, Diff = mean.reg.placebo),
              aes(x=Lag, y = Diff), 
              size = 0.4, 
              color = "grey")
}

cmp.ns.placebo.housing = cmp.ns.placebo.housing +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.4),
            aes(x=Lag, y = Diff),
            color = "blue") +
  labs(title = "Placebo Test: Housing & Transportation")

ggsave("./latest_missing_data_fixed/analysis/in_space_placebo/housing.jpg")


#### socioeconomic domain
cmp.ns.placebo.ses = ggplot()
svi_ses_results = read_csv("./latest_missing_data_fixed/analysis/in_space_placebo/subcat_svi/ses_results.csv")
for (i in 1:50) {
  results_summary = svi_ses_results %>% 
    filter(iter == i) %>%
    select(-iter)
  acee_pe_placebo = results_summary %>%
    pivot_longer(cols = starts_with("acee_"), 
                 names_to = "outcome_year",
                 names_prefix = "acee_",
                 values_to = "acee") %>%
    mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
    filter(lag > 0)
  jackfun.placebo = function(end.year) {
    jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
    reg.jack.placebo = lm(acee ~ ns(lag, df = 1),
                          weights = wtd_pop_size,
                          data = jack.data.placebo)
    coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
  }
  full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
  mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
    mean(sapply(full.reg.placebo, function(v) v[i]))
  })
  cmp.ns.placebo.ses = cmp.ns.placebo.ses + 
    geom_line(data = data.frame(Lag = all.lags, Diff = mean.reg.placebo),
              aes(x=Lag, y = Diff), 
              size = 0.4, 
              color = "grey")
}

cmp.ns.placebo.ses = cmp.ns.placebo.ses +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.1),
            aes(x=Lag, y = Diff),
            color = "blue") +
  labs(title = "Placebo Test: Socioeconomic")

ggsave("./latest_missing_data_fixed/analysis/in_space_placebo/ses.jpg")

#### Household Decomposition domain
cmp.ns.placebo.hhd = ggplot()
svi_hhd_results = read_csv("./latest_missing_data_fixed/analysis/in_space_placebo/subcat_svi/hhd_results.csv")
for (i in 1:50) {
  results_summary = svi_hhd_results %>% 
    filter(iter == i) %>%
    select(-iter)
  acee_pe_placebo = results_summary %>%
    pivot_longer(cols = starts_with("acee_"), 
                 names_to = "outcome_year",
                 names_prefix = "acee_",
                 values_to = "acee") %>%
    mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
    filter(lag > 0)
  jackfun.placebo = function(end.year) {
    jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
    reg.jack.placebo = lm(acee ~ ns(lag, df = 1),
                          weights = wtd_pop_size,
                          data = jack.data.placebo)
    coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
  }
  full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
  mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
    mean(sapply(full.reg.placebo, function(v) v[i]))
  })
  cmp.ns.placebo.hhd = cmp.ns.placebo.hhd + 
    geom_line(data = data.frame(Lag = all.lags, Diff = mean.reg.placebo),
              aes(x=Lag, y = Diff), 
              size = 0.4, 
              color = "grey")
}

cmp.ns.placebo.hhd = cmp.ns.placebo.hhd +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.2),
            aes(x=Lag, y = Diff),
            color = "blue") +
  labs(title = "Placebo Test: Household Decomposition")

ggsave("./latest_missing_data_fixed/analysis/in_space_placebo/hhd.jpg")

cmp.ns.placebo.minority = ggplot()
svi_minority_results = read_csv("./latest_missing_data_fixed/analysis/in_space_placebo/subcat_svi/minority_results.csv")
all.lags.ns.3 = ns(all.lags, df = 2)
for (i in 1:50) {
  results_summary = svi_minority_results %>% 
    filter(iter == i) %>%
    select(-iter)
  acee_pe_placebo = results_summary %>%
    pivot_longer(cols = starts_with("acee_"), 
                 names_to = "outcome_year",
                 names_prefix = "acee_",
                 values_to = "acee") %>%
    mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
    filter(lag > 0)
  jackfun.placebo = function(end.year) {
    jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
    reg.jack.placebo = lm(acee ~ ns(lag, df = 2),
                          weights = wtd_pop_size,
                          data = jack.data.placebo)
    coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns.3[, 1] + 
      coef(reg.jack.placebo)[3] * all.lags.ns.3[, 2]
  }
  full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
  mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
    mean(sapply(full.reg.placebo, function(v) v[i]))
  })
  cmp.ns.placebo.minority = cmp.ns.placebo.minority + 
    geom_line(data = data.frame(Lag = all.lags, Diff = mean.reg.placebo),
              aes(x=Lag, y = Diff), 
              size = 0.4, 
              color = "grey")
}

cmp.ns.placebo.minority = cmp.ns.placebo.minority +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.3),
            aes(x=Lag, y = Diff),
            color = "blue") +
  labs(title = "Placebo Test: Minority Status / Language")
ggsave("./latest_missing_data_fixed/analysis/in_space_placebo/minority.jpg")


################################################################################
### calculate the p-values associated with placebo analysis

## household domain
svi_hhd_results = read_csv("./latest_missing_data_fixed/analysis/in_space_placebo/subcat_svi/hhd_results.csv")
pval_df = data.frame(iter = numeric(),
                     chng = numeric())
for (i in 1:50) {
  results_summary = svi_hhd_results %>% 
    filter(iter == i) %>%
    select(-iter)
  acee_pe_placebo = results_summary %>%
    pivot_longer(cols = starts_with("acee_"), 
                 names_to = "outcome_year",
                 names_prefix = "acee_",
                 values_to = "acee") %>%
    mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
    filter(lag > 0)
  jackfun.placebo = function(end.year) {
    jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
    reg.jack.placebo = lm(acee ~ ns(lag, df = 1),
                          weights = wtd_pop_size,
                          data = jack.data.placebo)
    coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
  }
  full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
  mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
    mean(sapply(full.reg.placebo, function(v) v[i]))
  })
  reg_m = lm(mean.reg.placebo ~ all.lags)
  slope = coef(reg_m)[2]
  pval_df = rbind(pval_df, data.frame(iter = i,
                                      chng = slope))
}
sum(abs(pval_df$chng) >= coef(lm(mean.reg.2 ~ all.lags))[2])
6 / 51


## minority status / english proficiency
svi_minority_results = read_csv("./latest_missing_data_fixed/analysis/in_space_placebo/subcat_svi/minority_results.csv")
all.lags.ns.3 = ns(all.lags, df = 2)
pval_df_minority = data.frame(iter = numeric(),
                              chng = numeric())
for (i in 1:50) {
  results_summary = svi_minority_results %>% 
    filter(iter == i) %>%
    select(-iter)
  acee_pe_placebo = results_summary %>%
    pivot_longer(cols = starts_with("acee_"), 
                 names_to = "outcome_year",
                 names_prefix = "acee_",
                 values_to = "acee") %>%
    mutate(lag = as.numeric(outcome_year) - as.numeric(exposure_year)) %>%
    filter(lag > 0)
  jackfun.placebo = function(end.year) {
    jack.data.placebo = acee_pe_placebo %>% filter(outcome_year == end.year) 
    reg.jack.placebo = lm(acee ~ ns(lag, df = 2),
                          weights = wtd_pop_size,
                          data = jack.data.placebo)
    coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns.3[, 1] + 
      coef(reg.jack.placebo)[3] * all.lags.ns.3[, 2]
  }
  full.reg.placebo = lapply(all.outcome.years, jackfun.placebo)
  mean.reg.placebo = sapply(seq_along(full.reg.placebo[[1]]), function(i) {
    mean(sapply(full.reg.placebo, function(v) v[i]))
  })
  reg_m = lm(mean.reg.placebo ~ all.lags)
  slope = coef(reg_m)[2]
  pval_df_minority = rbind(pval_df_minority, 
                           data.frame(iter = i,
                                      chng = slope))
}
sum(abs(pval_df_minority$chng) >= abs(coef(lm(mean.reg.3 ~ all.lags))[2]))
