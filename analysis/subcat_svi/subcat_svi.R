# conduct weighted outcome analysis on subcategories of SVI
print(Sys.time())
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(dplyr)
library(mgcv)
library(splines)

####-----Calculate ACEE Point Estimates for subcategories of SVI-----####
for (T_0 in 2005:2018) {
  df = readRDS(paste0("./processed-data/df_", T_0, ".rds")) %>%
    filter(fips != 42033,
           fips != 35039)
  res = readRDS(list.files(file.path("./balancing/weights-data-upd"),
                           pattern = paste(T_0),
                           full.names = TRUE))
  df_weight = df %>% drop_na()
  df_weight$weight = df_weight$exposed *res$weights.1 + 
    (1-df_weight$exposed)*res$weights.0
  df_weight = df_weight %>% mutate(weight = case_when(is.nan(weight) ~ 1,
                                                      TRUE ~ weight)) # fix infinite weights to be 1
  control = df_weight$fips[df_weight$weight != 0 & df_weight$weight != 1]
  exposed = df_weight$fips[df_weight$weight == 1]
  df_final = df_weight[df_weight$fips %in% c(control, exposed), ] %>%
    select(fips, weight)
  
  # import outcome variables (i.e., SVI from 2010, 2014, 2016, 2018, 2020)
  svi_grouped_1 = read_csv("./raw-data/svi-data/svi_dat_ses.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_1 = c("fips", paste0("svi_", colnames(svi_grouped_1)[-1]))
  colnames(svi_grouped_1) = svi_colnames_1
  # keep only treated/control units
  svi_grouped_1 = svi_grouped_1[svi_grouped_1$fips %in% c(control, exposed), ]
  svi_grouped_1$exposed = 0
  svi_grouped_1[svi_grouped_1$fips %in% exposed, ]$exposed = 1
  svi_final_1 = merge(df_final, svi_grouped_1, by = "fips")
  # calculate summed weights
  svi_final_1 = svi_final_1 %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_1 = svi_final_1 %>%
    # for each unit, calculate weighted outcome as svi * weight / sum_weights
    mutate(across(starts_with("svi"), ~ . * weight / sum_weights)) %>%
    group_by(exposed) %>%
    # weighted outcome
    summarise(across(starts_with("svi"), sum)) %>%
    pivot_longer(cols = starts_with("svi"), names_to = "variable", values_to = "mean_value") %>%
    pivot_wider(names_from = exposed, values_from = mean_value) %>%
    mutate(difference = `1` - `0`) %>% 
    rename("exposed" = `1`,
           "control" = `0`) %>% 
    mutate(exposure_year = T_0,
           outcome_year = as.numeric(gsub("[^0-9]", "", variable)),
           lag = outcome_year - exposure_year,
           acee = difference,
           baseline = control)
  write.csv(acee_df_1, paste0("./analysis/subcat_svi/subcategories_svi/svi_ses/acee_", T_0, ".csv"))
  
  # need to streamline here
  svi_grouped_2 = read_csv("./raw-data/svi-data/svi_dat_hhd.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_2 = c("fips", paste0("svi_", colnames(svi_grouped_2)[-1]))
  colnames(svi_grouped_2) = svi_colnames_2
  # keep only treated/control units
  svi_grouped_2 = svi_grouped_2[svi_grouped_2$fips %in% c(control, exposed), ]
  svi_grouped_2$exposed = 0
  svi_grouped_2[svi_grouped_2$fips %in% exposed, ]$exposed = 1
  svi_final_2 = merge(df_final, svi_grouped_2, by = "fips")
  # calculate summed weights
  svi_final_2 = svi_final_2 %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_2 = svi_final_2 %>%
    # for each unit, calculate weighted outcome as svi * weight / sum_weights
    mutate(across(starts_with("svi"), ~ . * weight / sum_weights)) %>%
    group_by(exposed) %>%
    # weighted outcome
    summarise(across(starts_with("svi"), sum)) %>%
    pivot_longer(cols = starts_with("svi"), names_to = "variable", values_to = "mean_value") %>%
    pivot_wider(names_from = exposed, values_from = mean_value) %>%
    mutate(difference = `1` - `0`) %>% 
    rename("exposed" = `1`,
           "control" = `0`) %>% 
    mutate(exposure_year = T_0,
           outcome_year = as.numeric(gsub("[^0-9]", "", variable)),
           lag = outcome_year - exposure_year,
           acee = difference,
           baseline = control)
  write.csv(acee_df_2, paste0("./analysis/subcat_svi/subcategories_svi/svi_hhd/acee_", T_0, ".csv"))
  
  svi_grouped_3 = read_csv("./raw-data/svi-data/svi_dat_minority.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_3 = c("fips", paste0("svi_", colnames(svi_grouped_3)[-1]))
  colnames(svi_grouped_3) = svi_colnames_3
  # keep only treated/control units
  svi_grouped_3 = svi_grouped_3[svi_grouped_3$fips %in% c(control, exposed), ]
  svi_grouped_3$exposed = 0
  svi_grouped_3[svi_grouped_3$fips %in% exposed, ]$exposed = 1
  svi_final_3 = merge(df_final, svi_grouped_3, by = "fips")
  # calculate summed weights
  svi_final_3 = svi_final_3 %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_3 = svi_final_3 %>%
    # for each unit, calculate weighted outcome as svi * weight / sum_weights
    mutate(across(starts_with("svi"), ~ . * weight / sum_weights)) %>%
    group_by(exposed) %>%
    # weighted outcome
    summarise(across(starts_with("svi"), sum)) %>%
    pivot_longer(cols = starts_with("svi"), names_to = "variable", values_to = "mean_value") %>%
    pivot_wider(names_from = exposed, values_from = mean_value) %>%
    mutate(difference = `1` - `0`) %>% 
    rename("exposed" = `1`,
           "control" = `0`) %>% 
    mutate(exposure_year = T_0,
           outcome_year = as.numeric(gsub("[^0-9]", "", variable)),
           lag = outcome_year - exposure_year,
           acee = difference,
           baseline = control)
  write.csv(acee_df_3, paste0("./analysis/subcat_svi/subcategories_svi/svi_minority/acee_", T_0, ".csv"))
  
  svi_grouped_4 = read_csv("./raw-data/svi-data/svi_dat_housing.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_4 = c("fips", paste0("svi_", colnames(svi_grouped_4)[-1]))
  colnames(svi_grouped_4) = svi_colnames_4
  # keep only treated/control units
  svi_grouped_4 = svi_grouped_4[svi_grouped_4$fips %in% c(control, exposed), ]
  svi_grouped_4$exposed = 0
  svi_grouped_4[svi_grouped_4$fips %in% exposed, ]$exposed = 1
  svi_final_4 = merge(df_final, svi_grouped_4, by = "fips")
  # calculate summed weights
  svi_final_4 = svi_final_4 %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_4 = svi_final_4 %>%
    # for each unit, calculate weighted outcome as svi * weight / sum_weights
    mutate(across(starts_with("svi"), ~ . * weight / sum_weights)) %>%
    group_by(exposed) %>%
    # weighted outcome
    summarise(across(starts_with("svi"), sum)) %>%
    pivot_longer(cols = starts_with("svi"), names_to = "variable", values_to = "mean_value") %>%
    pivot_wider(names_from = exposed, values_from = mean_value) %>%
    mutate(difference = `1` - `0`) %>% 
    rename("exposed" = `1`,
           "control" = `0`) %>% 
    mutate(exposure_year = T_0,
           outcome_year = as.numeric(gsub("[^0-9]", "", variable)),
           lag = outcome_year - exposure_year,
           acee = difference,
           baseline = control)
  write.csv(acee_df_4, paste0("./analysis/subcat_svi/subcategories_svi/svi_housing/acee_", T_0, ".csv"))
}


#####---------Relative Difference in SVI vs Year Since Exposure---------#######
# extract acee estimates
acee_pe_1 = do.call(rbind, lapply(list.files(path = "./analysis/subcat_svi/subcategories_svi/svi_ses",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_2 = do.call(rbind, lapply(list.files(path = "./analysis/subcat_svi/subcategories_svi/svi_hhd",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_3 = do.call(rbind, lapply(list.files(path = "./analysis/subcat_svi/subcategories_svi/svi_minority",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_4 = do.call(rbind, lapply(list.files(path = "./analysis/subcat_svi/subcategories_svi/svi_housing",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_2.5 = acee_pe_2
acee_pe_2.5$acee = (acee_pe_2$acee + acee_pe_3$acee) / 2
wtd_pop_ctrl = do.call(rbind, lapply(list.files(path = "./analysis/preliminary/wtd_pop_ctrl/",
                                                pattern = "\\.csv$", 
                                                full.names = T), 
                                     read_csv))[, 2:5] %>%
  filter(lag == 0) %>%
  select(-outcome_year, -lag)
acee_pe_1 = merge(acee_pe_1, wtd_pop_ctrl, by = "exposure_year")
acee_pe_2 = merge(acee_pe_2, wtd_pop_ctrl, by = "exposure_year")
acee_pe_3 = merge(acee_pe_3, wtd_pop_ctrl, by = "exposure_year")
acee_pe_4 = merge(acee_pe_4, wtd_pop_ctrl, by = "exposure_year")
acee_pe_2.5 = merge(acee_pe_2.5, wtd_pop_ctrl, by = "exposure_year")
all.lags = 1:15
all.outcome.years = c(2010, 2014, 2016, 2018, 2020)
all.lags.ns = ns(all.lags, df = 1)
## for socioeconomic domain
jack.eval = function(df, jack.data) {
  reg.jack = lm(acee ~ ns(lag, df = df),
                weights = wtd_pop_ctrl,
                data = jack.data)
  (AIC(reg.jack) + BIC(reg.jack)) / 2
} 
jack.eval(1, acee_pe_1)
jack.eval(2, acee_pe_1)
jack.eval(3, acee_pe_1)
jack.eval(4, acee_pe_1)
jack.eval(5, acee_pe_1)
# for socioeconomic domain, optimal df = 1
jack.eval(1, acee_pe_2)
jack.eval(2, acee_pe_2)
jack.eval(3, acee_pe_2)
jack.eval(4, acee_pe_2)
jack.eval(5, acee_pe_2)
# for household domain, optimal df = 1
jack.eval(1, acee_pe_3)
jack.eval(2, acee_pe_3)
jack.eval(3, acee_pe_3)
jack.eval(4, acee_pe_3)
jack.eval(5, acee_pe_3)
# for minority domain, optimal df = 2
jack.eval(1, acee_pe_4)
jack.eval(2, acee_pe_4)
jack.eval(3, acee_pe_4)
jack.eval(4, acee_pe_4)
# for housing domain, optimal df = 1
jack.eval(1, acee_pe_2.5)
jack.eval(2, acee_pe_2.5)
jack.eval(3, acee_pe_2.5)
jack.eval(4, acee_pe_2.5)
jack.eval(5, acee_pe_2.5)
#先用df=1 

# for socioeconomic domain
acee_pe_1 = acee_pe_1 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
jackfun.1 = function(end.year) {
  jack.data = acee_pe_1 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.1 = lapply(all.outcome.years, jackfun.1)
mean.reg.1 = sapply(seq_along(full.reg.1[[1]]), function(i) {
  mean(sapply(full.reg.1, function(v) v[i]))
})
cmp.plot.1 = ggplot(acee_pe_1, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  ggtitle("Socioeconomic Domain") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.1),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

# for household domain
acee_pe_2 = acee_pe_2 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
jackfun.2 = function(end.year) {
  jack.data = acee_pe_2 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.2 = lapply(all.outcome.years, jackfun.2)
mean.reg.2 = sapply(seq_along(full.reg.2[[1]]), function(i) {
  mean(sapply(full.reg.2, function(v) v[i]))
})
cmp.plot.2 = ggplot(acee_pe_2, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  ggtitle("Household Composition Domain") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.2),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

# for minority status / language domain
acee_pe_3 = acee_pe_3 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
all.lags.ns.3 = ns(all.lags, df = 2)
jackfun.3 = function(end.year) {
  jack.data = acee_pe_3 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 2),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns.3[, 1] + 
    coef(reg.jack)[3] * all.lags.ns.3[, 2]
}
full.reg.3 = lapply(all.outcome.years, jackfun.3)
mean.reg.3 = sapply(seq_along(full.reg.3[[1]]), function(i) {
  mean(sapply(full.reg.3, function(v) v[i]))
})
cmp.plot.3 = ggplot(acee_pe_3, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  ggtitle("Minority Status / Language Domain") +
  scale_y_continuous(breaks=seq(0, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.3),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

# for housing / transportation domain
acee_pe_4 = acee_pe_4 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
jackfun.4 = function(end.year) {
  jack.data = acee_pe_4 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.4 = lapply(all.outcome.years, jackfun.4)
mean.reg.4 = sapply(seq_along(full.reg.4[[1]]), function(i) {
  mean(sapply(full.reg.4, function(v) v[i]))
})
cmp.plot.4 = ggplot(acee_pe_4, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  ggtitle("Housing / Transportation Domain") +
  scale_y_continuous(breaks=seq(0, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.4),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

# combined 
acee_pe_2.5 = acee_pe_2.5 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
jackfun.2.5 = function(end.year) {
  jack.data = acee_pe_2.5 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.2.5 = lapply(all.outcome.years, jackfun.2.5)
mean.reg.2.5 = sapply(seq_along(full.reg.2.5[[1]]), function(i) {
  mean(sapply(full.reg.2.5, function(v) v[i]))
})
cmp.plot.2.5 = ggplot(acee_pe_2.5, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  ggtitle("Household Composition & Minority & Language") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.2.5),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

# excluding 2020 
all.lags.b = 1:13
all.outcome.years.b = c(2010, 2014, 2016, 2018)
# all.lags.ns = ns(all.lags, df = 1)
# for household characteristics domain
acee_pe_2b = acee_pe_2 %>% filter(outcome_year != 2020)
jack.eval(1, acee_pe_2b)
jack.eval(2, acee_pe_2b)
jack.eval(3, acee_pe_2b)
jack.eval(4, acee_pe_2b)
jack.eval(5, acee_pe_2b)
# optimal df = 2
all.lags.ns.2b = ns(all.lags.b, df = 2)
jack.fun.2b = function(end.year) {
  jack.data = acee_pe_2b %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 2),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns.2b[, 1] + 
    coef(reg.jack)[3] * all.lags.ns.2b[, 2]
}
full.reg.2b = lapply(all.outcome.years.b, jack.fun.2b)
mean.reg.2b = sapply(seq_along(full.reg.2b[[1]]), function(i) {
  mean(sapply(full.reg.2b, function(v) v[i]))
})
cmp.plot.2b = ggplot(acee_pe_2b, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  ggtitle("Household Composition Domain") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags.b, Diff = mean.reg.2b),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

# for minority status & language proficiency
# excluding 2020
acee_pe_3b = acee_pe_3 %>% filter(outcome_year != 2020)
jack.eval(1, acee_pe_3b)
jack.eval(2, acee_pe_3b)
jack.eval(3, acee_pe_3b)
jack.eval(4, acee_pe_3b)
jack.eval(5, acee_pe_3b)
# optimal df = 2
all.lags.ns.3b = ns(all.lags.b, df = 2)
jack.fun.3b = function(end.year) {
  jack.data = acee_pe_3b %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 2),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns.3b[, 1] + 
    coef(reg.jack)[3] * all.lags.ns.3b[, 2]
}
full.reg.3b = lapply(all.outcome.years.b, jack.fun.3b)
mean.reg.3b = sapply(seq_along(full.reg.3b[[1]]), function(i) {
  mean(sapply(full.reg.3b, function(v) v[i]))
})
acee_pe_3b$outcome_year <- as.factor(acee_pe_3b$outcome_year)
cmp.plot.3b = ggplot(acee_pe_3b, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  ggtitle("Minority Status & Language (2020 Excluded)") +
  scale_y_continuous(breaks=seq(-1, 1, 0.01)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags.b, Diff = mean.reg.3b),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)


#########-------Confidence Intervals----------##########
jackreps.1 = t(sapply(1:length(all.outcome.years), 
                      function(ii) {
                        full.reg.1 = lapply(all.outcome.years[-ii], jackfun.1)
                        sapply(seq_along(full.reg.1[[1]]), function(i) {
                          mean(sapply(full.reg.1, function(v) v[i]))
                        })}))
jackreps.2 = t(sapply(1:length(all.outcome.years), 
                      function(ii) {
                        full.reg.2 = lapply(all.outcome.years[-ii], jackfun.2)
                        sapply(seq_along(full.reg.2[[1]]), function(i) {
                          mean(sapply(full.reg.2, function(v) v[i]))
                        })}))
jackreps.2b = t(sapply(1:length(all.outcome.years.b), 
                      function(ii) {
                        full.reg.2b = lapply(all.outcome.years.b[-ii], jack.fun.2b)
                        sapply(seq_along(full.reg.2b[[1]]), function(i) {
                          mean(sapply(full.reg.2b, function(v) v[i]))
                        })}))
jackreps.3 = t(sapply(1:length(all.outcome.years), 
                      function(ii) {
                        full.reg.3 = lapply(all.outcome.years[-ii], jackfun.3)
                        sapply(seq_along(full.reg.3[[1]]), function(i) {
                          mean(sapply(full.reg.3, function(v) v[i]))
                        })}))
jackreps.3b = t(sapply(1:length(all.outcome.years.b), 
                       function(ii) {
                         full.reg.3b = lapply(all.outcome.years.b[-ii], jack.fun.3b)
                         sapply(seq_along(full.reg.3b[[1]]), function(i) {
                           mean(sapply(full.reg.3b, function(v) v[i]))
                         })}))
jackreps.4 = t(sapply(1:length(all.outcome.years), 
                      function(ii) {
                        full.reg.4 = lapply(all.outcome.years[-ii], jackfun.4)
                        sapply(seq_along(full.reg.4[[1]]), function(i) {
                          mean(sapply(full.reg.4, function(v) v[i]))
                        })}))
jackreps.2.5 = t(sapply(1:length(all.outcome.years), 
                      function(ii) {
                        full.reg.2.5 = lapply(all.outcome.years[-ii], jackfun.2.5)
                        sapply(seq_along(full.reg.2.5[[1]]), function(i) {
                          mean(sapply(full.reg.2.5, function(v) v[i]))
                        })}))
domain.1 = c("Socioeconomic")
domain.2 = c("Household")
domain.2b = c("Household, excluded 2020")
domain.3 = c("Minority Status")
domain.4 = c("Household Type / Transportation")
domain.2.5 = c("Household, Minority, Language")
domain.3b = c("Minority & Language, excluded 2020")
jack.ci.plot = function(jackreps, mean.reg, domain, all.lags) {
  colnames(jackreps) = all.lags
  jackvar = apply(jackreps, 2, function(xx) { var(xx) * (length(xx) - 1)^2 / length(xx) })
  jackse = sqrt(jackvar)
  diff = mean.reg
  ub.diff = mean.reg + 1.96 * jackse
  lb.diff = mean.reg - 1.96 * jackse
  results <-  data.frame(year = 1:length(diff), diff = diff, lower = lb.diff, upper = ub.diff)
  # domain = ifelse(names(jackreps))
  ci.plot = ggplot(data = results, aes(x = year, y = diff)) +
    geom_ribbon(data = results, aes(ymin = lower, ymax = upper), fill = "grey70", alpha=0.3) +
    geom_line(aes(x= year, y= diff), lwd=1.2) +
    geom_line(aes(x= year, y= lower), linetype="dashed" , lwd=1.2) +
    geom_line(aes(x= year, y= upper), linetype="dashed" , lwd=1.2) +
    geom_hline(yintercept=0) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10, 12, 14, 16)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24),
          legend.position = "none",
          text = element_text(size=24),
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24)) +
    ggtitle(paste0("Effect of Tropical Cyclones on ", domain, " SVI")) + 
    xlab("Years since exposure") +
    ylab("Difference in outcome") 
  ci.plot
}
jack.ci.plot(jackreps.1, mean.reg.1, domain.1, all.lags)
jack.ci.plot(jackreps.2, mean.reg.2, domain.2, all.lags)
jack.ci.plot(jackreps.3, mean.reg.3, domain.3, all.lags)
jack.ci.plot(jackreps.4, mean.reg.4, domain.4, all.lags)
jack.ci.plot(jackreps.2.5, mean.reg.2.5, domain.2.5, all.lags)
jack.ci.plot(jackreps.2b, mean.reg.2b, domain.2b, all.lags.b)
jack.ci.plot(jackreps.3b, mean.reg.3b, domain.3b, all.lags.b)
# same images
path = "./analysis/subcat_svi/plots"
ggsave(paste0(path, "/cmp.plot.1.png"),
       cmp.plot.1)
ggsave(paste0(path, "/cmp.plot.2.png"),
       cmp.plot.2)
ggsave(paste0(path, "/cmp.plot.3.png"),
       cmp.plot.3)
ggsave(paste0(path, "/cmp.plot.4.png"),
       cmp.plot.4)
ggsave(paste0(path, "/cmp.plot.2.5.png"),
       cmp.plot.2.5)
ggsave(paste0(path, "/ci.plot.1.png"),
       jack.ci.plot(jackreps.1, mean.reg.1, domain.1, all.lags))
ggsave(paste0(path, "/ci.plot.2.png"),
       jack.ci.plot(jackreps.2, mean.reg.2, domain.2, all.lags))
ggsave(paste0(path, "/ci.plot.3.png"),
       jack.ci.plot(jackreps.3, mean.reg.3, domain.3, all.lags))
ggsave(paste0(path, "/ci.plot.4.png"),
       jack.ci.plot(jackreps.4, mean.reg.4, domain.4, all.lags))
ggsave(paste0(path, "/ci.plot.2.5.png"),
       jack.ci.plot(jackreps.2.5, mean.reg.2.5, domain.2.5, all.lags))

ggsave(paste0(path, "/cmp.plot.2b.png"),
       cmp.plot.2b)
ggsave(paste0(path, "/cmp.plot.3b.png"),
       cmp.plot.3b)
ggsave(paste0(path, "/ci.plot.2b.png"),
       jack.ci.plot(jackreps.2b, mean.reg.2b, domain.2b, all.lags.b))
ggsave(paste0(path, "/ci.plot.3b.png"),
       jack.ci.plot(jackreps.3b, mean.reg.3b, domain.3b, all.lags.b))
