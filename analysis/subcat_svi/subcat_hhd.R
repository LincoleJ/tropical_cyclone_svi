# calculate & plot acee for variables under "household decomposition" domain
print(Sys.time())
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(dplyr)
library(mgcv)
library(splines)
## for exposure years 2000, 2014, 2016, and 2018, the variables under this domain
## is consistent: percentage elderly, percentage youth, single parent household, 
## percentage disability
## for 2010, percentage disability is absent; for 2020, English proficiency is 
## added (from minority status domain)

####-----Calculate ACEE Point Estimates for subcategories of SVI-----####
for (T_0 in 2005:2018) {
  df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", T_0, ".rds")) %>%
    filter(fips != 42033,
           fips != 35039)
  res = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
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
  svi_grouped_age17 = read_csv("./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_age17.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_age17 = c("fips", paste0("svi_", colnames(svi_grouped_age17)[-1]))
  colnames(svi_grouped_age17) = svi_colnames_age17
  # keep only treated/control units
  svi_grouped_age17 = svi_grouped_age17[svi_grouped_age17$fips %in% c(control, exposed), ]
  svi_grouped_age17$exposed = 0
  svi_grouped_age17[svi_grouped_age17$fips %in% exposed, ]$exposed = 1
  svi_final_age17 = merge(df_final, svi_grouped_age17, by = "fips")
  # calculate summed weights
  svi_final_age17 = svi_final_age17 %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_age17 = svi_final_age17 %>%
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
  write.csv(acee_df_age17, paste0("./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_age17/acee_", T_0, ".csv"))
  
  svi_grouped_age65 = read_csv("./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_age65.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_age65 = c("fips", paste0("svi_", colnames(svi_grouped_age65)[-1]))
  colnames(svi_grouped_age65) = svi_colnames_age65
  # keep only treated/control units
  svi_grouped_age65 = svi_grouped_age65[svi_grouped_age65$fips %in% c(control, exposed), ]
  svi_grouped_age65$exposed = 0
  svi_grouped_age65[svi_grouped_age65$fips %in% exposed, ]$exposed = 1
  svi_final_age65 = merge(df_final, svi_grouped_age65, by = "fips")
  # calculate summed weights
  svi_final_age65 = svi_final_age65 %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_age65 = svi_final_age65 %>%
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
  write.csv(acee_df_age65, paste0("./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_age65/acee_", T_0, ".csv"))
  
  svi_grouped_disabl = read_csv("./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_disabl.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_disabl = c("fips", paste0("svi_", colnames(svi_grouped_disabl)[-1]))
  colnames(svi_grouped_disabl) = svi_colnames_disabl
  # keep only treated/control units
  svi_grouped_disabl = svi_grouped_disabl[svi_grouped_disabl$fips %in% c(control, exposed), ]
  svi_grouped_disabl$exposed = 0
  svi_grouped_disabl[svi_grouped_disabl$fips %in% exposed, ]$exposed = 1
  svi_final_disabl = merge(df_final, svi_grouped_disabl, by = "fips")
  # calculate summed weights
  svi_final_disabl = svi_final_disabl %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_disabl = svi_final_disabl %>%
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
  write.csv(acee_df_disabl, paste0("./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_disabl/acee_", T_0, ".csv"))
  
  svi_grouped_sngpnt = read_csv("./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_sngpnt.csv")[-c(1)] %>%
    filter(fips != 42033,
           fips != 35039)
  svi_colnames_sngpnt = c("fips", paste0("svi_", colnames(svi_grouped_sngpnt)[-1]))
  colnames(svi_grouped_sngpnt) = svi_colnames_sngpnt
  # keep only treated/control units
  svi_grouped_sngpnt = svi_grouped_sngpnt[svi_grouped_sngpnt$fips %in% c(control, exposed), ]
  svi_grouped_sngpnt$exposed = 0
  svi_grouped_sngpnt[svi_grouped_sngpnt$fips %in% exposed, ]$exposed = 1
  svi_final_sngpnt = merge(df_final, svi_grouped_sngpnt, by = "fips")
  # calculate summed weights
  svi_final_sngpnt = svi_final_sngpnt %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df_sngpnt = svi_final_sngpnt %>%
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
  write.csv(acee_df_sngpnt, paste0("./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_sngpnt/acee_", T_0, ".csv"))
}

# extract acee estimates
acee_pe_age17 = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_age17/",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_age65 = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_age65",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_disabl = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_disabl",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_sngpnt = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subvars_hhd/svi_sngpnt",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
wtd_pop_ctrl = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/preliminary/wtd_pop_ctrl/",
                                                pattern = "\\.csv$", 
                                                full.names = T), 
                                     read_csv))[, 2:5] %>%
  filter(lag == 0) %>%
  select(-outcome_year, -lag)
acee_pe_age17 = merge(acee_pe_age17, wtd_pop_ctrl, by = "exposure_year")
acee_pe_age65 = merge(acee_pe_age65, wtd_pop_ctrl, by = "exposure_year")
acee_pe_disabl = merge(acee_pe_disabl, wtd_pop_ctrl, by = "exposure_year")
acee_pe_sngpnt = merge(acee_pe_sngpnt, wtd_pop_ctrl, by = "exposure_year")


### evaluate optimal df for natural splines
jack.eval = function(df, jack.data) {
  reg.jack = lm(acee ~ ns(lag, df = df),
                weights = wtd_pop_ctrl,
                data = jack.data)
  (AIC(reg.jack) + BIC(reg.jack)) / 2
} 
jack.eval(1, acee_pe_age17)
jack.eval(2, acee_pe_age17)
jack.eval(3, acee_pe_age17)
jack.eval(4, acee_pe_age17)
jack.eval(5, acee_pe_age17)
# optimal df = 1 for pp age < 17
jack.eval(1, acee_pe_age65)
jack.eval(2, acee_pe_age65)
jack.eval(3, acee_pe_age65)
jack.eval(4, acee_pe_age65)
jack.eval(5, acee_pe_age65)
# optimal df = 1 for pp age > 65
jack.eval(1, acee_pe_disabl)
jack.eval(2, acee_pe_disabl)
jack.eval(3, acee_pe_disabl)
jack.eval(4, acee_pe_disabl)
jack.eval(5, acee_pe_disabl)
# optimal df = 1 for pp disabled
jack.eval(1, acee_pe_sngpnt)
jack.eval(2, acee_pe_sngpnt)
jack.eval(3, acee_pe_sngpnt)
jack.eval(4, acee_pe_sngpnt)
jack.eval(5, acee_pe_sngpnt)
# optimal df = 1 for pp single parent households

all.lags = 1:15
all.outcome.years = c(2010, 2014, 2016, 2018, 2020)
all.lags.ns = ns(all.lags, df = 1)

# for age less than 17
acee_pe_age17 = acee_pe_age17 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
jackfun.age17 = function(end.year) {
  jack.data = acee_pe_age17 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.age17 = lapply(all.outcome.years, jackfun.age17)
mean.reg.age17 = sapply(seq_along(full.reg.age17[[1]]), function(i) {
  mean(sapply(full.reg.age17, function(v) v[i]))
})
cmp.plot.age17 = ggplot(acee_pe_age17, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in Percentile") +
  xlab("Year since exposure") +
  ggtitle("Percentage Youngers (Age < 17)") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.age17),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

acee_pe_age65 = acee_pe_age65 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
jackfun.age65 = function(end.year) {
  jack.data = acee_pe_age65 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.age65 = lapply(all.outcome.years, jackfun.age65)
mean.reg.age65 = sapply(seq_along(full.reg.age65[[1]]), function(i) {
  mean(sapply(full.reg.age65, function(v) v[i]))
})
cmp.plot.age65 = ggplot(acee_pe_age65, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in Percentile") +
  xlab("Year since exposure") +
  ggtitle("Percentage Elderly (Age > 65)") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.age65),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

acee_pe_disabl = acee_pe_disabl %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
all.outcome.years.disabl = c(2014, 2016, 2018, 2020)
jackfun.disabl = function(end.year) {
  jack.data = acee_pe_age17 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.disabl = lapply(all.outcome.years.disabl, jackfun.disabl)
mean.reg.disabl = sapply(seq_along(full.reg.disabl[[1]]), function(i) {
  mean(sapply(full.reg.disabl, function(v) v[i]))
})
cmp.plot.disabl = ggplot(acee_pe_disabl, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in Percentile") +
  xlab("Year since exposure") +
  ggtitle("Percentage Disabled") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.disabl),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

acee_pe_sngpnt = acee_pe_sngpnt %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
jackfun.sngpnt = function(end.year) {
  jack.data = acee_pe_sngpnt %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg.sngpnt = lapply(all.outcome.years, jackfun.sngpnt)
mean.reg.sngpnt = sapply(seq_along(full.reg.sngpnt[[1]]), function(i) {
  mean(sapply(full.reg.sngpnt, function(v) v[i]))
})
cmp.plot.sngpnt = ggplot(acee_pe_sngpnt, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in Percentile") +
  xlab("Year since exposure") +
  ggtitle("Percentage Single-Parent Household") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.sngpnt),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)


##############------Confidence Intervals------#############
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
    ggtitle(paste0("On ", domain)) + 
    xlab("Years since exposure") +
    ylab("Difference in outcome") 
  ci.plot
}
jackreps.age17 = t(sapply(1:length(all.outcome.years), 
                          function(ii) {
                            full.reg.age17 = lapply(all.outcome.years[-ii], jackfun.age17)
                            sapply(seq_along(full.reg.age17[[1]]), function(i) {
                              mean(sapply(full.reg.age17, function(v) v[i]))
                            })}))
domain.age17 = c("Percentage Youngster (age < 17)")
jack.ci.plot(jackreps.age17, mean.reg.age17, domain.age17, all.lags)

jackreps.age65 = t(sapply(1:length(all.outcome.years), 
                          function(ii) {
                            full.reg.age65 = lapply(all.outcome.years[-ii], jackfun.age65)
                            sapply(seq_along(full.reg.age65[[1]]), function(i) {
                              mean(sapply(full.reg.age65, function(v) v[i]))
                            })}))
domain.age65 = c("Percentage Elderly (age > 65)")
jack.ci.plot(jackreps.age65, mean.reg.age65, domain.age65, all.lags)

jackreps.sngpnt = t(sapply(1:length(all.outcome.years), 
                          function(ii) {
                            full.reg.sngpnt = lapply(all.outcome.years[-ii], jackfun.sngpnt)
                            sapply(seq_along(full.reg.sngpnt[[1]]), function(i) {
                              mean(sapply(full.reg.sngpnt, function(v) v[i]))
                            })}))
domain.sngpnt = c("Percentage Single Parent Household")
jack.ci.plot(jackreps.sngpnt, mean.reg.sngpnt, domain.sngpnt, all.lags)

jackreps.disabl = t(sapply(1:length(all.outcome.years.disabl), 
                           function(ii) {
                             full.reg.disabl = lapply(all.outcome.years.disabl[-ii], jackfun.disabl)
                             sapply(seq_along(full.reg.disabl[[1]]), function(i) {
                               mean(sapply(full.reg.disabl, function(v) v[i]))
                             })}))
domain.disabl = c("Percentage Disabled")
jack.ci.plot(jackreps.disabl, mean.reg.disabl, domain.disabl, all.lags)

path = "./latest_missing_data_fixed/analysis/subcat_svi/plots"
ggsave(paste0(path, "/cmp.plot.age17.png"),
       cmp.plot.age17)
ggsave(paste0(path, "/cmp.plot.age65.png"),
       cmp.plot.age65)
ggsave(paste0(path, "/cmp.plot.disabl.png"),
       cmp.plot.disabl)
ggsave(paste0(path, "/cmp.plot.sngpnt.png"),
       cmp.plot.sngpnt)
ggsave(paste0(path, "/ci.plot.age17.png"),
       jack.ci.plot(jackreps.age17, mean.reg.age17, domain.age17, all.lags))
ggsave(paste0(path, "/ci.plot.age65.png"),
       jack.ci.plot(jackreps.age65, mean.reg.age65, domain.age65, all.lags))
ggsave(paste0(path, "/ci.plot.disabl.png"),
       jack.ci.plot(jackreps.disabl, mean.reg.disabl, domain.disabl, all.lags))
ggsave(paste0(path, "/ci.plot.sngpnt.png"),
       jack.ci.plot(jackreps.sngpnt, mean.reg.sngpnt, domain.sngpnt, all.lags))

