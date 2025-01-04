## create synthetic control regions using covariate balance weights
## calculate the relative difference of SVI between exposed region 
## and synthetic control region by year lag
## conduct sensitivity analysis by removing poorly balanced exposure years
print(Sys.time())
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(dplyr)
library(mgcv)
library(splines)

####-----Calculate ACEE Point Estimates and weighted population size-----####
for (T_0 in 2005:2018) {
  df = readRDS(paste0("./processed-data/df_", T_0, ".rds")) %>% 
    filter(fips != 35039,
           fips != 42033)
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
  svi_grouped = read_csv("./raw-data/svi-data/svi_dat.csv")[-c(1)] %>%
    filter(fips != 35039,
           fips != 42033)
  svi_colnames = c("fips", paste0("svi_", colnames(svi_grouped)[-1]))
  colnames(svi_grouped) = svi_colnames
  # keep only treated/control units
  svi_grouped = svi_grouped[svi_grouped$fips %in% c(control, exposed), ]
  svi_grouped$exposed = 0
  svi_grouped[svi_grouped$fips %in% exposed, ]$exposed = 1
  svi_final = merge(df_final, svi_grouped, by = "fips")
  # calculate summed weights
  svi_final = svi_final %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup()
  acee_df = svi_final %>%
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
  write.csv(acee_df, paste0("./analysis/preliminary/acee/acee_", T_0, ".csv"))
}

for (T_0 in 2005:2018) {
  # get county-level population
  df = readRDS(paste0("./processed-data/df_", T_0, ".rds")) %>%
    filter(fips != 35039,
           fips != 42033)
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
  
  # get county-level population at time of exposure
  pop = list.files(path = "./raw-data/census-data/cdc_population_monthly_infer/pct_hispanic",
                   pattern = "*.csv", full.name = TRUE) %>% 
    read_csv() %>% 
    filter(month == 6) %>% 
    select(-month) %>% 
    group_by(fips, year) %>%
    summarise(pop = sum(pop), .groups = "keep") %>%
    pivot_wider(names_from = year,
                values_from = pop,
                names_prefix = "pop_") %>%
    filter(fips != 35039,
           fips != 42033)
  pop_grouped = pop[pop$fips %in% c(control, exposed), ]
  pop_grouped$exposed = 0
  pop_grouped[pop_grouped$fips %in% exposed, ]$exposed = 1
  pop_final = merge(df_final, pop_grouped, by = "fips")
  # calculate summed weights
  pop_final = pop_final %>% group_by(exposed) %>%
    mutate(sum_weights = sum(weight)) %>%
    ungroup() 
  wtd_pop_df = pop_final %>%
    # for each unit, calculate weighted outcome as pop * weight / sum_weights
    mutate(across(starts_with("pop"), ~ . * weight / sum_weights)) %>%
    group_by(exposed) %>%
    # weighted outcome
    summarise(across(starts_with("pop"), sum)) %>%
    pivot_longer(cols = starts_with("pop"), names_to = "variable", values_to = "mean_value") %>%
    pivot_wider(names_from = exposed, values_from = mean_value) %>% 
    rename("exposed" = `1`,
           "control" = `0`) %>%
    mutate(exposure_year = T_0,
           outcome_year = as.numeric(gsub("[^0-9]", "", variable)),
           lag = outcome_year - exposure_year,
           wtd_pop_ctrl = control) %>%
    select(exposure_year, outcome_year, lag, wtd_pop_ctrl)
  write.csv(wtd_pop_df, paste0("./analysis/preliminary/wtd_pop_ctrl/wtd_pop_ctrl_", T_0, ".csv"))
}

#####---------Relative Difference in SVI vs Year Since Exposure---------#######
# extract acee estimates
acee_pe = do.call(rbind, lapply(list.files(path = "./analysis/preliminary/acee",
                                           pattern = "\\.csv$", 
                                           full.names = T), 
                                read_csv))[, 6:10] %>%
  filter(lag > 0)
wtd_pop_ctrl = do.call(rbind, lapply(list.files(path = "./analysis/preliminary/wtd_pop_ctrl/",
                                                pattern = "\\.csv$", 
                                                full.names = T), 
                                     read_csv))[, 2:5] %>%
  filter(lag == 0) %>%
  select(-outcome_year, -lag)
acee_pe = merge(acee_pe, wtd_pop_ctrl, by = "exposure_year")
# synthesized analysis using linear modeling
all.lags = 1:15
all.outcome.years = c(2010, 2014, 2016, 2018, 2020)
all.lags.ns = ns(all.lags, df = 1)

# obtain the optimal df for ns regression model
# optimal ns = 1 / linear model
jack.eval = function(df) {
  jack.data = acee_pe 
  reg.jack = lm(acee ~ ns(lag, df = df),
                weights = wtd_pop_ctrl,
                data = jack.data)
  (AIC(reg.jack) + BIC(reg.jack)) / 2
}

jackfun = function(end.year) {
  jack.data = acee_pe %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
}
full.reg = lapply(all.outcome.years, jackfun)
mean.reg = sapply(seq_along(full.reg[[1]]), function(i) {
  mean(sapply(full.reg, function(v) v[i]))
})

acee_pe = acee_pe %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))

cmp.plot = ggplot(acee_pe, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2)

#####---------C.I. for ACEE w/ jackknife--------#####
jackreps = t(sapply(1:length(all.outcome.years), 
                    function(ii) { 
                      full.reg = lapply(all.outcome.years[-ii], jackfun)
                      sapply(seq_along(full.reg[[1]]), function(i) {
                        mean(sapply(full.reg, function(v) v[i]))
                      })}))
colnames(jackreps) = all.lags
jackvar = apply(jackreps, 2, function(xx) { var(xx) * (length(xx) - 1)^2 / length(xx) })
jackse = sqrt(jackvar)
diff = mean.reg
ub.diff = mean.reg + 1.96 * jackse
lb.diff = mean.reg - 1.96 * jackse
results <-  data.frame(year = 1:length(diff), diff = diff, lower = lb.diff, upper = ub.diff)
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
  ggtitle(paste0("Effect of Tropical Cyclones on SVI" )) + 
  xlab("Years since exposure") +
  ylab("Difference in outcome") 
