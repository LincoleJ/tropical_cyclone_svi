# this file puts the C.I. figures together
print(Sys.time())
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(dplyr)
library(mgcv)
library(splines)
library(cowplot)

# extract acee estimates
acee_pe = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/preliminary/acee",
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
acee_pe = merge(acee_pe, wtd_pop_ctrl, by = "exposure_year")
all.lags = 1:15
all.outcome.years = c(2010, 2014, 2016, 2018, 2020)
all.lags.ns = ns(all.lags, df = 1)
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
  geom_ribbon(data = results, aes(ymin = lower, ymax = upper), fill = "grey", alpha=0.3) +
  geom_line(aes(x= year, y= diff), lwd=1.2, color = "lightcoral") +
  geom_line(aes(x= year, y= lower), lwd=0.3, color = "lightcoral") +
  geom_line(aes(x= year, y= upper), lwd=0.3, color = "lightcoral") +
  geom_line(aes(x= year, y= lower), linetype="dashed" , lwd=1, color = "lightcoral") +
  geom_line(aes(x= year, y= upper), linetype="dashed" , lwd=1, color = "lightcoral") +
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks=1:15) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "none",
        text = element_text(size=18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  ggtitle(paste0("Effects on Overall SVI" )) + 
  xlab("Years since exposure") +
  ylab("Weighted Mean Difference") 


#######################----- subcategories of SVI -----#######################
acee_pe_1 = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subcategories_svi/svi_ses",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_2 = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subcategories_svi/svi_hhd",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_3 = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subcategories_svi/svi_minority",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_4 = do.call(rbind, lapply(list.files(path = "./latest_missing_data_fixed/analysis/subcat_svi/subcategories_svi/svi_housing",
                                             pattern = "\\.csv$", 
                                             full.names = T), 
                                  read_csv))[, 6:10] %>%
  filter(lag > 0)
acee_pe_1 = merge(acee_pe_1, wtd_pop_ctrl, by = "exposure_year")
acee_pe_2 = merge(acee_pe_2, wtd_pop_ctrl, by = "exposure_year")
acee_pe_3 = merge(acee_pe_3, wtd_pop_ctrl, by = "exposure_year")
acee_pe_4 = merge(acee_pe_4, wtd_pop_ctrl, by = "exposure_year")

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

# For Household Composition Domain
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

# jackknife
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
jackreps.3 = t(sapply(1:length(all.outcome.years), 
                      function(ii) {
                        full.reg.3 = lapply(all.outcome.years[-ii], jackfun.3)
                        sapply(seq_along(full.reg.3[[1]]), function(i) {
                          mean(sapply(full.reg.3, function(v) v[i]))
                        })}))
jackreps.4 = t(sapply(1:length(all.outcome.years), 
                      function(ii) {
                        full.reg.4 = lapply(all.outcome.years[-ii], jackfun.4)
                        sapply(seq_along(full.reg.4[[1]]), function(i) {
                          mean(sapply(full.reg.4, function(v) v[i]))
                        })}))

domain.1 = c("Socioeconomic")
domain.2 = c("Household")
domain.3 = c("Minority Status")
domain.4 = c("Housing Type / Transportation")
jack.ci.plot = function(jackreps, mean.reg, domain, all.lags) {
  colnames(jackreps) = all.lags
  jackvar = apply(jackreps, 2, function(xx) { var(xx) * (length(xx) - 1)^2 / length(xx) })
  jackse = sqrt(jackvar)
  diff = mean.reg
  ub.diff = mean.reg + 1.96 * jackse
  lb.diff = mean.reg - 1.96 * jackse
  results <-  data.frame(year = 1:length(diff), diff = diff, lower = lb.diff, upper = ub.diff)
  ci.plot = ggplot(data = results, aes(x = year, y = diff)) +
    geom_ribbon(data = results, aes(ymin = lower, ymax = upper), fill = "grey", alpha=0.3) +
    geom_line(aes(x= year, y= diff), lwd=1.2, color = "lightcoral") +
    geom_line(aes(x= year, y= lower), lwd=0.3, color = "lightcoral") +
    geom_line(aes(x= year, y= upper), lwd=0.3, color = "lightcoral") +
    geom_line(aes(x= year, y= lower), linetype="dashed" , lwd=1, color = "lightcoral") +
    geom_line(aes(x= year, y= upper), linetype="dashed" , lwd=1, color = "lightcoral") +
    geom_hline(yintercept=0) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10, 12, 14, 16)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24),
          legend.position = "none",
          text = element_text(size=18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)) +
    ggtitle(paste0(domain, " SVI")) + 
    xlab("Years since exposure") +
    ylab("Weighted Mean Difference") 
  ci.plot
}
ci.plot.socio = jack.ci.plot(jackreps.1, mean.reg.1, domain.1, all.lags)
ci.plot.hhd = jack.ci.plot(jackreps.2, mean.reg.2, domain.2, all.lags)
ci.plot.minority = jack.ci.plot(jackreps.3, mean.reg.3, domain.3, all.lags)
ci.plot.housing = jack.ci.plot(jackreps.4, mean.reg.4, domain.4, all.lags)

# put plots together in one single plot
ept_p = ggplot() + theme_void()
plot.sum = plot_grid(
  plot_grid(ept_p, ci.plot, ept_p, 
            ncol = 3, 
            rel_widths = c(1, 3, 1), 
            rel_heights = c(1, 1, 1)),
  plot_grid(ci.plot.socio, ci.plot.hhd, ncol = 2),
  plot_grid(ci.plot.minority, ci.plot.housing, ncol = 2),
  nrow = 3,
  rel_heights = c(1.2, 1, 1)
)
ggsave("./latest_missing_data_fixed/figures_tables/ci_plots_conglomerate/ci_sum.jpg",
       plot.sum, width = 11.8, height = 12.8)
