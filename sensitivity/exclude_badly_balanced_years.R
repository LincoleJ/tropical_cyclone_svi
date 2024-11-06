library(tidyverse)
library(ggplot2)
library(dplyr)
library(mgcv)
library(splines)

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

acee_pe_sens_1 = acee_pe %>% filter(exposure_year != 2009,
                                    exposure_year != 2014,
                                    exposure_year != 2016)
# find optimal df for outcome model
jack.eval.sens.1 = function(df) {
  jack.data = acee_pe_sens_1
  reg.jack = lm(acee ~ ns(lag, df = df),
                weights = wtd_pop_ctrl,
                data = jack.data)
  (AIC(reg.jack) + BIC(reg.jack)) / 2
}
jack.eval.sens.1(1)
jack.eval.sens.1(2)
jack.eval.sens.1(3)
jack.eval.sens.1(4)
# df = 1 is optimal
all.lags = 1:15
all.outcome.years = c(2010, 2014, 2016, 2018, 2020)
all.lags.ns.sens.1 = ns(all.lags, df = 1)
jackfun_sens_1 = function(end.year) {
  jack.data = acee_pe_sens_1 %>% filter(outcome_year == end.year) 
  reg.jack = lm(acee ~ ns(lag, df = 1),
                weights = wtd_pop_ctrl,
                data = jack.data)
  coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns.sens.1
}
full.reg.sens.1 = lapply(all.outcome.years, jackfun_sens_1)
mean.reg.sens.1 = sapply(seq_along(full.reg.sens.1[[1]]), function(i) {
  mean(sapply(full.reg.sens.1, function(v) v[i]))
})
acee_pe_sens_1 = acee_pe_sens_1 %>% arrange(outcome_year) %>%
  mutate(outcome_year = as.character(outcome_year))
cmp.plot.sens.1 = ggplot(acee_pe_sens_1, aes(x = lag, y = acee)) +
  geom_point(aes(colour = outcome_year, size = wtd_pop_ctrl)) +
  scale_colour_manual(name = "Outcome year",
                      values=c("#FFD300", "#FF00B6", "#009FFF", "#783FC1", "#00FFBE")) +
  theme_grey(base_size = 14) +
  ylab("Relative Difference in SVI") +
  xlab("Year since exposure") +
  scale_y_continuous(breaks=seq(-1, 1, 0.05)) +
  scale_x_continuous(breaks=1:15) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(data=data.frame(Lag = all.lags, Diff = mean.reg.sens.1),
            aes(x=Lag, y = Diff), size = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.title = element_text(size = 20),
        text = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  geom_line(aes(y = 0), lty = 2) + 
  ggtitle("Well-Balanced Years Only")

# confidence interval
jackreps.sens.1 = t(sapply(1:length(all.outcome.years), 
                           function(ii) { 
                             full.reg.sens.1 = lapply(all.outcome.years[-ii], jackfun_sens_1)
                             sapply(seq_along(full.reg.sens.1[[1]]), function(i) {
                               mean(sapply(full.reg.sens.1, function(v) v[i]))
                             })}))
colnames(jackreps.sens.1) = all.lags
jackvar.sens.1 = apply(jackreps.sens.1, 2, function(xx) { var(xx) * (length(xx) - 1)^2 / length(xx) })
jackse.sens.1 = sqrt(jackvar.sens.1)
diff.sens.1 = mean.reg.sens.1
ub.diff.sens.1 = mean.reg.sens.1 + 1.96 * jackse.sens.1
lb.diff.sens.1 = mean.reg.sens.1 - 1.96 * jackse.sens.1
results.sens.1 <-  data.frame(year = 1:length(diff.sens.1), diff = diff.sens.1, 
                              lower = lb.diff.sens.1, upper = ub.diff.sens.1)
ci.plot.sens.1 = ggplot(data = results.sens.1, aes(x = year, y = diff)) +
  geom_ribbon(data = results.sens.1, aes(ymin = lower, ymax = upper), fill = "grey70", alpha=0.3) +
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
  ggtitle(paste0("Effect of Tropical Cyclones on SVI")) + 
  xlab("Years since exposure") +
  ylab("Difference in outcome") 

ggsave("./latest_missing_data_fixed/figures_tables/ci_exclu_poor_bal.jpg",
       width = 14, height = 8)
