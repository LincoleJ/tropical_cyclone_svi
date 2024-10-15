# the trajectory of SVI of synthetic control units and exposed units in 2005
rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)

# get
T_0 = 2005
df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", T_0, ".rds")) %>% 
  filter(fips != 35039,
         fips != 42033)
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
  pivot_longer(cols = starts_with("svi"), names_to = "year", values_to = "svi") %>%
  mutate(exposed = case_when(exposed == 0 ~ "synth",
                             exposed == 1 ~ "exposed")) %>%
  mutate(year = as.numeric(gsub("svi_", "", year)))

interpolated_acee = data.frame(exposed = c("synth", "exposed", "synth", "exposed",
                                           "synth", "exposed", "synth", "exposed"),
                               year = c(2005, 2005, 
                                        1995, 1995), 
                                        # 2006, 2006, 
                                        # 2007, 2007),
                               svi = c(0.7795626, 0.7795626 + 0.0034081, 
                                       0.76, 0.76 + 0.0034081)) 
                                       # 0.00195626 + 0.7795626, 0.8079707, 
                                       # 0.7834751, 0.8029707))
interpolated_acee = data.frame(exposed = c("synth", "exposed", "synth", "exposed",
                                           "synth", "exposed"),
                               year = c(2005, 2005,
                                        2000, 2000,
                                        1995, 1995), 
                               # 2006, 2006, 
                               # 2007, 2007),
                               svi = c(0.7795626, 0.7795626 + 0.001, 
                                       0.7697813, 0.7697813 + 0.001,
                                       0.76, 0.76 + 0.001)) 
acee_df_in = rbind(acee_df, interpolated_acee)

after_exposure = acee_df_in %>% filter(year >= 2005)
before_exposure = acee_df_in %>% filter(year <= 2005)
scm_evolve = ggplot() +
  geom_smooth(data = acee_df_in, 
              aes(x = year, y = svi, group = exposed, color = exposed),
              se = FALSE) +
  scale_color_manual(values = c("synth" = "green", "exposed" = "red")) +
  scale_x_continuous(breaks = c(1999, 2005, 2013),  # Define breaks
                     labels = c("Pre-Exposed", "Exposed", "Post-Exposed")) +
  scale_y_continuous(limits = c(0.7, 0.85)) +
  theme_bw() +
  geom_vline(xintercept=2005.2, color = "darkgrey", 
             linetype = "dashed") +
  annotate("text", x = 2013, y = 0.8, label = "Exposed Unit", size = 5, color = "green") +
  annotate("text", x = 2013, y = 0.73, label = "Synthetic Control Unit", size = 5, color = "red") +
  geom_vline(xintercept=2004.7, color = "khaki2", linewidth = 2.5) +
  geom_vline(xintercept=2005.7, color = "khaki2", linewidth = 2.5) +
  labs(x = NULL, y = "SVI") +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) 

ggsave("./latest_missing_data_fixed/figures_tables/scm_overview/scm_overview.jpg", scm_evolve)
