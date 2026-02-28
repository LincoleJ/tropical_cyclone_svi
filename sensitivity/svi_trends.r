# Descriptive analysis: Raw SVI trends + Synthetic Control, aligned by lag
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(dplyr)

# Import SVI data
svi_grouped = read_csv("./raw-data/svi-data/svi_dat.csv")[-c(1)] %>%
  filter(fips != 35039, fips != 42033)
svi_colnames = c("fips", paste0("svi_", colnames(svi_grouped)[-1]))
colnames(svi_grouped) = svi_colnames

# Calculate raw SVI means for each exposure year cohort
raw_svi_list <- list()

for (T_0 in 2005:2018) {
  df = readRDS(paste0("./processed-data/df_", T_0, ".rds")) %>% 
    filter(fips != 35039, fips != 42033) %>%
    drop_na()
  
  exposed_fips = df$fips[df$exposed == 1]
  control_fips = df$fips[df$exposed == 0]
  
  svi_exposed = svi_grouped %>% filter(fips %in% exposed_fips)
  svi_control = svi_grouped %>% filter(fips %in% control_fips)
  
  raw_svi <- data.frame(
    exposure_year = T_0,
    outcome_year = c(2000, 2010, 2014, 2016, 2018, 2020),
    mean_exposed = c(mean(svi_exposed$svi_2000, na.rm = TRUE),
                     mean(svi_exposed$svi_2010, na.rm = TRUE),
                     mean(svi_exposed$svi_2014, na.rm = TRUE),
                     mean(svi_exposed$svi_2016, na.rm = TRUE),
                     mean(svi_exposed$svi_2018, na.rm = TRUE),
                     mean(svi_exposed$svi_2020, na.rm = TRUE)),
    mean_control = c(mean(svi_control$svi_2000, na.rm = TRUE),
                     mean(svi_control$svi_2010, na.rm = TRUE),
                     mean(svi_control$svi_2014, na.rm = TRUE),
                     mean(svi_control$svi_2016, na.rm = TRUE),
                     mean(svi_control$svi_2018, na.rm = TRUE),
                     mean(svi_control$svi_2020, na.rm = TRUE))
  )
  
  raw_svi_list[[as.character(T_0)]] <- raw_svi
}

raw_svi_df <- bind_rows(raw_svi_list)

# Load synthetic control data from ACEE files
acee_files <- list.files(path = "./analysis/preliminary/acee",
                         pattern = "\\.csv$", 
                         full.names = TRUE)
acee_df <- do.call(rbind, lapply(acee_files, read_csv)) %>%
  select(exposure_year, outcome_year, exposed, control) %>%
  rename(sc_exposed = exposed,
         sc_control = control)

# Merge raw and synthetic control
combined_df <- raw_svi_df %>%
  left_join(acee_df, by = c("exposure_year", "outcome_year")) %>%
  mutate(lag = outcome_year - exposure_year)

# Reshape for plotting
combined_long <- combined_df %>%
  pivot_longer(cols = c(mean_exposed, mean_control, sc_control),
               names_to = "group",
               values_to = "svi") %>%
  mutate(group = case_when(
    group == "mean_exposed" ~ "Exposed",
    group == "mean_control" ~ "Control (Raw)",
    group == "sc_control" ~ "Control (Synthetic)"
  )) %>%
  drop_na(svi)

# Plot: 7 rows x 2 columns, x-axis = lag, aligned at 0
svi_plot <- ggplot(combined_long, aes(x = lag, y = svi, color = group, linetype = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~exposure_year, nrow = 7, ncol = 2) +
  scale_color_manual(values = c("Exposed" = "red", 
                                "Control (Raw)" = "blue", 
                                "Control (Synthetic)" = "darkgreen")) +
  scale_linetype_manual(values = c("Exposed" = "solid",
                                   "Control (Raw)" = "dashed",
                                   "Control (Synthetic)" = "solid")) +
  scale_x_continuous(breaks = seq(-10, 15, by = 5)) +
  labs(x = "Years Since Exposure", 
       y = "Mean SVI",
       color = NULL,
       linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("./sensitivity/raw_vs_sc_svi_aligned.jpg", svi_plot,
       width = 6, height = 12)
