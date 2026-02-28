# Katrina Case Study: Aggregate Synthetic Control
# Donor pool: Counties not exposed to any TC in 2005

rm(list = ls())
library(tidyverse)
library(dplyr)
library(Hmisc)
library(hurricaneexposure)
library(hurricaneexposuredata)
source("./balancing/cbps_ATT.R")

# Helper functions
exclude_non_contiguous_states <- function(df) {
  contiguous_states <- c("01", "04", "05", "06", "08", "09", "10", "11", "12", 
                         "13", "16", "17", "18", "19", "20", "21", "22", "23", 
                         "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", 
                         "34", "35", "36", "37", "38", "39", "40", "41", "42", "44", 
                         "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", 
                         "56")
  contiguous_rows <- df[substr(df$fips, 1, 2) %in% contiguous_states, ]
  return(contiguous_rows)
}

fix_county_codes <- function(dat) {
  dat <- dat %>%
    mutate(fips = case_when(
      fips == 'XXXXX' ~ 'XXXXX',
      fips == '08001' | fips == '08013' | fips == '08059' | fips == '08123' ~ '08014',
      fips == '12025' ~ '12086',
      fips == '30031' | fips == '30067' ~ '30113',
      fips == '46113' ~ '46102',
      fips == '51560' ~ '51005',
      fips == '51780' ~ '51083',
      fips == '51515' ~ '51019',
      TRUE ~ fips
    ))
  return(dat)
}

dir.create("./sensitivity/sandy_case_study", recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------------------------
# Step 1: Get Sandy-specific exposure
#------------------------------------------------------------------------------

ct_names <- county_centers[, c(1, 2)]

sandy_exposure <- county_wind(
  counties = ct_names$fips,
  start_year = 2012,
  end_year = 2012,
  wind_limit = 17.4911
) %>%
  filter(storm_id == "Sandy-2012") %>%
  fix_county_codes()

sandy_fips <- unique(sandy_exposure$fips)
cat("Counties exposed to Sandy:", length(sandy_fips), "\n")

#------------------------------------------------------------------------------
# Step 2: Get all counties exposed to ANY TC in 2012 (to exclude from controls)
#------------------------------------------------------------------------------

all_2012_exposure <- county_wind(
  counties = ct_names$fips,
  start_year = 2012,
  end_year = 2012,
  wind_limit = 17.4911
) %>%
  select(fips, storm_id) %>%
  fix_county_codes()

all_2012_exposed_fips <- unique(all_2012_exposure$fips)
cat("Counties exposed to any TC in 2012:", length(all_2012_exposed_fips), "\n")

#------------------------------------------------------------------------------
# Step 3: Load 2005 data, define Katrina exposure and control pool
#------------------------------------------------------------------------------

df_full <- readRDS("./processed-data/df_2012.rds") %>%
  filter(fips != 42033, fips != 35039)

# Control pool: not exposed to ANY TC in 2012
control_fips <- df_full$fips[!df_full$fips %in% all_2012_exposed_fips]
cat("Control counties (no 2012 TC exposure):", length(control_fips), "\n")

# Keep Katrina-exposed + unexposed-in-2012 controls
df <- df_full %>%
  filter(fips %in% sandy_fips | fips %in% control_fips)

df$exposed <- ifelse(df$fips %in% sandy_fips, 1, 0)

cat("Sandy exposed in analysis:", sum(df$exposed), "\n")
cat("Controls in analysis:", sum(df$exposed == 0), "\n")

#------------------------------------------------------------------------------
# Step 4: Run CBPS
#------------------------------------------------------------------------------

X <- df %>% drop_na()
W <- X$exposed
fips_vec <- X$fips
X$exposed <- NULL
X$fips <- NULL

X.mean <- colMeans(X)
X.sd <- apply(X, 2, sd)
X.sd[X.sd == 0] <- 1
X.scl <- scale(X, center = X.mean, scale = X.sd)

res_regu.list <- lapply(1:8, function(n) {
  res <- cbps_att(as.matrix(X.scl),
                  W,
                  theta.init = rep(0, ncol(X) + 1),
                  control = list(trace = 10, maxit = 5000),
                  lambda = rep(10^{n - 7}, ncol(X)))
  return(res)
})

converge_set <- sapply(res_regu.list, function(res) res$convergence)
res <- res_regu.list[[min(which(converge_set == 0))]]
lambda <- 10^{min(which(converge_set == 0)) - 7}

cat("Converged with lambda:", lambda, "\n")

saveRDS(res, "./sensitivity/sandy_case_study/weights_sandy.RDS")

#------------------------------------------------------------------------------
# Step 5: Calculate aggregate ACEE
#------------------------------------------------------------------------------

df_weight <- df %>% drop_na()
df_weight$fips <- fips_vec
df_weight$weight <- df_weight$exposed * res$weights.1 + 
  (1 - df_weight$exposed) * res$weights.0
df_weight <- df_weight %>%
  mutate(weight = case_when(is.nan(weight) ~ 1, TRUE ~ weight))

# Load SVI
svi_dat <- read_csv("./raw-data/svi-data/svi_dat.csv") %>%
  select(-1) %>%
  filter(fips != 42033, fips != 35039) %>%
  exclude_non_contiguous_states()

colnames(svi_dat) <- c("fips", paste0("svi_", colnames(svi_dat)[-1]))

# Merge
df_final <- df_weight %>% select(fips, exposed, weight)
svi_final <- merge(df_final, svi_dat, by = "fips")

# Aggregate ACEE for each outcome year
outcome_years <- c("2000", "2010", "2014", "2016", "2018", "2020")

acee_results <- data.frame()

for (outcome_year in outcome_years) {
  svi_col <- paste0("svi_", outcome_year)
  
  # Weighted mean for exposed
  exposed_mean <- svi_final %>%
    filter(exposed == 1) %>%
    summarise(mean_svi = sum(weight * !!sym(svi_col), na.rm = TRUE) / sum(weight)) %>%
    pull(mean_svi)
  
  # Weighted mean for synthetic control
  control_mean <- svi_final %>%
    filter(exposed == 0, weight > 0) %>%
    summarise(mean_svi = sum(weight * !!sym(svi_col), na.rm = TRUE) / sum(weight)) %>%
    pull(mean_svi)
  
  acee <- exposed_mean - control_mean
  lag <- as.numeric(outcome_year) - 2012
  
  acee_results <- rbind(acee_results, data.frame(
    outcome_year = as.numeric(outcome_year),
    exposure_year = 2012,
    lag = lag,
    exposed_svi = exposed_mean,
    synthetic_control_svi = control_mean,
    acee = acee
  ))
}

write_csv(acee_results, "./sensitivity/sandy_case_study/sandy_acee.csv")

cat("\nSandy ACEE Results:\n")
print(acee_results)

#------------------------------------------------------------------------------
# Step 6: Plot
#------------------------------------------------------------------------------

plot_data <- acee_results %>%
  select(outcome_year, exposed_svi, synthetic_control_svi) %>%
  pivot_longer(cols = c(exposed_svi, synthetic_control_svi),
               names_to = "series", values_to = "svi") %>%
  mutate(series = ifelse(series == "exposed_svi", "Exposed", "Synthetic Control"))

p <- ggplot(plot_data, aes(x = outcome_year, y = svi, 
                           linetype = series, color = series)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 2012, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Exposed" = "red", "Synthetic Control" = "blue")) +
  scale_linetype_manual(values = c("Exposed" = "solid", "Synthetic Control" = "dashed")) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2014, 2016, 2018, 2020)) +
  labs(x = "Year", y = "SVI", color = NULL, linetype = NULL,
       title = "Hurricane Sandy (2012): SVI Trajectories",
       subtitle = "Donor pool: Counties with no TC exposure in 2012") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("./sensitivity/sandy_case_study/sandy_trajectory.png", 
       p, width = 8, height = 5, dpi = 150)


p_main <- ggplot(plot_data, aes(x = outcome_year, y = svi, 
                                linetype = series, color = series)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 2012, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Exposed" = "red", "Synthetic Control" = "blue")) +
  scale_linetype_manual(values = c("Exposed" = "solid", "Synthetic Control" = "dashed")) +
  scale_x_continuous(breaks = c(2000, 2010, 2014, 2016, 2018, 2020)) +
  labs(x = NULL, y = "SVI", color = NULL, linetype = NULL,
       title = "Hurricane Sandy (2012): SVI Trajectories",
       subtitle = "Donor pool: Counties with no TC exposure in 2012") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Calculate differences
diff_data <- plot_data %>%
  pivot_wider(names_from = series, values_from = svi) %>%
  mutate(difference = Exposed - `Synthetic Control`)

# Difference plot
p_diff <- ggplot(diff_data, aes(x = outcome_year, y = difference)) +
  geom_line(linewidth = 1, color = "black") +
  geom_point(size = 2.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 2012, linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = c(2000, 2010, 2014, 2016, 2018, 2020)) +
  labs(x = "Year", y = "ATT") +
  theme_bw() +
  theme(plot.margin = margin(0, 5.5, 5.5, 5.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Stack them with patchwork
p_combined <- p_main / p_diff + 
  plot_layout(heights = c(5, 1))

ggsave("./sensitivity/sandy_case_study/sandy_trajectory_with_att.png", 
       p_combined, width = 8, height = 6, dpi = 600)

#------------------------------------------------------------------------------
# Step 7: Check covariate balance
#------------------------------------------------------------------------------

covariates <- c("exposure_", "svi_", "tmean_", "ppt_mean_", "pct_hsg_", 
                "pct_black_", "pct_hispanic_", "pci_",
                "pct_pov_", "mf_ratio_", "pct_youth_", "pct_elderly_", 
                "pop_density_", "pop_size_")

df_cb <- df %>% drop_na()
df_cb$fips <- fips_vec
df_cb$weight <- df_cb$exposed * res$weights.1 + 
  (1 - df_cb$exposed) * res$weights.0
df_cb <- df_cb %>%
  mutate(weight = case_when(is.nan(weight) ~ 1, TRUE ~ weight))

post_bal <- NULL
pre_bal <- NULL
ncov <- length(covariates)

for (i in 1:ncov) {
  covariate <- covariates[i]
  
  post_bal[i] <- (mean(rowMeans(subset(df_cb, exposed == 1)[, grep(covariate, colnames(df))])) - 
                    sum(rowMeans(subset(df_cb, exposed == 0)[, grep(covariate, colnames(df))]) * 
                          subset(df_cb, exposed == 0)$weight) / sum(subset(df_cb, exposed == 0)$weight)) /
    sqrt(wtd.var(rowMeans(subset(df_cb, exposed == 0)[, grep(covariate, colnames(df))]), 
                 subset(df_cb, exposed == 0)$weight))
  
  pre_bal[i] <- (mean(rowMeans(subset(df_cb, exposed == 1)[, grep(covariate, colnames(df))])) - 
                   mean(rowMeans(subset(df_cb, exposed == 0)[, grep(covariate, colnames(df))]))) /
    sqrt(wtd.var(rowMeans(subset(df_cb, exposed == 0)[, grep(covariate, colnames(df))]), 
                 subset(df_cb, exposed == 0)$weight))
}

post_bal <- replace(post_bal, is.infinite(post_bal), 0)
pre_bal <- replace(pre_bal, is.infinite(pre_bal), 0)

balance_df <- data.frame(
  covariate = covariates,
  pre_balance_smd = pre_bal,
  post_balance_smd = post_bal
)

cat("\nCovariate Balance:\n")
print(balance_df)
cat("\nMean ASMD pre-balance:", mean(abs(pre_bal)), "\n")
cat("Mean ASMD post-balance:", mean(abs(post_bal)), "\n")

write_csv(balance_df, "./analysis/katrina_case_study/katrina_covariate_balance.csv")

cat("\nDone!\n")