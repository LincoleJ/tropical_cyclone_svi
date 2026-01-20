rm(list=ls())
library(dplyr)
library(tidyverse)
library(Hmisc)
library(splines)

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

# Population data
pop_grouped <- list.files(path = "./raw-data/census-data/cdc_population_monthly_infer/pct_hispanic",
                          pattern = "*.csv", full.name = TRUE) %>% 
  read_csv() %>% 
  filter(month == 6) %>% 
  select(-month) %>% 
  group_by(fips, year) %>%
  summarise(pop = sum(pop), .groups = "keep") %>%
  pivot_wider(names_from = year, values_from = pop, names_prefix = "pop_") %>%
  exclude_non_contiguous_states() %>%
  filter(fips != 42033, fips != 35039)

all.lags <- 1:15
all.outcome.years <- c(2010, 2014, 2016, 2018, 2020)
all.lags.ns <- ns(all.lags, df = 1)

########################################################################
# STEP 1: Calculate TRUE exposed effects for each household variable
########################################################################

hhd_vars <- c("age17", "age65", "disabl", "sngpnt")
true_effects <- list()

for (var_name in hhd_vars) {
  
  cat("Calculating true exposed effect for:", var_name, "\n")
  
  # Read SVI data for this variable
  svi_file_name <- paste0("./raw-data/svi-data/svi_hhd_var/svi_dat_", 
                          var_name, ".csv")
  svi_grouped <- read_csv(svi_file_name)[, -1] %>%
    filter(fips != 42033, fips != 35039)
  svi_colnames <- c("fips", paste0("svi_", colnames(svi_grouped)[-1]))
  colnames(svi_grouped) <- svi_colnames
  
  # Determine outcome years
  if (var_name == "disabl") {
    outcome_years <- c(2014, 2016, 2018, 2020)
  } else {
    outcome_years <- c(2010, 2014, 2016, 2018, 2020)
  }
  
  # Calculate ACEE for each exposure year using TRUE weights
  acee_true <- data.frame()
  
  for (T_0 in 2005:2018) {
    
    df <- readRDS(paste0("./processed-data/df_", T_0, ".rds")) %>%
      filter(fips != 35039, fips != 42033)
    
    res <- readRDS(list.files(file.path("./balancing/weights-data-upd"),
                              pattern = paste(T_0),
                              full.names = TRUE))
    
    df_weight <- df %>% drop_na()
    df_weight$weight <- df_weight$exposed * res$weights.1 + 
      (1 - df_weight$exposed) * res$weights.0
    df_weight <- df_weight %>% 
      mutate(weight = case_when(is.nan(weight) ~ 1, TRUE ~ weight))
    
    control <- df_weight$fips[df_weight$weight != 0 & df_weight$weight != 1]
    exposed <- df_weight$fips[df_weight$weight == 1]
    
    df_final <- df_weight[df_weight$fips %in% c(control, exposed), ] %>%
      select(fips, weight)
    
    svi <- svi_grouped[svi_grouped$fips %in% c(control, exposed), ]
    svi_final <- merge(df_final, svi, by = "fips")
    
    # Calculate ACEE for each outcome year
    for (yr in outcome_years) {
      col_name <- paste0("svi_", yr)
      if (col_name %in% names(svi_final)) {
        exposed_mean <- mean(subset(svi_final, weight == 1)[[col_name]], na.rm = TRUE)
        ctrl_data <- subset(svi_final, weight != 1)
        ctrl_weighted <- sum(ctrl_data$weight * ctrl_data[[col_name]], na.rm = TRUE) / 
          sum(ctrl_data$weight)
        
        # Get weighted pop size at exposure
        pop_size <- pop_grouped[pop_grouped$fips %in% c(control, exposed), ]
        pop_final <- merge(df_final, pop_size, by = "fips")
        pop_col <- paste0("pop_", T_0)
        ctrl_pop <- subset(pop_final, weight != 1)
        wtd_pop_size <- sum(ctrl_pop[[pop_col]] * ctrl_pop$weight, na.rm = TRUE) / 
          sum(ctrl_pop$weight)
        
        acee_true <- rbind(acee_true, data.frame(
          exposure_year = T_0,
          outcome_year = yr,
          lag = yr - T_0,
          acee = exposed_mean - ctrl_weighted,
          wtd_pop_ctrl = wtd_pop_size
        ))
      }
    }
  }
  
  acee_true <- acee_true %>% filter(lag > 0)
  
  # Fit model to get true effect curve
  jackfun_true <- function(end.year) {
    jack.data <- acee_true %>% filter(outcome_year == end.year)
    if (nrow(jack.data) < 2) return(rep(NA, length(all.lags)))
    reg.jack <- lm(acee ~ ns(lag, df = 1),
                   weights = wtd_pop_ctrl,
                   data = jack.data)
    coef(reg.jack)[1] + coef(reg.jack)[2] * all.lags.ns
  }
  
  full.reg.true <- lapply(outcome_years, jackfun_true)
  mean.reg.true <- sapply(seq_along(full.reg.true[[1]]), function(i) {
    mean(sapply(full.reg.true, function(v) v[i]), na.rm = TRUE)
  })
  
  true_effects[[var_name]] <- mean.reg.true
  
  # Also save the slope for p-value calculation
  true_slope <- coef(lm(mean.reg.true ~ all.lags))[2]
  attr(true_effects[[var_name]], "slope") <- true_slope
}

########################################################################
# STEP 2: Load placebo results and create plots with TRUE blue line
########################################################################

weights_data <- read_csv("./analysis/in_space_placebo/weights_data_(set_seed).csv")[, -1]
pvals = list()
for (var_name in hhd_vars) {
  
  cat("Creating placebo plot for:", var_name, "\n")
  
  # Load placebo results (you should have generated these already)
  placebo_results <- read_csv(paste0("./analysis/in_space_placebo/hhd_var/", var_name, "_results.csv"))[, -1]
  
  if (var_name == "disabl") {
    outcome_years <- c(2014, 2016, 2018, 2020)
  } else {
    outcome_years <- c(2010, 2014, 2016, 2018, 2020)
  }
  
  # Build plot with grey placebo lines
  cmp.ns <- ggplot()
  
  # Store slopes for p-value calculation
  placebo_slopes <- numeric(50)
  
  for (i in 1:50) {
    results_summary <- placebo_results %>% 
      filter(iter == i) %>%
      select(-iter)
    
    acee_pe_placebo <- results_summary %>%
      pivot_longer(cols = starts_with("acee_"), 
                   names_to = "outcome_year",
                   names_prefix = "acee_",
                   values_to = "acee") %>%
      mutate(outcome_year = as.numeric(outcome_year),
             lag = outcome_year - exposure_year) %>%
      filter(lag > 0, outcome_year %in% outcome_years) %>%
      drop_na(acee)
    
    jackfun.placebo <- function(end.year) {
      jack.data.placebo <- acee_pe_placebo %>% filter(outcome_year == end.year) 
      if (nrow(jack.data.placebo) < 2) return(rep(NA, length(all.lags)))
      reg.jack.placebo <- lm(acee ~ ns(lag, df = 1),
                             weights = wtd_pop_size,
                             data = jack.data.placebo)
      coef(reg.jack.placebo)[1] + coef(reg.jack.placebo)[2] * all.lags.ns
    }
    
    full.reg.placebo <- lapply(outcome_years, jackfun.placebo)
    mean.reg.placebo <- sapply(seq_along(full.reg.placebo[[1]]), function(j) {
      mean(sapply(full.reg.placebo, function(v) v[j]), na.rm = TRUE)
    })
    
    # Store slope for p-value
    placebo_slopes[i] <- coef(lm(mean.reg.placebo ~ all.lags))[2]
    
    # Add grey line
    cmp.ns <- cmp.ns + 
      geom_line(data = data.frame(Lag = all.lags, Diff = mean.reg.placebo),
                aes(x = Lag, y = Diff), 
                size = 0.4,
                color = "grey")
  }
  
  # Add TRUE exposed effect as blue line
  true_effect <- true_effects[[var_name]]
  true_slope <- attr(true_effects[[var_name]], "slope")
  
  final_plot <- cmp.ns + 
    geom_line(data = data.frame(Lag = all.lags, Diff = true_effect),
              aes(x = Lag, y = Diff), 
              color = "blue", size = 1.25) +
    labs(x = "Lag", y = "Difference in outcome",
         title = paste("Placebo Test:", var_name)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13, 15)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  ggsave(paste0("./analysis/in_space_placebo/hhd_var/placebo_", var_name, "_corrected.jpg"), 
         final_plot, width = 8, height = 6)
}