print(Sys.time())
rm(list = ls())
library("sf")
library("tidyverse")
library("mltools")
library("data.table")
library("fst")
library("readr")
library("tidyverse")
library("purrr")

# read census data
census_df = readRDS("./latest_missing_data_fixed/raw-data/census-data/census_df.rds")
census_postfixes = as.numeric(gsub(".*_([0-9]+)$", "\\1", colnames(census_df[, -1])))

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

#--- 2010-2018 hurricane data ---#
for (y in 2005:2018) {
  
  # import hurricane exposure data
  hurr_dat <- read_csv("./latest_missing_data_fixed/raw-data/tc-data/hurr_dat.csv")[-1] %>% #select(-county_name) %>%
    pivot_longer(cols = !fips, names_to = "year", values_to = "exposure") %>% filter(year <= y) %>%
    janitor::clean_names() %>%
    pivot_wider(names_from = year, values_from = exposure) %>%
    rename_with(.fn = function(.x){paste0("exposure_", .x)}, .cols = -c(fips))
  
  # import previous svi data
  svi_dat <- read_csv("./latest_missing_data_fixed/raw-data/svi-data/svi_dat.csv")[-1] %>% 
    pivot_longer(cols = !fips, names_to = "year", values_to = "svi") %>%
    janitor::clean_names() %>% filter(year < y) %>%
    pivot_wider(names_from = year, values_from = svi) %>%
    rename_with(.fn = function(.x){paste0("svi_", .x)}, .cols = -c(fips))
  
  # import meteorological data
  weather_dat <- read_csv("./latest_missing_data_fixed/raw-data/weather-data/seasonal_means.csv")[-1] %>% 
    filter(year < y) %>% janitor::clean_names() %>%
    pivot_wider(names_from = year, values_from = c(tmean, ppt_mean))
  
  # import census data (> 1 year prior)
  selected_indices = which(census_postfixes < y) + 1
  all_indices = c(1, selected_indices)
  census_dat = census_df[, all_indices]
  
  # combine datasets
  list_df <- list(hurr_dat, svi_dat, weather_dat, census_dat) 
  treated = paste0("exposure_", y)
  df <- list_df %>% reduce(full_join, by = "fips") %>% 
    mutate(exposed = case_when(.data[[treated]] == 1 | .data[[treated]] == 2 ~ 1,
                               .data[[treated]] == 0 ~ 0)) %>%
    select(-.data[[treated]])
  
  # exclude all counties not in contiguous U.S. states
  df = exclude_non_contiguous_states(df)
  
  # save as RDS file
  saveRDS(df, paste0("./latest_missing_data_fixed/processed-data/df_", y,".rds"))
}

# check missing data
missing_id = apply(df_2018, 1, function(x) any(is.na(x)))
missing_entries = df_2018[missing_id, ]