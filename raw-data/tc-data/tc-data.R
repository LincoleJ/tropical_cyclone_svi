# install hurricaneexposure packages
library(drat)
addRepo("geanders")
install.packages("hurricaneexposuredata")

rm(list=ls())
library(hurricaneexposure)
library(hurricaneexposuredata)
library(data.table)
library(dplyr)
library(tidyverse)

# fix counties to be consistent across time, mutually exclusive, and collectively exhaustive
fix_county_codes = function(dat){
  dat = dat %>%
    mutate(fips = case_when(
      # Alaska FIPS codes to do here when have time
      fips== 'XXXXX' ~ 'XXXXX',
      fips== '08001' | fips== '08013' | fips== '08059' | fips== '08123' ~ '08014', # 08001, 08013, 08059, 08123 -> 08014
      fips== '12025' ~ '12086', #  12025 -> 12086
      fips== '30031' | fips== '30067'~ '30113', # 30113 -> 30031, 30067
      fips== '46113' ~ '46102', # 46113 -> 46102
      fips== '51560' ~ '51005', # 51560 -> 51005
      fips== '51780' ~ '51083', # 51780 -> 51083
      fips== '51515' ~ '51019', # 51515 -> 51019 
      TRUE ~ fips
    ))
}

# `fips` variable codes the five-digit code unique for each county.
# `year` variable codes the year of record from 2000 to 2019.
# `exposure` variable is coded as a binary variable:
  # 1 ~ there is a hurricanes (> 64 knots) at that exposed year
  # 0 ~ no tropical cyclones at all (no wind > 33 knots) at that exposed year.
  # 0.5 ~ there is at least one occurrence of tropical storm (>=33 but < 64, 
  # but no hurricane) at that exposed year.

ct_names <- county_centers[,c(1,2)] # county fips code

# extract hurricane events
# >= 64 in knots <--> >= 32.9244 in m/s
hurr_raw <- county_wind(counties = ct_names$fips,
                        start_year = 1995,
                        end_year = 2018,
                        wind_limit = 32.9244) %>%
  mutate(year = substr(closest_date, 0, 4),
         event_type = "hurricane") %>%
  select(fips, storm_id, year, event_type)

# extract tropical storm events
# >= 34 & <= 64 in knots <--> >= 17.4911 and <= 32.9244 in m/s
ts_raw <- county_wind(counties = ct_names$fips,
                      start_year = 1995,
                      end_year = 2018,
                      wind_limit = 17.4911) %>%
  filter(vmax_sust <= 32.9244) %>%
  mutate(year = substr(closest_date, 0, 4),
         event_type = "tropical storm") %>%
  select(fips, storm_id, year, event_type)

tc_dat <- rbind(hurr_raw, ts_raw) %>% 
  select(-storm_id) 

# check if the observations add up
# nrow(tc_dat) == nrow(ts_raw) + nrow(hurr_raw)

hurr_dat <- full_join(ct_names, tc_dat) %>%
  arrange(year) %>%
  mutate(exposure = case_when(event_type == "hurricane" ~ 1,
                              event_type == "tropical storm" ~ 0)) %>%
  select(-event_type) %>%
  pivot_wider(names_from = year,
              values_from = exposure,
              # value >= 1 ~ hurricane, == 0 ~ tropical storm, NA ~ no TC at all
              values_fn = sum) %>%
  select(-"NA") %>%
  pivot_longer(cols = "1995":"2018",
               names_to = "year",
               values_to = "exposure") %>%
  mutate(exposure = case_when(exposure >= 1 ~ 1,
                              exposure == 0 ~ 0.5,
                              is.na(exposure) ~ 0)) %>%
  mutate(exposure = case_when(exposure == 0 ~ 0,
                              exposure == 0.5 ~ 1,
                              exposure == 1 ~ 2)) %>%
  mutate(exposure = case_when(exposure == 0 ~ 0,
                              exposure == 1 | exposure == 2 ~ 1)) %>%
  pivot_wider(names_from = year,
              values_from = exposure)

# fix fips codes
hurr_dat_1 = hurr_dat %>% 
  select(-county_name) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(across(2:24, max)) %>%
  distinct()

write.csv(hurr_dat_1, file = "./latest_missing_data_fixed/raw-data/tc-data/hurr_dat.csv")

missing_entries <- hurr_dat %>% filter_all(any_vars(is.na(.)))
View(missing_entries)
# no missing data