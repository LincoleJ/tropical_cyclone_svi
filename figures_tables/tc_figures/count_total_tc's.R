# count all tropical cyclones in contiguous u.s. states
rm(list=ls())
library(hurricaneexposure)
library(hurricaneexposuredata)
library(data.table)
library(dplyr)
library(tidyverse)

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

ct_names <- county_centers[,c(1,2)] # county fips code

# extract tropical cyclone events
tc_raw <- county_wind(counties = ct_names$fips,
                      start_year = 2005,
                      end_year = 2018,
                      wind_limit = 17.4911) %>%
  mutate(year = substr(closest_date, 0, 4),
         event_type = "tropical storm") %>%
  select(fips, storm_id, year, event_type)

# number of tropical cyclones
length(unique(tc_raw$storm_id))

# number of TCs that contributed to yearly exposures
tc_raw %>% group_by(year) %>%
  select(fips, storm_id) %>% View()

# number of (yearly) tropical cyclone exposures
nrow(tc_raw %>% distinct(year, fips))
nrow(tc_raw %>% distinct(fips))

# exposure count on the most exposed counties
n_exposed = tc_raw %>% group_by(fips) %>% 
  select(fips, year) %>% distinct() %>% summarise(exposures = n())
n_exposed[order(-n_exposed$exposures), ]
