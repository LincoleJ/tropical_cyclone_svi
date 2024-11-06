rm(list=ls())
library(dplyr)
library(tidyverse)

# test code
readRDS("./latest_missing_data_fixed/raw-data/weather-data/tmeans_raw_data/weighted_area_raster_fips_tmean_daily_1995.rds")

########################--------helper functions--------########################
tmeans <- function(year) {
  file_name <- paste0("./latest_missing_data_fixed/raw-data/weather-data/tmeans_raw_data/weighted_area_raster_fips_tmean_daily_",
                      year,
                      ".rds")
  tmean <- readRDS(file_name) %>%
    filter(month == "06" | month == "07" | month =="08") %>%
    #                 month == "12" | month == "01" | month == "02" ~ "winter")) %>%
    group_by(fips) %>%
    summarise(tmean = mean(tmean)) %>%
    mutate(year = year) %>%
    select(fips, year, everything())
  
  # add year as a new row to the variable
  return(tmean)
}

ppt_means <- function(year) {
  file_name <- paste0("./latest_missing_data_fixed/raw-data/weather-data/ppt_raw_data/weighted_area_raster_fips_ppt_daily_",
                      year,
                      ".rds")
  ppt_mean <- readRDS(file_name) %>%
    filter(month == "06" | month == "07" | month =="08") %>%
    # month == "12" | month == "01" | month == "02" ~ "winter")) %>%
    group_by(fips) %>%
    summarise(ppt_mean = mean(ppt)) %>%
    mutate(year = year) %>%
    select(fips, year, everything())
  
  # add year as a new row to the variable
  return(ppt_mean)
}

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

# read raw data
seasonal_tmeans <- rbind(tmeans(1995), tmeans(1996), tmeans(1997), tmeans(1998),
                         tmeans(1999), tmeans(2000), tmeans(2001), tmeans(2002),
                         tmeans(2003), tmeans(2004), tmeans(2005), tmeans(2006),
                         tmeans(2007), tmeans(2008), tmeans(2009), tmeans(2010),
                         tmeans(2011), tmeans(2012), tmeans(2013), tmeans(2014),
                         tmeans(2015), tmeans(2016), tmeans(2017), tmeans(2018),
                         tmeans(2019), tmeans(2020))

seasonal_ppt_means <- rbind(ppt_means(1995), ppt_means(1996), ppt_means(1997), ppt_means(1998),
                            ppt_means(1999), ppt_means(2000), ppt_means(2001), ppt_means(2002),
                            ppt_means(2003), ppt_means(2004), ppt_means(2005), ppt_means(2006),
                            ppt_means(2007), ppt_means(2008), ppt_means(2009), ppt_means(2010),
                            ppt_means(2011), ppt_means(2012), ppt_means(2013), ppt_means(2014),
                            ppt_means(2015), ppt_means(2016), ppt_means(2017), ppt_means(2018),
                            ppt_means(2019), ppt_means(2020))

seasonal_means <- full_join(seasonal_tmeans, seasonal_ppt_means, by = c("fips" = "fips", 
                                                                        "year" = "year")) %>%
  fix_county_codes() %>%
  group_by(fips, year) %>%
  mutate(tmean = mean(tmean),
         ppt_mean = mean(ppt_mean)) %>%
  unique()

write.csv(seasonal_means, file = "./latest_missing_data_fixed/raw-data/weather-data/seasonal_means.csv")

missing_entries <- seasonal_means %>% filter_all(any_vars(is.na(.)))
View(missing_entries)
# no missing data