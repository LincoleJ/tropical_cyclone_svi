# plot tropical cyclone exposure from 1995 to 2018 (24 year course)
# as a binary variable for yearly exposure. 
# number of years with tropical cyclone exposure by county for 1995-2018

rm(list=ls())
library(readr)
library(tidyvserse)
library(dplyr)
library(ggplot2)
library(sf)
library(tigris)

# ensure fips code consistency over time
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

states_sf = states(cb = TRUE) %>%
  filter(GEOID != "02" & GEOID != "72" & GEOID != "69" & 
           GEOID != "78" & GEOID != "66" & GEOID != "60" &
           GEOID != "15")

counties_sf = counties(cb = TRUE) %>%
  filter(GEOID != "02" & GEOID != "72" & GEOID != "69" & 
           GEOID != "78" & GEOID != "66" & GEOID != "60" &
           GEOID != "15") %>%
  rename(fips = GEOID) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  summarise(geometry = st_union(geometry))

hurr_dat = read_csv("./latest_missing_data_fixed/raw-data/tc-data/hurr_dat.csv")[, -1] %>%
  fix_county_codes() %>%
  exclude_non_contiguous_states()
tc_exposure_sum = hurr_dat %>% mutate(total_exposure = rowSums(.[, 12:25])) %>%
  select(fips, total_exposure)
tc_sum_sf = merge(counties_sf, tc_exposure_sum, by.x = "fips", by.y = "fips")
tc_sum_sf = st_transform(tc_sum_sf, crs = 4326)
sf_counties = tc_sum_sf %>% select(-total_exposure)
tc_sum_sf$total_exposure = as.numeric(tc_sum_sf$total_exposure)
map_tc = ggplot() +
#  geom_sf(data = sf_counties, fill = "white", color = "gray99") +
  geom_sf(data = tc_sum_sf, aes(fill = total_exposure)) +
  scale_fill_gradient(name = "total number of exposure years",
                      low = "lightyellow",
                      high = "darkblue",
                      breaks = c(0, 1, 3, 5, 7, 9)) + 
  # add state borderline
  geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.35) +
  coord_sf(datum = NA) +
  theme_minimal() + 
  theme(legend.position = "bottom") +
  labs(title = "Tropical Cyclone Exposure Years") +
  theme(plot.title = element_text(size = 24, hjust = 0.5))

ggsave("./latest_missing_data_fixed/figures_tables/tc_figures/total_exposure.jpg",
       map_tc)


# ccombine the tc plot with SV
