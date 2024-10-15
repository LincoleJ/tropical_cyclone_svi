# plot svi for contiguous U.S. in 2020, then for 2018, 2016, 2014, 2010 (figure 1)
# we plot the SVI as a categorical variable, definition consistent with
# https://www.atsdr.cdc.gov/placeandhealth/svi/faq_svi.html
## "we classify data using quartiles 
## (0 to 0.2500, 0.2501 to 0.5000, 0.5001 to 0.7500, 0.7501 to 1.0) 
## and indicate that the classification goes from least vulnerable to most vulnerable"
rm(list=ls())
print(Sys.time())
library("sf")
library("tidyverse")
library("dplyr")
library("tigris")
library("readr")
options(tigris_class = "sf")
library("ggplot2")

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

counties_sf = counties(cb = TRUE) %>%
  filter(GEOID != "02" & GEOID != "72" & GEOID != "69" & 
           GEOID != "78" & GEOID != "66" & GEOID != "60" &
           GEOID != "15") %>%
  rename(fips = GEOID) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  summarise(geometry = st_union(geometry))

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

svi_dat = read_csv("./latest_missing_data_fixed/raw-data/svi-data/svi_dat.csv")[, -1] %>%
  fix_county_codes() %>%
  exclude_non_contiguous_states()
  
outcome_years = c(2010, 2014, 2016, 2018, 2020)
for (i in 1:length(outcome_years)) {
  outcome_year = outcome_years[i]
  svi = svi_dat[, c("fips", outcome_year)]
  colnames(svi) = c("fips", "svi")
  
  svi_sf = merge(counties_sf, svi, by.x = "fips", by.y = "fips")
  svi_sf = st_transform(svi_sf, crs = 4326)
  sf_counties = svi_sf %>% select(-svi)
  svi_sf = svi_sf %>% mutate(vul_level = case_when(svi > 0.75 ~ "High",
                                                   svi > 0.5 & svi <= 0.75 ~ "Medium-High",
                                                   svi > 0.25 & svi <= 0.5 ~ "Low-Medium",
                                                   svi <= 0.25 ~ "Low"),
                             vul_level = factor(vul_level, 
                                                levels = c("Low", "Low-Medium", "Medium-High", "High"))) 
  svi_map = ggplot() +
  #  geom_sf(data = sf_counties, fill = "white", color = "black") +
    geom_sf(data = svi_sf, aes(fill = vul_level)) +
    scale_fill_manual(values = c("Low" = "lightyellow",
                                 "Low-Medium" = "palegreen",
                                 "Medium-High" = "turquoise3",
                                 "High" = "steelblue"),
                      labels = c("Low", "Low-Medium", "Medium-High", "High"),
                      name = "Level of Social Vulnerability") +
    geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.3) +
    coord_sf(datum = NA) +
    theme_minimal() +
    labs(title = paste0("Social Vulnerability Index in ", outcome_year)) +
    theme(plot.title = element_text(size = 24, hjust = 0.5)) + 
    theme(legend.position = "bottom") 
  ggsave(paste0("./latest_missing_data_fixed/figures_tables/svi", outcome_year, "_svi.jpg"))
}

# test code for 2020
outcome_year = 2020
svi = svi_dat[, c("fips", outcome_year)]
colnames(svi) = c("fips", "svi")

svi_sf = merge(counties_sf, svi, by.x = "fips", by.y = "fips")
svi_sf = st_transform(svi_sf, crs = 4326)
sf_counties = svi_sf %>% select(-svi)
svi_sf = svi_sf %>% mutate(vul_level = case_when(svi > 0.75 ~ "High",
                                                 svi > 0.5 & svi <= 0.75 ~ "Medium-High",
                                                 svi > 0.25 & svi <= 0.5 ~ "Low-Medium",
                                                 svi <= 0.25 ~ "Low"),
                           vul_level = factor(vul_level, 
                                              levels = c("Low", "Low-Medium", "Medium-High", "High"))) 
svi_map = ggplot() +
  #  geom_sf(data = sf_counties, fill = "white", color = "black") +
  geom_sf(data = svi_sf, aes(fill = vul_level)) +
  scale_fill_manual(values = c("Low" = "lightyellow",
                               "Low-Medium" = "palegreen",
                               "Medium-High" = "turquoise3",
                               "High" = "steelblue"),
                    labels = c("Low", "Low-Medium", "Medium-High", "High"),
                    name = "Level of Social Vulnerability") +
 # geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.3) +
  coord_sf(datum = NA) +
  theme_minimal() +
#  labs(title = paste0("Social Vulnerability Index in ", outcome_year)) +
  theme(plot.title = element_text(size = 24, hjust = 0.5)) + 
  theme(legend.position = "bottom") 
ggsave(paste0("./latest_missing_data_fixed/figures_tables/svi_figures/", outcome_year, "_svi.jpg"),
       svi_map)


# combine SVI plot with TC plot
library(cowplot)
plot_grid(map_tc, svi_map, align = "v", ncol = 1, labels = c("A", "B"))
