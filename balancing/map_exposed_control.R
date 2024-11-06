rm(list=ls())
library("sf")
library("tidyverse")
library("dplyr")
library("tigris")
options(tigris_class = "sf")
library("ggplot2")
library("scales")

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

# states_sf <- states(cb = TRUE) %>%
#   filter(GEOID != "02" & GEOID != "72" & GEOID != "69" & 
#            GEOID != "78" & GEOID != "66" & GEOID != "60" &
#            GEOID != "15")
counties_sf = counties(cb = TRUE) %>%
  filter(GEOID != "02" & GEOID != "72" & GEOID != "69" & 
           GEOID != "78" & GEOID != "66" & GEOID != "60" &
           GEOID != "15") %>%
  rename(fips = GEOID) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  summarise(geometry = st_union(geometry))


for (year in 2005:2018) {
  df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds"))
  weights = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
                               pattern = paste(year),
                               full.names = TRUE))
  df_weight = df %>% drop_na() %>%
    filter(fips != 35039,
           fips != 42033)
  df_weight$weight = df_weight$exposed *weights$weights.1 + 
    (1-df_weight$exposed)*weights$weights.0
  df_weight = df_weight %>% mutate(weight = case_when(is.nan(weight) ~ 1,
                                                      TRUE ~ weight)) # fix infinite weights to be 1
  df_weight2 = merge(counties_sf, df_weight, by.x = "fips", by.y = "fips")
  df_weight2 = st_transform(df_weight2, crs = 4326)
  df_weight2 = df_weight2 %>% 
    mutate(control = case_when(weight != 1 & weight != 0 ~ weight))
  sf_counties = df_weight2 %>% select(fips, geometry)
  exposed_counties = df_weight2 %>% filter(weight == 1)
  exposed_counties$exposed = " "
  control_counties = df_weight2 %>% filter(!is.na(control)) 
  #high_wt_ctrls = control_counties %>% filter(control >= 0.1 * max(control_counties$control))
  sc_wt_cff = max(control_counties$control) * 0.1
  map_p = ggplot() +
    geom_sf(data = sf_counties, fill = "white", color = "black") +
    geom_sf(data = exposed_counties, aes(fill = exposed)) +
    scale_fill_manual(name = "Exposed", values = c(" " = "red")) +
    ggnewscale::new_scale_fill() +
    geom_sf(data = control_counties, aes(fill = control)) +
    scale_fill_gradientn(name = "Control",
                         colors = c("aliceblue", "deepskyblue", "darkblue"),
                         values = scales::rescale(c(min(control_counties$control), 
                                                    sc_wt_cff,
                                                    max(control_counties$control)))) + 
    coord_sf(datum = NA)  +
    theme_minimal() + 
    labs(title = year) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 20, face = "bold", hjust=0.5),
          legend.box = "vertical",
          legend.box.just = "left",
          legend.spacing.y = unit(0.05, "cm")) 
  ggsave(file.path("./latest_missing_data_fixed/balancing/map_figures", paste0("map_region", year, ".jpeg")), 
         map_p)
}

# test code 
year = 2015
df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds"))
weights = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
                             pattern = paste(year),
                             full.names = TRUE))
df_weight = df %>% drop_na() %>%
  filter(fips != 35039,
         fips != 42033)
df_weight$weight = df_weight$exposed *weights$weights.1 + 
  (1-df_weight$exposed)*weights$weights.0
df_weight = df_weight %>% mutate(weight = case_when(is.nan(weight) ~ 1,
                                                    TRUE ~ weight)) # fix infinite weights to be 1
df_weight2 = merge(counties_sf, df_weight, by.x = "fips", by.y = "fips")
df_weight2 = st_transform(df_weight2, crs = 4326)
df_weight2 = df_weight2 %>% 
  mutate(control = case_when(weight != 1 & weight != 0 ~ weight))
sf_counties = df_weight2 %>% select(fips, geometry)
exposed_counties = df_weight2 %>% filter(weight == 1)
exposed_counties$exposed = " "
control_counties = df_weight2 %>% filter(!is.na(control)) 
sc_wt_cff = max(control_counties$control) * 0.1
map_p = ggplot() +
  geom_sf(data = sf_counties, fill = "white", color = "black") +
  geom_sf(data = exposed_counties, aes(fill = exposed)) +
  scale_fill_manual(name = "Exposed", values = c(" " = "red")) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = control_counties, aes(fill = control)) +
  scale_fill_gradientn(name = "Control",
                       colors = c("aliceblue", "deepskyblue", "darkblue"),
                       values = scales::rescale(c(min(control_counties$control), 
                                                  sc_wt_cff,
                                                  max(control_counties$control)))) + 
  coord_sf(datum = NA)  +
  theme_minimal() + 
  labs(title = year) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 24, hjust=0.5),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(0.05, "cm")) 


#### test code for show exposed regions
year = 2012
df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds"))
weights = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
                             pattern = paste(year),
                             full.names = TRUE))
df_weight = df %>% drop_na() %>%
  filter(fips != 35039,
         fips != 42033)
df_weight$weight = df_weight$exposed *weights$weights.1 + 
  (1-df_weight$exposed)*weights$weights.0
df_weight = df_weight %>% mutate(weight = case_when(is.nan(weight) ~ 1,
                                                    TRUE ~ weight)) # fix infinite weights to be 1
df_weight2 = merge(counties_sf, df_weight, by.x = "fips", by.y = "fips")
df_weight2 = st_transform(df_weight2, crs = 4326)
df_weight2 = df_weight2 %>% 
  mutate(control = case_when(weight != 1 & weight != 0 ~ weight))
sf_counties = df_weight2 %>% select(fips, geometry)
exposed_counties = df_weight2 %>% filter(weight == 1)
control_counties = df_weight2 %>% filter(!is.na(control)) 
#high_wt_ctrls = control_counties %>% filter(control >= 0.1 * max(control_counties$control))
sc_wt_cff = max(control_counties$control) * 0.1
map_p = ggplot() +
  geom_sf(data = sf_counties, fill = "white", color = "black") +
  geom_sf(data = exposed_counties, fill = "red") +
  geom_sf(data = control_counties, aes(fill = control)) +
  scale_fill_gradientn(name = "Control",
                       colors = c("aliceblue", "deepskyblue", "darkblue"),
                       values = scales::rescale(c(min(control_counties$control), 
                                                  sc_wt_cff,
                                                  max(control_counties$control)))) + 
  coord_sf(datum = NA)  +
  theme_minimal() + 
  labs(title = year) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5)) 

##############
###### check distribution of weights
for (year in 2005:2018) {
  df = readRDS(paste0("./latest_missing_data_fixed/processed-data/df_", year, ".rds"))
  weights = readRDS(list.files(file.path("./latest_missing_data_fixed/balancing/weights-data-upd"),
                               pattern = paste(year),
                               full.names = TRUE))
  df_weight = df %>% drop_na() %>%
    filter(fips != 35039,
           fips != 42033)
  df_weight$weight = df_weight$exposed *weights$weights.1 +
    (1-df_weight$exposed)*weights$weights.0
  df_weight = df_weight %>% mutate(weight = case_when(is.nan(weight) ~ 1,
                                                      TRUE ~ weight)) %>% 
    mutate(control = case_when(weight != 1 & weight != 0 ~ weight))
  control_counties = df_weight %>% 
    filter(!is.na(control)) %>%
    mutate(normalized_weight = weight / max(weight))
  png(paste0("./latest_missing_data_fixed/balancing/sensitivity/dist_weights_ctrls/", year, "_control_weights.png"))
  hist(control_counties$normalized_weight)
  dev.off()
}
## with exception of 2009, for other exposure years the majority of weights are
## within 10% of the maximum weight
## as such, to represent the maps 