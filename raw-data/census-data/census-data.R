rm(list=ls())
library(purrr)
library(readr)
library(tidyverse)
library(dplyr)

############################----helper functions----############################
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

#####################-------covariates annual infer-------#####################
###
### per capita income
# fix county codes to standardize BEA CAINC1 for per capita income
# source: https://apps.bea.gov/itable/?ReqID=70&step=1&_gl=1*3ckuj0*_ga*MTE5ODgzNTEwNC4xNzE3Njk1MjQ0*_ga_J4698JNNFT*MTcxNzY5NTI0My4xLjEuMTcxNzY5NTY2Mi40OC4wLjA.#eyJhcHBpZCI6NzAsInN0ZXBzIjpbMSwyOSwyNSwzMSwyNiwyNywzMF0sImRhdGEiOltbIlRhYmxlSWQiLCIyMCJdLFsiTWFqb3JfQXJlYSIsIjQiXSxbIlN0YXRlIixbIjUxMDAwIl1dLFsiQXJlYSIsWyJYWCJdXSxbIlN0YXRpc3RpYyIsWyItMSJdXSxbIlVuaXRfb2ZfbWVhc3VyZSIsIkxldmVscyJdLFsiWWVhciIsWyIyMDIyIl1dLFsiWWVhckJlZ2luIiwiLTEiXSxbIlllYXJfRW5kIiwiLTEiXV19
va_fips_mapping = list(
  "51901" = c("51003", "51540"),
  "51903" = c("51005", "51580"),
  "51907" = c("51015", "51790", "51820"),
  "51911" = c("51031", "51680"),
  "51913" = c("51035", "51640"),
  "51918" = c("51053", "51570", "51730"),
  "51919" = c("51059", "51600", "51610"),
  "51921" = c("51069", "51840"),
  "51923" = c("51081", "51595"),
  "51929" = c("51089", "51690"),
  "51931" = c("51095", "51830"),
  "51933" = c("51121", "51750"),
  "51939" = c("51143", "51590"),
  "51941" = c("51149", "51670"),
  "51942" = c("51153", "51683", "51685"),
  "51944" = c("51161", "51775"),
  "51945" = c("51163", "51530", "51678"),
  "51947" = c("51165", "51660"),
  "51949" = c("51175", "51620"),
  "51951" = c("51177", "51630"),
  "51953" = c("51191", "51520"), 
  "51955" = c("51195", "51720"),
  "51958" = c("51199", "51735")
)

per_capita_income = list.files(path = "./latest_missing_data_fixed/raw-data/census-data/covariates_annual_infer/per_capita_income",
                               pattern = "*.csv", full.name = TRUE) %>% 
  read_csv() %>%
  rename(pci = income_adjusted_to_2000) %>%
  drop_na() 
pci_unchanged = per_capita_income %>% filter(!fips %in% names(va_fips_mapping))
pci_to_change = per_capita_income %>% 
  filter(fips %in% names(va_fips_mapping))  # for some reason 51770 is available in original processed dataset
mapping_df <- data.frame(
  fips = rep(names(va_fips_mapping), sapply(va_fips_mapping, length)),
  new_fips = unlist(va_fips_mapping)
) 
pci_changed = pci_to_change %>%
  left_join(mapping_df, by = "fips") %>%
  mutate(fips = new_fips) %>%
  select(-new_fips)
pci = rbind(pci_changed, pci_unchanged) %>%
  pivot_wider(names_from = year,
              values_from = pci,
              names_prefix = "pci_")

# poverty
pct_pov = list.files(path = "./latest_missing_data_fixed/raw-data/census-data/covariates_annual_infer/poverty",
                     pattern = "*.csv", full.name = TRUE) %>% read_csv() %>%
  drop_na() %>%
  rename(pct_pov = poverty_rate) %>%
  pivot_wider(names_from = year,
              values_from = pct_pov,
              names_prefix = "pct_pov_") %>%
  mutate(pct_pov_1996 = (pct_pov_1995 + pct_pov_1997) / 2) 

# percent high school graduate 
pct_hsg = list.files(path = "./latest_missing_data_fixed/raw-data/census-data/covariates_annual_infer/education",
                     pattern = "*.csv", full.names = TRUE) %>% read_csv() %>%
  drop_na() %>%
  rename(pct_hsg = percent_high_school_or_above) %>%
  pivot_wider(names_from = year,
              values_from = pct_hsg,
              names_prefix = "pct_hsg_")

#########################-------cdc monthly infer-------#########################

# pct_black
pct_black = list.files(path = "./raw-data/census-data/cdc_population_monthly_infer/pct_black",
                       pattern = "*.csv", full.name = TRUE) %>% 
  read_csv() %>% filter(month == 6) %>% select(-month) %>% 
  group_by(year, fips, race) %>% 
  summarise(pop = sum(pop)) %>% 
  group_by(year, fips) %>% 
  mutate(total_pop = sum(pop)) %>% ungroup() %>% 
  filter(race == "Black or African American") %>%
  mutate(pct_black = pop / total_pop) %>%
  select(year, fips, pct_black) %>% drop_na() %>%
  pivot_wider(names_from = year,
              values_from = pct_black,
              names_prefix = "pct_black_")

# percent hispanic & male/female ratio
sex_ethnicity = list.files(path = "./raw-data/census-data/cdc_population_monthly_infer/pct_hispanic",
                           pattern = "*.csv", full.name = TRUE) %>% 
  read_csv() %>% 
  filter(month == 6) %>% select(-month) 

# pct_hispanic
pct_hispanic = sex_ethnicity %>% 
  group_by(year, fips, ethnicity) %>% summarise(pop = sum(pop), .groups = "keep") %>%
  group_by(fips, year) %>%
  mutate(total_pop = sum(pop)) %>% ungroup() %>%
  filter(ethnicity == "Hispanic or Latino") %>%
  mutate(pct_hispanic = pop / total_pop) %>% 
  drop_na() %>%
  select(year, fips, pct_hispanic) %>%
  pivot_wider(names_from = year, 
              values_from = pct_hispanic,
              names_prefix = "pct_hispanic_")

# male-female ratio
mf_ratio = sex_ethnicity %>% group_by(year, fips, sex) %>%
  summarise(pop = sum(pop), .groups = "keep") %>%
  pivot_wider(names_from = sex,
              values_from = pop,
              names_prefix = "pop_") %>%
  mutate(mf_ratio = pop_1 / pop_2) %>%
  select(-pop_1, -pop_2) %>%
  pivot_wider(names_from = year,
              values_from = mf_ratio,
              names_prefix = "mf_ratio_")

# pct_blw_15 & pct_abv_65
age = list.files(path = "./raw-data/census-data/cdc_population_monthly_infer/age", 
                 pattern = "*.csv", full.name = TRUE) %>% read_csv() %>%
  filter(month == 6) %>% select(-month) %>%
  group_by(year, fips, age) %>%
  summarise(pop = sum(pop), .groups = "keep") %>%
  group_by(year, fips) %>%
  mutate(total_pop = sum(pop)) %>% ungroup()

pct_blw_15 = age %>% filter(age == 0 | age == 5) %>%
  group_by(year, fips) %>% mutate(blw_15 = sum(pop)) %>% ungroup() %>%
  select(-age, -pop) %>% distinct() %>%
  mutate(pct_blw_15 = blw_15 / total_pop) %>%
  select(fips, year, pct_blw_15) %>%
  pivot_wider(names_from = year,
              values_from = pct_blw_15,
              names_prefix = "pct_youth_")

pct_abv_65 = age %>% filter(age == 65 | age == 75 | age == 85) %>%
  group_by(year, fips) %>% mutate(abv_65 = sum(pop)) %>% ungroup() %>%
  select(-age, -pop) %>% distinct() %>%
  mutate(pct_abv_65 = abv_65 / total_pop) %>%
  select(fips, year, pct_abv_65) %>%
  pivot_wider(names_from = year,
              values_from = pct_abv_65,
              names_prefix = "pct_elderly_")

# population density
pop = sex_ethnicity %>% group_by(fips, year) %>%
  summarise(pop = sum(pop), .groups = "keep") 
cty_area <- counties(cb = TRUE) %>%
  # ALAND - land area in square meters
  # area - land area in square miles
  mutate(area = ALAND / 2589988,
         FIPS = GEOID) %>%
  select(FIPS, area) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(area = sum(area)) %>%
  data.frame() %>%
  select(-geometry) %>%
  unique()

pop_density = merge(pop, cty_area, by = "fips") %>% 
  mutate(pop_density = pop / area) %>%
  select(fips, pop_density, year) %>%
  pivot_wider(names_from = year,
              values_from = pop_density,
              names_prefix = "pop_density_")

# keep county-level population for linear model in outcome analysis
pop_size = pop %>% 
  pivot_wider(names_from = year,
              values_from = pop,
              names_prefix = "pop_size_")

###############------------Conglomerating Data---------------###################
list_df = list(pct_hsg, pci, pct_pov, pct_black, pct_hispanic,
               mf_ratio, pct_blw_15, pct_abv_65, pop_density, pop_size)
census_df = reduce(list_df, full_join, by = "fips")

missing_id = apply(census_df, 1, function(x) any(is.na(x)))
missing_entries = census_df[missing_id, ]

census_df = exclude_non_contiguous_states(census_df) # %>% drop_na()

# check missing data
missing_id = apply(census_df, 1, function(x) any(is.na(x)))
missing_entries = census_df[missing_id, ]
# only missing entry: 08014, per capita income, from 1995 to 2001
# checked cleaned data on annual_covariates_infer, available for all counties 
# before applying fix_county_codes
# recalculate the values for 08014 from unadjusted datasets for 1995-2001

dat = data.frame()
for(year_current in 1995:2001) {
  
  print(year_current)
  
  filename_in = paste0("./latest_missing_data_fixed/raw-data/census-data/covariates_annual_infer/pci_unfixed_08014/",
                       "per_capita_income_adjusted_to_2000_", year_current, ".csv")
  
  dat_year = read_csv(filename_in) %>% 
    filter(fips == "08001" | fips == "08013" | fips == "08059" | fips == "08123")
  
  dat = data.table::rbindlist(list(dat, dat_year))
}

dat = fix_county_codes(dat)

dat = dat %>% group_by(year, fips) %>%
  summarise(income_adjusted_to_2000 = mean(income_adjusted_to_2000))

fixed_county_values_08014 = dat$income_adjusted_to_2000
ridx = which(census_df$fips == "08014")
census_df[ridx, "pci_1995"] = fixed_county_values_08014[1]
census_df[ridx, "pci_1996"] = fixed_county_values_08014[2]
census_df[ridx, "pci_1997"] = fixed_county_values_08014[3]
census_df[ridx, "pci_1998"] = fixed_county_values_08014[4]
census_df[ridx, "pci_1999"] = fixed_county_values_08014[5]
census_df[ridx, "pci_2000"] = fixed_county_values_08014[6]
census_df[ridx, "pci_2001"] = fixed_county_values_08014[7]

saveRDS(census_df, paste("./latest_missing_data_fixed/raw-data/census-data/census_df.rds"))
