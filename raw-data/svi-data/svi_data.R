print(Sys.time())
rm(list=ls())
library(tidyverse)
library(dplyr)

# processes available CDC/ATSDR SVI data from 2000, 2010; biennially 2014-2020
# the final `svi_year` data-sets contain county-level percentile rankings 
# (averaged among census-tracts) for the *total* 
# percentile ranking across the four domains.

#######################---------helper functions---------#######################
# read_svi returns the raw data of the given svi year.
read_svi <- function(year) {
  file_name <- paste0("./latest_missing_data_fixed/raw-data/svi-data/raw_data/SVI_",
                      year,
                      "_US_county.csv")
  dat <- read_csv(file_name)
  return(dat)
}

# fix counties to be consistent across time, 
# mutually exclusive, and collectively exhaustive
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

#######################----------read raw data---------########################
# for all years of SVI, for duplicate counties due to merger 
# (i.e., 30113 and 08014, the mean SVI was substituted)
svi_2022 <- read_svi(2022) %>%
  select(FIPS, RPL_THEMES) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEMES = mean(RPL_THEMES)) %>%
  rename("2022" = RPL_THEMES) %>%
  unique()

svi_2000 <- read_svi(2000) %>%
  select(STCOFIPS, USTP) %>%
  rename(fips = STCOFIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(USTP = mean(USTP)) %>%
  rename("2000" = USTP) %>%
  unique()

svi_2010 <- read_svi(2010) %>%
  select(FIPS, R_PL_THEMES) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(R_PL_THEMES = mean(R_PL_THEMES)) %>%
  rename("2010" = R_PL_THEMES) %>%
  unique()

svi_2014 <- read_svi(2014) %>%
  select(FIPS, RPL_THEMES) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEMES = mean(RPL_THEMES)) %>%
  rename("2014" = RPL_THEMES) %>%
  unique()

svi_2016 <- read_svi(2016) %>%
  select(FIPS, RPL_THEMES) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEMES = mean(RPL_THEMES)) %>%
  rename("2016" = RPL_THEMES) %>%
  unique()

svi_2018 <- read_svi(2018) %>%
  select(FIPS, RPL_THEMES) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEMES = mean(RPL_THEMES)) %>%
  rename("2018" = RPL_THEMES) %>%
  unique()

svi_2020 <- read_svi(2020) %>%
  select(FIPS, RPL_THEMES) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEMES = mean(RPL_THEMES)) %>%
  rename("2020" = RPL_THEMES) %>%
  unique()

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

list_svi <- list(svi_2000, svi_2010, svi_2014, svi_2016, svi_2018, svi_2020, svi_2022)
svi_dat <- list_svi %>% reduce(full_join, by = "fips") %>% exclude_non_contiguous_states()
missing_rows <- svi_dat[rowSums(is.na(svi_dat)) > 0, ]

# fips code changes, effective 8.13.2024
# https://developer.ap.org/ap-elections-api/docs/CT_FIPS_Codes_forPlanningRegions.htm
# the change in SVI is almost impossible to map out due to inter-changes of planning
# regions within each FIPS code
fix_county_codes_ct = function(dat){
  dat = dat %>%
    mutate(fips = case_when(
      # Alaska FIPS codes to do here when have time
      fips== 'XXXXX' ~ 'XXXXX',
      fips == "09013" | fips == "09003" ~ "09110", 
      fips == "09001" ~ "09120",
      fips == "09007" ~ "09130",
      TRUE ~ fips
    ))
}

write.csv(svi_dat, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_dat_full.csv")

#####################---------SVI Subcategories---------########################
# 2000
# XXG1TP - Socioeconomic Domain Total Percentile Ranking
# XXG2TP - Household Composition & Disability Total Percentile Ranking
# XXG3TP - Minority Status/Language Domain Total Percentile Ranking
# XXG4TP - Housing/Transportation Domain Total Percentile Ranking

# 2010
# R_PL_THEME1 - Percentile ranking for Socioeconomic theme
# R_PL_THEME2 - Percentile ranking for Household Composition theme
# R_PL_THEME3 - Percentile ranking for Minority Status/Language theme
# R_PL_THEME4	- Percentile ranking for Housing/Transportation theme

# 2014
# Socioeconomic theme – RPL_THEME1
# Household Composition and Disability – RPL_THEME2
# Minority Status & Language – RPL_THEME3
# Housing & Transportation – RPL_THEME4

# 2016
# Socioeconomic – RPL_THEME1
# Household Composition & Disability – RPL_THEME2
# Minority Status & Language – RPL_THEME3
# Housing Type & Transportation – RPL_THEME4

# 2018
# Socioeconomic – RPL_THEME1
# Household Composition & Disability – RPL_THEME2
# Minority Status & Language – RPL_THEME3
# Housing Type & Transportation – RPL_THEME4

# 2020
# Socioeconomic Status – RPL_THEME1
# Household Characteristics – RPL_THEME2
# Racial & Ethnic Minority Status – RPL_THEME3
# Housing Type & Transportation – RPL_THEME4

# housing domain
svi_2000_housing <- read_svi(2000) %>%
  select(STCOFIPS, USG4TP) %>%
  rename(fips = STCOFIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(USG4TP = mean(USG4TP)) %>%
  rename("2000" = USG4TP) %>%
  unique()

svi_2010_housing <- read_svi(2010) %>%
  select(FIPS, R_PL_THEME4) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(R_PL_THEME4 = mean(R_PL_THEME4)) %>%
  rename("2010" = R_PL_THEME4) %>%
  unique()

svi_2014_housing <- read_svi(2014) %>%
  select(FIPS, RPL_THEME4) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME4 = mean(RPL_THEME4)) %>%
  rename("2014" = RPL_THEME4) %>%
  unique()

svi_2016_housing <- read_svi(2016) %>%
  select(FIPS, RPL_THEME4) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME4 = mean(RPL_THEME4)) %>%
  rename("2016" = RPL_THEME4) %>%
  unique()

svi_2018_housing <- read_svi(2018) %>%
  select(FIPS, RPL_THEME4) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME4 = mean(RPL_THEME4)) %>%
  rename("2018" = RPL_THEME4) %>%
  unique()

svi_2020_housing <- read_svi(2020) %>%
  select(FIPS, RPL_THEME4) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME4 = mean(RPL_THEME4)) %>%
  rename("2020" = RPL_THEME4) %>%
  unique()

list_svi_housing <- list(svi_2000_housing, svi_2010_housing, svi_2014_housing, 
                         svi_2016_housing, svi_2018_housing, svi_2020_housing)
svi_dat_housing <- list_svi_housing %>% reduce(full_join, by = "fips")
write.csv(svi_dat_housing, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_dat_housing.csv")

# minority domain
svi_2000_minority <- read_svi(2000) %>%
  select(STCOFIPS, USG3TP) %>%
  rename(fips = STCOFIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(USG3TP = mean(USG3TP)) %>%
  rename("2000" = USG3TP) %>%
  unique()

svi_2010_minority <- read_svi(2010) %>%
  select(FIPS, R_PL_THEME3) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(R_PL_THEME3 = mean(R_PL_THEME3)) %>%
  rename("2010" = R_PL_THEME3) %>%
  unique()

svi_2014_minority <- read_svi(2014) %>%
  select(FIPS, RPL_THEME3) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME3 = mean(RPL_THEME3)) %>%
  rename("2014" = RPL_THEME3) %>%
  unique()

svi_2016_minority <- read_svi(2016) %>%
  select(FIPS, RPL_THEME3) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME3 = mean(RPL_THEME3)) %>%
  rename("2016" = RPL_THEME3) %>%
  unique()

svi_2018_minority <- read_svi(2018) %>%
  select(FIPS, RPL_THEME3) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME3 = mean(RPL_THEME3)) %>%
  rename("2018" = RPL_THEME3) %>%
  unique()

svi_2020_minority <- read_svi(2020) %>%
  select(FIPS, RPL_THEME3) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME3 = mean(RPL_THEME3)) %>%
  rename("2020" = RPL_THEME3) %>%
  unique()

list_svi_minority <- list(svi_2000_minority, svi_2010_minority, svi_2014_minority,
                          svi_2016_minority, svi_2018_minority, svi_2020_minority)
svi_dat_minority <- list_svi_minority %>% reduce(full_join, by = "fips")
write.csv(svi_dat_minority, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_dat_minority.csv")

# socioeconomics domain
svi_2000_ses <- read_svi(2000) %>%
  select(STCOFIPS, USG1TP) %>%
  rename(fips = STCOFIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(USG1TP = mean(USG1TP)) %>%
  rename("2000" = USG1TP) %>%
  unique()

svi_2010_ses <- read_svi(2010) %>%
  select(FIPS, R_PL_THEME1) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(R_PL_THEME1 = mean(R_PL_THEME1)) %>%
  rename("2010" = R_PL_THEME1) %>%
  unique()

svi_2014_ses <- read_svi(2014) %>%
  select(FIPS, RPL_THEME1) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME1 = mean(RPL_THEME1)) %>%
  rename("2014" = RPL_THEME1) %>%
  unique()

svi_2016_ses <- read_svi(2016) %>%
  select(FIPS, RPL_THEME1) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME1 = mean(RPL_THEME1)) %>%
  rename("2016" = RPL_THEME1) %>%
  unique()

svi_2018_ses <- read_svi(2018) %>%
  select(FIPS, RPL_THEME1) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME1 = mean(RPL_THEME1)) %>%
  rename("2018" = RPL_THEME1) %>%
  unique()

svi_2020_ses <- read_svi(2020) %>%
  select(FIPS, RPL_THEME1) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME1 = mean(RPL_THEME1)) %>%
  rename("2020" = RPL_THEME1) %>%
  unique()

list_svi_ses <- list(svi_2000_ses, svi_2010_ses, svi_2014_ses, 
                     svi_2016_ses, svi_2018_ses, svi_2020_ses)
svi_dat_ses <- list_svi_ses %>% reduce(full_join, by = "fips")
write.csv(svi_dat_ses, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_dat_ses.csv")

# household domain
svi_2000_hhd <- read_svi(2000) %>%
  select(STCOFIPS, USG2TP) %>%
  rename(fips = STCOFIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(USG2TP = mean(USG2TP)) %>%
  rename("2000" = USG2TP) %>%
  unique()

svi_2010_hhd <- read_svi(2010) %>%
  select(FIPS, R_PL_THEME2) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(R_PL_THEME2 = mean(R_PL_THEME2)) %>%
  rename("2010" = R_PL_THEME2) %>%
  unique()

svi_2014_hhd <- read_svi(2014) %>%
  select(FIPS, RPL_THEME2) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME2 = mean(RPL_THEME2)) %>%
  rename("2014" = RPL_THEME2) %>%
  unique()

svi_2016_hhd <- read_svi(2016) %>%
  select(FIPS, RPL_THEME2) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME2 = mean(RPL_THEME2)) %>%
  rename("2016" = RPL_THEME2) %>%
  unique()

svi_2018_hhd <- read_svi(2018) %>%
  select(FIPS, RPL_THEME2) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME2 = mean(RPL_THEME2)) %>%
  rename("2018" = RPL_THEME2) %>%
  unique()

svi_2020_hhd <- read_svi(2020) %>%
  select(FIPS, RPL_THEME2) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(RPL_THEME2 = mean(RPL_THEME2)) %>%
  rename("2020" = RPL_THEME2) %>%
  unique()

list_svi_hhd <- list(svi_2000_hhd, svi_2010_hhd, svi_2014_hhd, 
                     svi_2016_hhd, svi_2018_hhd, svi_2020_hhd)
svi_dat_hhd <- list_svi_hhd %>% reduce(full_join, by = "fips")
write.csv(svi_dat_hhd, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_dat_hhd.csv")

missing_entries <- svi_dat %>% filter_all(any_vars(is.na(.)))
View(missing_entries)
# missing data only in Alaska, which is irrelevant for our analysis

################################################################################
#####################-----------extract sub-hhd-----------#####################
################################################################################
# definition of variables for 2010 SVI Household domain
# PL_AGE65: percentile of proportion of persons aged 65 or older
# PL_AGE17: Percentile of the proportion of persons aged 17 and younger
# PL_SNGPRNT: Percentile of the proportion of single parent households with children under 18
svi_2010 = read_svi(2010)
svi_2010_sngprnt <- svi_2010 %>%
  select(FIPS, PL_SNGPRNT) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(PL_SNGPRNT = mean(PL_SNGPRNT)) %>%
  rename("2010" = PL_SNGPRNT) %>%
  unique()
svi_2010_age65 <- svi_2010 %>%
  select(FIPS, PL_AGE65) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(PL_AGE65 = mean(PL_AGE65)) %>%
  rename("2010" = PL_AGE65) %>%
  unique()
svi_2010_age17 <- svi_2010 %>%
  select(FIPS, PL_AGE17) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(PL_AGE17 = mean(PL_AGE17)) %>%
  rename("2010" = PL_AGE17) %>%
  unique()

# definition of variables for 2014 SVI Household domain
# EPL_AGE65: percentile peercentage of persons aged 65 and older estimate
# EPL_AGE17: percentile percentage of persons aged 17 and younger estimate
# EPL_DISABL: percentile percentage of civilian noninstitutionalized population with a disability estimate
# EPL_SNGPNT: percentile percentage of single parent households with children under 18 estimate
svi_2014 = read_svi(2014)
svi_2014_disabl <- svi_2014 %>%
  select(FIPS, EPL_DISABL) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_DISABL = mean(EPL_DISABL)) %>%
  rename("2014" = EPL_DISABL) %>%
  unique()
svi_2014_sngpnt <- svi_2014 %>%
  select(FIPS, EPL_SNGPNT) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_SNGPNT = mean(EPL_SNGPNT)) %>%
  rename("2014" = EPL_SNGPNT) %>%
  unique()
svi_2014_age65 <- svi_2014 %>%
  select(FIPS, EPL_AGE65) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE65 = mean(EPL_AGE65)) %>%
  rename("2014" = EPL_AGE65) %>%
  unique()
svi_2014_age17 <- svi_2014 %>%
  select(FIPS, EPL_AGE17) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE17 = mean(EPL_AGE17)) %>%
  rename("2014" = EPL_AGE17) %>%
  unique()

# definition of varaibles for 2016 SVI Household domain
# EPL_AGE65
# EPL_AGE17
# EPL_DISABL
# EPL_SNGPNT
svi_2016 = read_svi(2016)
svi_2016_disabl <- svi_2016 %>%
  select(FIPS, EPL_DISABL) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_DISABL = mean(EPL_DISABL)) %>%
  rename("2016" = EPL_DISABL) %>%
  unique()
svi_2016_sngpnt <- svi_2016 %>%
  select(FIPS, EPL_SNGPNT) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_SNGPNT = mean(EPL_SNGPNT)) %>%
  rename("2016" = EPL_SNGPNT) %>%
  unique()
svi_2016_age65 <- svi_2016 %>%
  select(FIPS, EPL_AGE65) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE65 = mean(EPL_AGE65)) %>%
  rename("2016" = EPL_AGE65) %>%
  unique()
svi_2016_age17 <- svi_2016 %>%
  select(FIPS, EPL_AGE17) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE17 = mean(EPL_AGE17)) %>%
  rename("2016" = EPL_AGE17) %>%
  unique()

# definition of variables for 2018 SVI Household domain
# EPL_AGE65: percentile percentage of persons aged 65 and older estimate
# EPL_AGE17: Percentile percentage of persons aged 17 and younger estimate
# EPL_DISABL: percentile percentage of civilian noninstiuttionalized population with a disability estimate
# EPL_SNGPNT: percentage percentage of single parent households with children under 18 estimate
svi_2018 = read_svi(2018)
svi_2018_disabl <- svi_2018 %>%
  select(FIPS, EPL_DISABL) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_DISABL = mean(EPL_DISABL)) %>%
  rename("2018" = EPL_DISABL) %>%
  unique()
svi_2018_sngpnt <- svi_2018 %>%
  select(FIPS, EPL_SNGPNT) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_SNGPNT = mean(EPL_SNGPNT)) %>%
  rename("2018" = EPL_SNGPNT) %>%
  unique()
svi_2018_age65 <- svi_2018 %>%
  select(FIPS, EPL_AGE65) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE65 = mean(EPL_AGE65)) %>%
  rename("2018" = EPL_AGE65) %>%
  unique()
svi_2018_age17 <- svi_2018 %>%
  select(FIPS, EPL_AGE17) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE17 = mean(EPL_AGE17)) %>%
  rename("2018" = EPL_AGE17) %>%
  unique()

# definition of variables for 2020 SVI Household domain
# EPL_AGE65
# EPL_AGE17
# EPL_DISABL
# EPL_SNGPNT
# EPL_LIMENG: percentile of percentage of persons (age 5+) who speak English "less than well"

svi_2020 = read_svi(2020)
svi_2020_disabl <- svi_2020 %>%
  select(FIPS, EPL_DISABL) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_DISABL = mean(EPL_DISABL)) %>%
  rename("2020" = EPL_DISABL) %>%
  unique()
svi_2020_sngpnt <- svi_2020 %>%
  select(FIPS, EPL_SNGPNT) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_SNGPNT = mean(EPL_SNGPNT)) %>%
  rename("2020" = EPL_SNGPNT) %>%
  unique()
svi_2020_age65 <- svi_2020 %>%
  select(FIPS, EPL_AGE65) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE65 = mean(EPL_AGE65)) %>%
  rename("2020" = EPL_AGE65) %>%
  unique()
svi_2020_age17 <- svi_2020 %>%
  select(FIPS, EPL_AGE17) %>%
  rename(fips = FIPS) %>%
  fix_county_codes() %>%
  group_by(fips) %>%
  mutate(EPL_AGE17 = mean(EPL_AGE17)) %>%
  rename("2020" = EPL_AGE17) %>%
  unique()

list_svi_age17 <- list(svi_2010_age17, svi_2014_age17, 
                     svi_2016_age17, svi_2018_age17, svi_2020_age17)
svi_dat_age17 <- list_svi_age17 %>% reduce(full_join, by = "fips")
write.csv(svi_dat_age17, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_age17.csv")

list_svi_age65 <- list(svi_2010_age65, svi_2014_age65, 
                       svi_2016_age65, svi_2018_age65, svi_2020_age65)
svi_dat_age65 <- list_svi_age65 %>% reduce(full_join, by = "fips")
write.csv(svi_dat_age65, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_age65.csv")

list_svi_disabl <- list(svi_2014_disabl, 
                       svi_2016_disabl, svi_2018_disabl, svi_2020_disabl)
svi_dat_disabl <- list_svi_disabl %>% reduce(full_join, by = "fips")
write.csv(svi_dat_disabl, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_disabl.csv")

list_svi_sngpnt <- list(svi_2010_sngprnt, svi_2014_sngpnt, 
                       svi_2016_sngpnt, svi_2018_sngpnt, svi_2020_sngpnt)
svi_dat_sngpnt <- list_svi_sngpnt %>% reduce(full_join, by = "fips")
write.csv(svi_dat_sngpnt, file = "./latest_missing_data_fixed/raw-data/svi-data/svi_hhd_var/svi_dat_sngpnt.csv")
