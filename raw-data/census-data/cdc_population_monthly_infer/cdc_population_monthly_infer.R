library(readr)
library(tidyverse)
library(dplyr)

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
  read_csv() %>% filter(month == 6) %>% select(-month) 
# pct_hispanic
pct_hispanic = sex_ethnicity %>% 
  group_by(year, fips, ethnicity) %>% summarise(pop = sum(pop)) #%>%
  # group_by(fips, year) %>%
  # mutate(total_pop = sum(pop)) %>% ungroup() %>%
  # filter(ethnicity == "Hispanic or Latino") %>%
  # mutate(pct_hispanic = pop / total_pop) %>% drop_na() %>%
  # pivot_wider(names_from = year, 
  #             values_from = ethnicity,
  #             names_prefix = "pct_hispanic_")
