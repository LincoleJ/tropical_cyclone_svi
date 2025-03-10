# Evolution of Long-Term Social Vulnerability After Tropical Cyclones in the United States

![image](https://img.shields.io/badge/Study_Status-Submitted-red)

This is the data repository for publicly available code and data to conduct analyses in the paper titled "Long-Term Effects of Tropical Cyclones on Social Vulnerability in the United States."

We use a synthetic control approach to analyze fourteen years of tropical cyclone exposure data across contiguous U.S., and provide evidence that tropical cyclone exposures first induces a spike in social vulnerability immediately after tropical cyclones, but over time, affected regions tend to recover and experience gentrification, and exhibit even lower social vulnerability in the long run compared to similar regions not affected by tropical cyclones.


### Code:

1. [raw-data](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/raw-data) process downloaded geospatial and demographic data from various data source.
   
2. [processed-data](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/processed-data) puts them into a tabular for statistical analysis purpose.
   
3. [balancing](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/balancing) apply covariate balancing synthetic control approach to obtain control weights to create the synthetic control region.

4. [analysis](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/analysis) conduct outcome analysis on covariate balanced data to generate main results and result graphs.
   
5. [figures_tables](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/figures_tables) generate figures and tables in the main text and supplementary materials.


### Data Source:

| Data | Sources |
| --- | --- |
| Tropical Cyclone Exposure | [HURDAT](https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html) | 
| Social Vulnerability Index | [CDC/ATSDR](https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html) |
| Metereological | [PRISM](https://prism.oregonstate.edu/) |
| Percentage Minority Population | [CDC Wonder](https://wonder.cdc.gov/) |
| Gender Ratio | [CDC Wonder](https://wonder.cdc.gov/) |
| Per Capita Income | [BEA](https://www.bea.gov/data/income-saving/personal-income-county-metro-and-other-areas) |
| Percentage Below Poverty Threshold | [USCB](https://www.census.gov/topics/income-poverty/poverty/data/datasets.html) |
| Percentage High School Graduate | [USDA](https://data.ers.usda.gov/reports.aspx?ID=17829)

### Data: 

All data needed to evaluate the conclusions in the paper are present in the paper and/or the Supplementary Information and Online Repository. Those interested in the original data can contact the corresponding author. 

All the analyses are run on local R.


### Contact us: 

* Email: xw2892@cumc.columbia.edu; lj2575@cumc.columbia.edu; yoshira@ucla.edu; rmp2198@columbia.edu
