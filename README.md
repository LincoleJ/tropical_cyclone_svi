# Long-Term Social Vulnerability Changes After Tropical Cyclones in Contiguous United States

![image](https://github.com/user-attachments/assets/c9d95746-171c-456b-8445-feab8ab6484c)

This is the data repository for publicly available code and data to conduct analyses in the paper titled "Long-Term Effects of Tropical Cyclones on Social Vulnerability in the United States."

We use a synthetic control approach to analyze fourteen years of tropical cyclone exposure data across contiguous U.S., and provide evidence that tropical cyclone exposures first induces a spike in social vulnerability immediately after tropical cyclones, but over time, affected regions tend to recover and experience gentrification, leaning toward lower social vulnerability compared to similar regions not affected by tropical cyclones.

### Code:

1. [raw-data](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/raw-data) process downloaded geospatial and demographic data from various data source.
   
3. [processed-data](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/processed-data) puts them into a tabular for statistical analysis purpose.
   
4. [balancing](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/balancing) apply covariate balancing synthetic control approach to obtain control weights to create the synthetic control region.

5. [analysis](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/analysis) conduct outcome analysis on covariate balanced data to generate main results and result graphs.
   
7. [figures_tables](https://github.com/LincoleJ/tropical_cyclone_svi/tree/main/figures_tables) generate figures and tables in the main text and supplementary materials.
