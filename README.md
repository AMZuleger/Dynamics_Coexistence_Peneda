# Ecological Dynamics and Coexistence Patterns of Wild and Domestic Mammals in an Abandoned Landscape

Annika M. Zuleger, Andrea Perino, Henrique M. Pereira

## Abstract

The issue of agricultural land abandonment in Southern Europe has raised concerns about its impact on biodiversity. While abandoned areas can lead to positive developments like creating new habitats and restoring native vegetation, they can also result in human-wildlife conflicts, particularly in areas with extensive farming and free-ranging livestock. To understand habitat selection and use of livestock and wild ungulates, it is essential to study their spatial and temporal distribution patterns. In this context, we conducted a long-term large mammal monitoring project using camera traps in the Peneda-Gerês National Park in Northern Portugal. Our primary focus was on exploring habitat preferences, occupancy dynamics, and potential spatial use correlations between domestic and wild species, utilizing dynamic occupancy models. Most wild species exhibited stable area use patterns, while domestic species experienced marginal declines, and the Iberian ibex displayed signs of repopulation. We observed distinct effects of habitat variables on occupancy, colonization, and extinction, revealing species-specific patterns of habitat utilization. Human disturbance had a notable impact on domestic species but did not affect wild ones. Camera sensitivity emerged as a critical factor, enhancing detection probability for all species. Additionally, habitat and weather variables exerted varying effects on detection probabilities, underscoring the necessity of accounting for these factors in modeling the detection process. We found shared habitat preferences between cattle and horses, both positively correlated with wolves, suggesting potential human-wildlife conflicts. Despite extensive spatial overlap, domestic and wild species seem to exhibit ecological independence due to distinct strategies and low predation pressure. Overall, the study emphasizes the multifaceted factors influencing habitat use. The observed species associations contribute to understanding ecological relationships and potential resource competition, emphasizing the importance of considering environmental variables for effective wildlife conservation and management.

## Structure
* R_Ecological_Dynamics_and_Coexistence_Patterns_rv.R --> R Script to perform all analysis from the publication
* This is the only file needed to perform the entire analysis. Code will automaticall download data from GitHub and produce all Results and Figures provided in the subfolders.
  

### Data

* Presence-absence tables for each species (per grid cell and week) from 2015 to 2022 (e.g. Domestic cattle.csv)
* SiteCovs_2015_grid.csv --> intial site covariates for the first primary sampling period
* yearlySiteCovs_grid.rds --> yearly site covariates for each primary sampling period
* ObsCovs.rds --> observation covariates for each secondary sampling period

### Results (produced through R script)

* Occupancy estimates per species for each grid cell per primary sampling period obtained from final models (e.g. Occ_cattle_site.csv)
* occ_all_grid.csv --> Occupancy estimates for all species per primary sampling period and grid cell obtained from final models
* occ_all_site_long.csv --> long format version of occ_all_grid.csv 
* occ_pred_all.csv --> Occupancy estimates and standard errors for each species per primary sampling period obtained from final models through bootstrapping
* var_effect.csv --> effect size of variables obtained from final models
  

### Figures (produced thorugh R script)

* Figure_3_Variable_effects.png --> Effect sizes of the variables influencing initial Occupancy, Colonization probability, extinction probability and detection probability according to the best model after backwards elimination process (p > 0.157).
* Figure_4_Occupancy_trends.png --> Occupancy trends of the study species across years based on the final model obtained. A) Smoothed occupancy estimates and standard errors for the entire study area obtained from non-parametric bootstrapping with 1,000 iterations. B) Occupancy trends based on linear regression fitted to site-specific occupancy estimates.
* Figure_5_Interaction_plot.png --> Effect of the occurrence of A) domestic cattle, B) domestic horses and C) wolf on the occupancy of the study species obtained from a linear regression based on occupancy estimates per sampling location and year obtained from the final models. 


## Authors

Annika M. Zuleger*, Andrea Perino and Henrique M. Pereira

*Corresponding author:
Annika Mikaela Zuleger
German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena- Leipzig
Puschstraße 4, 04103 Leipzig, Germany
Email: annika_mikaela.zuleger@idiv.de
