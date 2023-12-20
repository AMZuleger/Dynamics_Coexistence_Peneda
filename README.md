# Ecological Dynamics and Coexistence Patterns of Wild and Domestic Mammals in an Abandoned Landscape

Annika M. Zuleger, Andrea Perino, Henrique M. Pereira

## Abstract

The issue of agricultural land abandonment in Southern Europe has raised concerns about its impact on biodiversity. While abandoned areas can lead to positive developments like creating new habitats and restoring native vegetation, they can also result in human-wildlife conflicts, particularly concerning low-intensity farming with free-ranging domestic herbivores. Limited knowledge exists on livestock-wild ungulate interactions in abandoned European landscapes, making it essential to study the spatial and temporal distributional patterns of organisms to understand habitat selection and species distributions. To assess these impacts, a long-term biodiversity monitoring project was conducted using camera traps in the Peneda-Gerês National Park in Northern Portugal, focusing on exploring habitat preferences, dynamics of occupancy, and potential spatial interactions between domestic and wild species in the area using dynamic occupancy models. The study found that most species' occupancy remained stable across years, except for domestic cattle and horses, which declined slightly, while the Iberian ibex showed signs of repopulating the area. Initial occupancy as well as colonization probabilities were mostly unaffected by covariates for most of the study species, but we found species specific effects of habitat variables on extinction probabilities. Moreover, camera sensitivity improved detection probability of all species. Different habitat and weather factors also influenced species' detectability, highlighting the importance of accounting for those variables when modelling the detection process. We found no spatial interactions between domestic and wild herbivores, but we observed signs of shared area use between domestic cattle and wolves as well as between roe deer and wild boar and potential avoidance of roe deer and wolves in the area. Further investigation is needed to assess the interaction between livestock and wild ungulates in terms of resource competition and ecological carrying capacity.

## Structure
* R_Ecological_Dynamics_and_Coexistence_Patterns.R --> R Script to perform all analysis from the publication

### Data

* Presence-absence tables for each species (per grid cell and week) from 2015 to 2022 (e.g. Domestic cattle.csv)
* SiteCovs_2015_grid.csv --> intial site covariates for the first primary sampling period
* yearlySiteCovs_grid.rds --> yearly site covariates for each primary sampling period
* ObsCovs.rds --> observation covariates for each secondary sampling period

### Results

Created through R-Script
* Occupancy estimates per species for each grid cell per primary sampling period obtained from final models (e.g. Occ_cattle_site.csv)
* occ_all_grid.csv --> Occupancy estimates for all species per primary sampling period and grid cell obtained from final models
* occ_all_site_long.csv --> long format version of occ_all_grid.csv 
* occ_pred_all.csv --> Occupancy estimates and standard errors for each species per primary sampling period obtained from final models through bootstrapping
* var_effect.csv --> effect size of variables obtained from final models

### Figures

Created through R-Script
* Figure_3_Occupancy_trends.png --> Occupancy trends of the study species across years based on the final model obtained. A) Smoothed occupancy estimates and standard errors for the entire study area obtained from non-parametric bootstrapping with 1,000 iterations. B) Occupancy trends based on linear regression fitted to site-specific occupancy estimates.
* Figure_4_Variable_effects.png --> Effect sizes of the variables influencing initial Occupancy, Colonization probability, extinction probability and detection probability according to the best model following the automated QAIC selection procedure.
* Figure_4_Variable_effects_new.png --> manually updated version of Figure_4_Variable_effects.png
* Figure_5_Interaction_plot.png --> Effect of the occurrence of A) domestic cattle, B) domestic horses and C) wolf on the occupancy of the study species obtained from a linear regression based on occupancy estimates per sampling location and year obtained from the final models. 


## Authors

Annika M. Zuleger*, Andrea Perino and Henrique M. Pereira

*Corresponding author:
Annika Mikaela Zuleger
German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena- Leipzig
Puschstraße 4, 04103 Leipzig, Germany
Email: annika_mikaela.zuleger@idiv.de
