# Ecological Dynamics and Coexistence Patterns of Wild and Domestic Mammals in an Abandoned Landscape

Annika M. Zuleger, Andrea Perino, Henrique M. Pereira

This repository contains all scripts and data used for analysis and figures for the paper **Zuleger, Annika M., Perino, Andrea, Pereira, Henrique M. (2024): _Ecological Dynamics and Coexistence Patterns of Wild and Domestic Mammals in an Abandoned Landscape._ Wildlife Biology. DOI: 10.1002/wlb3.01319**

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


## Authors

Annika M. Zuleger*, Andrea Perino and Henrique M. Pereira

*Corresponding author:
Annika Mikaela Zuleger
German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena- Leipzig
Puschstraße 4, 04103 Leipzig, Germany
Email: annika_mikaela.zuleger@idiv.de
