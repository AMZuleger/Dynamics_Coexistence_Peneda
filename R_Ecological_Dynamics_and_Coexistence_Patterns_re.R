################################################################################
################ Ecological Dynamics and Coexistence Patterns ##################
#####of Wild and Domestic Mammals in an Abandoned Landscape ####################
################################################################################

############# Zuleger, A.M.; Perino, A.; Pereira, H.M. (2024) ##################

########################## Load required data ##################################

set.seed(123)

library(tidyr)
library(psych)
library(unmarked)
library(ggplot2)
library(reshape)
library(MASS)
library(gridExtra)

## Load required data from GitHub

t <- tempdir()
setwd(t) # set temporal wd

# download data from github
url <- "https://github.com/AMZuleger/Dynamics_Coexistence_Peneda/archive/refs/heads/main.zip"
download.file(url,destfile="Dynamics_Coexistence_Peneda-main.zip")
unzip(zipfile="Dynamics_Coexistence_Peneda-main.zip")

setwd(dir = file.path(t,"Dynamics_Coexistence_Peneda-main/Data"))

## Site covariates for season 1 (2015)

SiteCovs_2015 <- read.csv("SiteCovs_2015_grid.csv",header=T,sep=";")

## Yearly site covariates

yearlySiteCovs <- readRDS("yearlySiteCovs_grid.rds")

## Observation covariates

ObsCovs <- readRDS("ObsCovs.rds")

## Presence-absence tables per species

Cattle <- read.csv("Domestic cattle.csv",header=T,sep=";")
Horse <- read.csv("Domestic horse.csv",header=T,sep=";")
Deer <- read.csv("European roe deer.csv",header=T,sep=";")
Boar <- read.csv("Wild boar.csv",header=T,sep=";")
Fox <- read.csv("Red fox.csv",header=T,sep=";")
Wolf <- read.csv("Gray wolf.csv",header=T,sep=";")
Ibex <- read.csv("Iberian ibex.csv",header=T,sep=";")

################ Check for correlation between covariates ######################

for (i in 1:length(yearlySiteCovs)) {
  data <- gather(yearlySiteCovs[[i]])[2]
  names(data) <- names(yearlySiteCovs[i])
  if (i==1) {site_covs <- data} else {site_covs <- cbind(site_covs,data)}
}

for (i in 1:length(ObsCovs)) {
  data <- gather(ObsCovs[[i]])[2]
  names(data) <- names(ObsCovs[i])
  if (i==1) {obs_covs <- data} else {obs_covs <- cbind(obs_covs,data)}
}

corPlot(site_covs[c(1:3,5:6,8,11,14:15)],method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Site Covariates",xlas=2)
# removed CostDist_urban and verbatimElevation
corPlot(site_covs[c(1:2,6,8,11,14:15)],method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Site Covariates",xlas=2)

corPlot(obs_covs[1:5],method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Observation Covariates",xlas=2)
# keep Temp_min and Rad_max
corPlot(obs_covs[, c(1,4)],method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Observation Covariates",xlas=2)

############# Function for backwards variable selection  #######################

# Start with global model and reduce variables per parameter, starting with detection, then occupancy, colonization and extinction

parameter <- c("psi","gamma","epsilon","p")
parm <- c("psi","col","ext","det")

# Function for backward elimination for a specific parameter
backward_eliminate_parameter <- function(model, pn, threshold = 0.157) {
  formulas <- list(
    psi = model@psiformula,
    gamma = model@gamformula,
    epsilon = model@epsformula,
    p = model@detformula
  )
  formula <- formulas[[parameter[pn]]]
  while (TRUE) {
    model_summary <- summary(model)
    pvalues <- model_summary[[parm[pn]]][,4]
    names(pvalues) <- rownames(model_summary[[parm[pn]]])
    max_pvalue <- max(pvalues[-1], na.rm = TRUE)
    
    if (max_pvalue > threshold) {
      if (length(attr(terms(formula), "term.labels")) > 1) {
        # Update the formula by removing the predictor
        removed_predictor <- names(pvalues)[which(pvalues==max_pvalue)]
        updated_formula <- update(formula, paste("~ . -", removed_predictor))
      } else {
        # If only the intercept is left, set formula to ~1
        updated_formula <- as.formula("~1")
      }
      cat("Removing", removed_predictor, "from", parameter[pn], "formula\n")
      print(updated_formula)
      
      psiform <- if (parameter[pn] == "psi") updated_formula else formulas$psi
      gammaform <- if (parameter[pn] == "gamma") updated_formula else formulas$gamma
      epsilonform <- if (parameter[pn] == "epsilon") updated_formula else formulas$epsilon
      pform <- if (parameter[pn] == "p") updated_formula else formulas$p
      
      model <- colext(
        psiformula = psiform,
        gammaformula = gammaform,
        epsilonformula = epsilonform,
        pformula = pform,
        data = model@data
      )
      formulas[[parameter[pn]]] <- updated_formula
      formula <- updated_formula
    } else {
      cat("No more variables to remove\n")
      break
    }
  }
  return(list(model = model, formulas = formulas))
}

# Function to perform backwards elimination for each parameter 
iterate_backward_elimination <- function(model, threshold = 0.157) {
  repeat {
    detection <- backward_eliminate_parameter(model, 4, threshold)  # Passing index directly
    refined_model <- detection$model
    
    initial_psi <- backward_eliminate_parameter(refined_model, 1, threshold)  # Passing index directly
    refined_model <- initial_psi$model
    
    colonization <- backward_eliminate_parameter(refined_model, 2, threshold)  # Passing index directly
    refined_model <- colonization$model
    
    extinction <- backward_eliminate_parameter(refined_model, 3, threshold)  # Passing index directly
    refined_model <- extinction$model
    
    # Check if any parameters need re-evaluation, excluding the intercept
    new_summary <- summary(refined_model)
    psi_pvalues <- new_summary$psi[, 4]
    psi_pvalues <- psi_pvalues[names(psi_pvalues) != "(Intercept)"]
    gamma_pvalues <- new_summary$col[, 4]
    gamma_pvalues <- gamma_pvalues[names(gamma_pvalues) != "(Intercept)"]
    epsilon_pvalues <- new_summary$ext[, 4]
    epsilon_pvalues <- epsilon_pvalues[names(epsilon_pvalues) != "(Intercept)"]
    p_pvalues <- new_summary$det[, 4]
    p_pvalues <- p_pvalues[names(p_pvalues) != "(Intercept)"]
    
    # If all non-intercept p-values are below the threshold, break the loop
    if (all(c(psi_pvalues, gamma_pvalues, epsilon_pvalues, p_pvalues) <= threshold)) {
      break
    }
    
    model <- refined_model
  }
  return(refined_model)
}

################################################################################
######################### Analysis per species #################################
################################################################################

########################### Domestic Cattle ####################################

umfm_cattle <- unmarkedMultFrame(y = `Cattle`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_cattle) 

fm_cattle_global <- colext(~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                           ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                           ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                           ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri), umfm_cattle)
summary(fm_cattle_global)

# Backwards selection
final_model_cattle <- iterate_backward_elimination(fm_cattle_global, threshold = 0.157)
summary(final_model_cattle)

final_model_cattle <- colext(
  psiformula = final_model_cattle@psiformula,
  gammaformula = final_model_cattle@gamformula,
  epsilonformula = final_model_cattle@epsformula,
  pformula = final_model_cattle@detformula,
  data = umfm_cattle)

summary(final_model_cattle)

# Get variable estimates
variables_cattle <- data.frame(matrix(ncol = 10, nrow = 9))
colnames(variables_cattle) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_cattle$Species <- "Domestic cattle"
variables_cattle$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","OakForest","PineForest","Shrub","Urban_agri")

variables_cattle[1,3] <- backTransform(final_model_cattle,type="psi")@estimate
variables_cattle[1,4] <- SE(backTransform(final_model_cattle,type="psi"))

variables_cattle[1,5] <- backTransform(linearComb(final_model_cattle,c(1,0),type="col"))@estimate
variables_cattle[1,6] <- SE(backTransform(linearComb(final_model_cattle,c(1,0),type="col")))
variables_cattle[6,5] <- backTransform(linearComb(final_model_cattle,c(1,1),type="col"))@estimate
variables_cattle[6,6] <- SE(backTransform(linearComb(final_model_cattle,c(1,1),type="col")))

variables_cattle[1,7] <- backTransform(linearComb(final_model_cattle,c(1,0,0,0),type="ext"))@estimate
variables_cattle[1,8] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,0,0),type="ext")))
variables_cattle[5,7] <- backTransform(linearComb(final_model_cattle,c(1,1,0,0),type="ext"))@estimate
variables_cattle[5,8] <- SE(backTransform(linearComb(final_model_cattle,c(1,1,0,0),type="ext")))
variables_cattle[6,7] <- backTransform(linearComb(final_model_cattle,c(1,0,1,0),type="ext"))@estimate
variables_cattle[6,8] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,1,0),type="ext")))
variables_cattle[8,7] <- backTransform(linearComb(final_model_cattle,c(1,0,0,1),type="ext"))@estimate
variables_cattle[8,8] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,0,1),type="ext")))

variables_cattle[1,9] <- backTransform(linearComb(final_model_cattle,c(1,0,0,0,0,0),type="det"))@estimate
variables_cattle[1,10] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,0,0,0,0),type="det")))
variables_cattle[2,9] <- backTransform(linearComb(final_model_cattle,c(1,1,0,0,0,0),type="det"))@estimate
variables_cattle[2,10] <- SE(backTransform(linearComb(final_model_cattle,c(1,1,0,0,0,0),type="det")))
variables_cattle[3,9] <- backTransform(linearComb(final_model_cattle,c(1,0,1,0,0,0),type="det"))@estimate
variables_cattle[3,10] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,1,0,0,0),type="det")))
variables_cattle[7,9] <- backTransform(linearComb(final_model_cattle,c(1,0,0,1,0,0),type="det"))@estimate
variables_cattle[7,10] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,0,1,0,0),type="det")))
variables_cattle[8,9] <- backTransform(linearComb(final_model_cattle,c(1,0,0,0,1,0),type="det"))@estimate
variables_cattle[8,10] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,0,0,1,0),type="det")))
variables_cattle[9,9] <- backTransform(linearComb(final_model_cattle,c(1,0,0,0,0,1),type="det"))@estimate
variables_cattle[9,10] <- SE(backTransform(linearComb(final_model_cattle,c(1,0,0,0,0,1),type="det")))

# Get yearly occupancy estimates
m3 <- nonparboot(final_model_cattle,B = 1000)
occ_pred_cattle<- data.frame(year = c(2015:2022),
                             smoothed_occ = smoothed(final_model_cattle)[2,],
                             SE = m3@smoothed.mean.bsse[2,])
occ_pred_cattle$Species <- 'Domestic cattle'

cattle_year <- lm(smoothed_occ ~ year, occ_pred_cattle)
summary(cattle_year)

# Get site-specific occupancy estimates obtained from model output
occ_cattle_site <- data.frame(t(final_model_cattle@smoothed[2, , ]))
colnames(occ_cattle_site) <- unique(yearlySiteCovs(umfm_cattle)$year)
occ_cattle_site$grid <- unique(Cattle$grid) 

#write.table(occ_cattle_site, "Occ_cattle_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_cattle_site_long <- data.frame(melt(occ_cattle_site, id.vars="grid"))
names(occ_cattle_site_long) <- c("Grid","Year","Occupancy")
occ_cattle_site_long$Year <- as.numeric(as.character(occ_cattle_site_long$Year))
occ_cattle_site_long$Species <- "Domestic cattle"

########################### Domestic Horse #####################################

umfm_horse <- unmarkedMultFrame(y = `Horse`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_horse) 

fm_horse_global <- colext(~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                          ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                          ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                          ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),umfm_horse)

summary(fm_horse_global)

# Backwards elimination
final_model_horse <- iterate_backward_elimination(fm_horse_global, threshold = 0.157)
summary(final_model_horse)

final_model_horse <- colext(
  psiformula = final_model_horse@psiformula,
  gammaformula = final_model_horse@gamformula,
  epsilonformula = final_model_horse@epsformula,
  pformula = final_model_horse@detformula,
  data = umfm_horse)

summary(final_model_horse)

# Get variable estimates
variables_horse <- data.frame(matrix(ncol = 10, nrow = 9))
colnames(variables_horse) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_horse$Species <- "Domestic horse"
variables_horse$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","OakForest","PineForest","Shrub","Urban_agri")

variables_horse[1,3] <- backTransform(final_model_horse,type="psi")@estimate
variables_horse[1,4] <- SE(backTransform(final_model_horse,type="psi"))

variables_horse[1,5] <- backTransform(final_model_horse,type="col")@estimate
variables_horse[1,6] <- SE(backTransform(final_model_horse,type="col"))

variables_horse[1,7] <- backTransform(linearComb(final_model_horse,c(1,0,0),type="ext"))@estimate
variables_horse[1,8] <- SE(backTransform(linearComb(final_model_horse,c(1,0,0),type="ext")))
variables_horse[6,7] <- backTransform(linearComb(final_model_horse,c(1,1,0),type="ext"))@estimate
variables_horse[6,8] <- SE(backTransform(linearComb(final_model_horse,c(1,1,0),type="ext")))
variables_horse[9,7] <- backTransform(linearComb(final_model_horse,c(1,0,1),type="ext"))@estimate
variables_horse[9,8] <- SE(backTransform(linearComb(final_model_horse,c(1,0,1),type="ext")))

variables_horse[1,9] <- backTransform(linearComb(final_model_horse,c(1,0,0,0,0,0,0,0),type="det"))@estimate
variables_horse[1,10] <- SE(backTransform(linearComb(final_model_horse,c(1,0,0,0,0,0,0,0),type="det")))
variables_horse[2,9] <- backTransform(linearComb(final_model_horse,c(1,1,0,0,0,0,0,0),type="det"))@estimate
variables_horse[2,10] <- SE(backTransform(linearComb(final_model_horse,c(1,1,0,0,0,0,0,0),type="det")))
variables_horse[3,9] <- backTransform(linearComb(final_model_horse,c(1,0,1,0,0,0,0,0),type="det"))@estimate
variables_horse[3,10] <- SE(backTransform(linearComb(final_model_horse,c(1,0,1,0,0,0,0,0),type="det")))
variables_horse[4,9] <- backTransform(linearComb(final_model_horse,c(1,0,0,1,0,0,0,0),type="det"))@estimate
variables_horse[4,10] <- SE(backTransform(linearComb(final_model_horse,c(1,0,0,1,0,0,0,0),type="det")))
variables_horse[6,9] <- backTransform(linearComb(final_model_horse,c(1,0,0,0,1,0,0,0),type="det"))@estimate
variables_horse[6,10] <- SE(backTransform(linearComb(final_model_horse,c(1,0,0,0,1,0,0,0),type="det")))
variables_horse[7,9] <- backTransform(linearComb(final_model_horse,c(1,0,0,0,0,1,0,0),type="det"))@estimate
variables_horse[7,10] <- SE(backTransform(linearComb(final_model_horse,c(1,0,0,0,0,1,0,0),type="det")))
variables_horse[8,9] <- backTransform(linearComb(final_model_horse,c(1,0,0,0,0,0,1,0),type="det"))@estimate
variables_horse[8,10] <- SE(backTransform(linearComb(final_model_horse,c(1,0,0,0,0,0,1,0),type="det")))
variables_horse[9,9] <- backTransform(linearComb(final_model_horse,c(1,0,0,0,0,0,0,1),type="det"))@estimate
variables_horse[9,10] <- SE(backTransform(linearComb(final_model_horse,c(1,0,0,0,0,0,0,1),type="det")))

# Get yearly occupancy estimates
m3 <- nonparboot(final_model_horse,B = 1000)
occ_pred_horse<- data.frame(year = c(2015:2022),
                            smoothed_occ = smoothed(final_model_horse)[2,],
                            SE = m3@smoothed.mean.bsse[2,])
occ_pred_horse$Species <- 'Domestic horse'

horse_year <- lm(smoothed_occ ~ year, occ_pred_horse)
summary(horse_year)

# Get site-specific occupancy estimates obtained from model output
occ_horse_site <- data.frame(t(final_model_horse@smoothed[2, , ]))
colnames(occ_horse_site) <- unique(yearlySiteCovs(umfm_horse)$year)
occ_horse_site$grid <- unique(Horse$grid) 

#write.table(occ_horse_site, "Occ_horse_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_horse_site_long <- data.frame(melt(occ_horse_site, id.vars="grid"))
names(occ_horse_site_long) <- c("Grid","Year","Occupancy")
occ_horse_site_long$Year <- as.numeric(as.character(occ_horse_site_long$Year))
occ_horse_site_long$Species <- "Domestic horse"


############################# Roe deer #########################################

umfm_deer <- unmarkedMultFrame(y = `Deer`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_deer) 

fm_deer_global <- colext(~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri), umfm_deer)

summary(fm_deer_global)

# Backwards elimination
final_model_deer <- iterate_backward_elimination(fm_deer_global, threshold = 0.157)
summary(final_model_deer)

final_model_deer <- colext(
  psiformula = final_model_deer@psiformula,
  gammaformula = final_model_deer@gamformula,
  epsilonformula = final_model_deer@epsformula,
  pformula = final_model_deer@detformula,
  data = umfm_deer)

summary(final_model_deer)

# Get variable estimates
variables_deer <- data.frame(matrix(ncol = 10, nrow = 9))
colnames(variables_deer) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_deer$Species <- "Roe deer"
variables_deer$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","OakForest","PineForest","Shrub","Urban_agri")

variables_deer[1,3] <- backTransform(final_model_deer,type="psi")@estimate
variables_deer[1,4] <- SE(backTransform(final_model_deer,type="psi"))

variables_deer[1,5] <- backTransform(final_model_deer,type="col")@estimate
variables_deer[1,6] <- SE(backTransform(final_model_deer,type="col"))

variables_deer[1,7] <- backTransform(linearComb(final_model_deer,c(1,0),type="ext"))@estimate
variables_deer[1,8] <- SE(backTransform(linearComb(final_model_deer,c(1,0),type="ext")))
variables_deer[8,7] <- backTransform(linearComb(final_model_deer,c(1,1),type="ext"))@estimate
variables_deer[8,8] <- SE(backTransform(linearComb(final_model_deer,c(1,1),type="ext")))

variables_deer[1,9] <- backTransform(linearComb(final_model_deer,c(1,0,0,0,0,0,0),type="det"))@estimate
variables_deer[1,10] <- SE(backTransform(linearComb(final_model_deer,c(1,0,0,0,0,0,0),type="det")))
variables_deer[2,9] <- backTransform(linearComb(final_model_deer,c(1,1,0,0,0,0,0),type="det"))@estimate
variables_deer[2,10] <- SE(backTransform(linearComb(final_model_deer,c(1,1,0,0,0,0,0),type="det")))
variables_deer[3,9] <- backTransform(linearComb(final_model_deer,c(1,0,1,0,0,0,0),type="det"))@estimate
variables_deer[3,10] <- SE(backTransform(linearComb(final_model_deer,c(1,0,1,0,0,0,0),type="det")))
variables_deer[4,9] <- backTransform(linearComb(final_model_deer,c(1,0,0,1,0,0,0),type="det"))@estimate
variables_deer[4,10] <- SE(backTransform(linearComb(final_model_deer,c(1,0,0,1,0,0,0),type="det")))
variables_deer[6,9] <- backTransform(linearComb(final_model_deer,c(1,0,0,0,1,0,0),type="det"))@estimate
variables_deer[6,10] <- SE(backTransform(linearComb(final_model_deer,c(1,0,0,0,1,0,0),type="det")))
variables_deer[8,9] <- backTransform(linearComb(final_model_deer,c(1,0,0,0,0,1,0),type="det"))@estimate
variables_deer[8,10] <- SE(backTransform(linearComb(final_model_deer,c(1,0,0,0,0,1,0),type="det")))
variables_deer[9,9] <- backTransform(linearComb(final_model_deer,c(1,0,0,0,0,0,1),type="det"))@estimate
variables_deer[9,10] <- SE(backTransform(linearComb(final_model_deer,c(1,0,0,0,0,0,1),type="det")))

# Get yearly occupancy estimates
m3 <- nonparboot(final_model_deer,B = 1000)
occ_pred_deer<- data.frame(year = c(2015:2022),
                           smoothed_occ = smoothed(final_model_deer)[2,],
                           SE = m3@smoothed.mean.bsse[2,])
occ_pred_deer$Species <- 'Roe deer'

deer_year <- lm(smoothed_occ ~ year, occ_pred_deer)
summary(deer_year)

# Get site-specific occupancy estimates obtained from model output
occ_deer_site <- data.frame(t(final_model_deer@smoothed[2, , ]))
colnames(occ_deer_site) <- unique(yearlySiteCovs(umfm_deer)$year)
occ_deer_site$grid <- unique(Deer$grid) 

#write.table(occ_deer_site, "Occ_deer_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_deer_site_long <- data.frame(melt(occ_deer_site, id.vars="grid"))
names(occ_deer_site_long) <- c("Grid","Year","Occupancy")
occ_deer_site_long$Year <- as.numeric(as.character(occ_deer_site_long$Year))
occ_deer_site_long$Species <- "Roe deer"


############################# Wild boar ########################################

umfm_boar <- unmarkedMultFrame(y = `Boar`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_boar) 

fm_boar_global <- colext(~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri), umfm_boar)

summary(fm_boar_global)

# Backwards elimination
final_model_boar <- iterate_backward_elimination(fm_boar_global, threshold = 0.157)
summary(final_model_boar)

final_model_boar <- colext(
  psiformula = final_model_boar@psiformula,
  gammaformula = final_model_boar@gamformula,
  epsilonformula = final_model_boar@epsformula,
  pformula = final_model_boar@detformula,
  data = umfm_boar)

summary(final_model_boar)

# Get variable estimates
variables_boar <- data.frame(matrix(ncol = 10, nrow = 9))
colnames(variables_boar) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_boar$Species <- "Wild boar"
variables_boar$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","OakForest","PineForest","Shrub","Urban_agri")

variables_boar[1,3] <- backTransform(linearComb(final_model_boar,c(1,0),type="psi"))@estimate
variables_boar[1,4] <- SE(backTransform(linearComb(final_model_boar,c(1,0),type="psi")))
variables_boar[8,3] <- backTransform(linearComb(final_model_boar,c(1,1),type="psi"))@estimate
variables_boar[8,4] <- SE(backTransform(linearComb(final_model_boar,c(1,1),type="psi")))

variables_boar[1,5] <- backTransform(final_model_boar,type="col")@estimate
variables_boar[1,6] <- SE(backTransform(final_model_boar,type="col"))

variables_boar[1,7] <- backTransform(final_model_boar,type="ext")@estimate
variables_boar[1,8] <- SE(backTransform(final_model_boar,type="ext"))

variables_boar[1,9] <- backTransform(linearComb(final_model_boar,c(1,0,0,0,0,0),type="det"))@estimate
variables_boar[1,10] <- SE(backTransform(linearComb(final_model_boar,c(1,0,0,0,0,0),type="det")))
variables_boar[2,9] <- backTransform(linearComb(final_model_boar,c(1,1,0,0,0,0),type="det"))@estimate
variables_boar[2,10] <- SE(backTransform(linearComb(final_model_boar,c(1,1,0,0,0,0),type="det")))
variables_boar[3,9] <- backTransform(linearComb(final_model_boar,c(1,0,1,0,0,0),type="det"))@estimate
variables_boar[3,10] <- SE(backTransform(linearComb(final_model_boar,c(1,0,1,0,0,0),type="det")))
variables_boar[6,9] <- backTransform(linearComb(final_model_boar,c(1,0,0,1,0,0),type="det"))@estimate
variables_boar[6,10] <- SE(backTransform(linearComb(final_model_boar,c(1,0,0,1,0,0),type="det")))
variables_boar[8,9] <- backTransform(linearComb(final_model_boar,c(1,0,0,0,1,0),type="det"))@estimate
variables_boar[8,10] <- SE(backTransform(linearComb(final_model_boar,c(1,0,0,0,1,0),type="det")))
variables_boar[9,9] <- backTransform(linearComb(final_model_boar,c(1,0,0,0,0,1),type="det"))@estimate
variables_boar[9,10] <- SE(backTransform(linearComb(final_model_boar,c(1,0,0,0,0,1),type="det")))

# Get yearly occupancy estimates
m3 <- nonparboot(final_model_boar,B = 1000)
occ_pred_boar<- data.frame(year = c(2015:2022),
                           smoothed_occ = smoothed(final_model_boar)[2,],
                           SE = m3@smoothed.mean.bsse[2,])
occ_pred_boar$Species <- 'Wild boar'

boar_year <- lm(smoothed_occ ~ year, occ_pred_boar)
summary(boar_year)

# Get site-specific occupancy estimates obtained from model output
occ_boar_site <- data.frame(t(final_model_boar@smoothed[2, , ]))
colnames(occ_boar_site) <- unique(yearlySiteCovs(umfm_boar)$year)
occ_boar_site$grid <- unique(Boar$grid) 

#write.table(occ_boar_site, "Occ_boar_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_boar_site_long <- data.frame(melt(occ_boar_site, id.vars="grid"))
names(occ_boar_site_long) <- c("Grid","Year","Occupancy")
occ_boar_site_long$Year <- as.numeric(as.character(occ_boar_site_long$Year))
occ_boar_site_long$Species <- "Wild boar"


############################## Red fox #########################################

umfm_fox <- unmarkedMultFrame(y = `Fox`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_fox) 

fm_fox_global <- colext(~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                        ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                        ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                        ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri), umfm_fox)

summary(fm_fox_global)

# Backwards elimination
final_model_fox <- iterate_backward_elimination(fm_fox_global, threshold = 0.157)
summary(final_model_fox)

final_model_fox <- colext(
  psiformula = final_model_fox@psiformula,
  gammaformula = final_model_fox@gamformula,
  epsilonformula = final_model_fox@epsformula,
  pformula = final_model_fox@detformula,
  data = umfm_fox)

summary(final_model_fox)

# Get variable estimates
variables_fox <- data.frame(matrix(ncol = 10, nrow = 9))
colnames(variables_fox) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_fox$Species <- "Red fox"
variables_fox$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","OakForest","PineForest","Shrub","Urban_agri")

variables_fox[1,3] <- backTransform(linearComb(final_model_fox,c(1,0,0,0),type="psi"))@estimate
variables_fox[1,4] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0,0),type="psi")))
variables_fox[5,3] <- backTransform(linearComb(final_model_fox,c(1,1,0,0),type="psi"))@estimate
variables_fox[5,4] <- SE(backTransform(linearComb(final_model_fox,c(1,1,0,0),type="psi")))
variables_fox[8,3] <- backTransform(linearComb(final_model_fox,c(1,0,1,0),type="psi"))@estimate
variables_fox[8,4] <- SE(backTransform(linearComb(final_model_fox,c(1,0,1,0),type="psi")))
variables_fox[9,3] <- backTransform(linearComb(final_model_fox,c(1,0,0,1),type="psi"))@estimate
variables_fox[9,4] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0,1),type="psi")))

variables_fox[1,5] <- backTransform(linearComb(final_model_fox,c(1,0,0),type="col"))@estimate
variables_fox[1,6] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0),type="col")))
variables_fox[6,5] <- backTransform(linearComb(final_model_fox,c(1,1,0),type="col"))@estimate
variables_fox[6,6] <- SE(backTransform(linearComb(final_model_fox,c(1,1,0),type="col")))
variables_fox[9,5] <- backTransform(linearComb(final_model_fox,c(1,0,1),type="col"))@estimate
variables_fox[9,6] <- SE(backTransform(linearComb(final_model_fox,c(1,0,1),type="col")))

variables_fox[1,7] <- backTransform(linearComb(final_model_fox,c(1,0,0),type="ext"))@estimate
variables_fox[1,8] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0),type="ext")))
variables_fox[6,7] <- backTransform(linearComb(final_model_fox,c(1,1,0),type="ext"))@estimate
variables_fox[6,8] <- SE(backTransform(linearComb(final_model_fox,c(1,1,0),type="ext")))
variables_fox[7,7] <- backTransform(linearComb(final_model_fox,c(1,0,1),type="ext"))@estimate
variables_fox[7,8] <- SE(backTransform(linearComb(final_model_fox,c(1,0,1),type="ext")))

variables_fox[1,9] <- backTransform(linearComb(final_model_fox,c(1,0,0,0,0,0),type="det"))@estimate
variables_fox[1,10] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0,0,0,0),type="det")))
variables_fox[2,9] <- backTransform(linearComb(final_model_fox,c(1,1,0,0,0,0),type="det"))@estimate
variables_fox[2,10] <- SE(backTransform(linearComb(final_model_fox,c(1,1,0,0,0,0),type="det")))
variables_fox[3,9] <- backTransform(linearComb(final_model_fox,c(1,0,1,0,0,0),type="det"))@estimate
variables_fox[3,10] <- SE(backTransform(linearComb(final_model_fox,c(1,0,1,0,0,0),type="det")))
variables_fox[4,9] <- backTransform(linearComb(final_model_fox,c(1,0,0,1,0,0),type="det"))@estimate
variables_fox[4,10] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0,1,0,0),type="det")))
variables_fox[6,9] <- backTransform(linearComb(final_model_fox,c(1,0,0,0,1,0),type="det"))@estimate
variables_fox[6,10] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0,0,1,0),type="det")))
variables_fox[8,9] <- backTransform(linearComb(final_model_fox,c(1,0,0,0,0,1),type="det"))@estimate
variables_fox[8,10] <- SE(backTransform(linearComb(final_model_fox,c(1,0,0,0,0,1),type="det")))

# Get yearly occupancy estimates
m3 <- nonparboot(final_model_fox,B = 1000)
occ_pred_fox<- data.frame(year = c(2015:2022),
                          smoothed_occ = smoothed(final_model_fox)[2,],
                          SE = m3@smoothed.mean.bsse[2,])
occ_pred_fox$Species <- 'Red fox'

fox_year <- lm(smoothed_occ ~ year, occ_pred_fox)
summary(fox_year)

# Get site-specific occupancy estimates obtained from model output
occ_fox_site <- data.frame(t(final_model_fox@smoothed[2, , ]))
colnames(occ_fox_site) <- unique(yearlySiteCovs(umfm_fox)$year)
occ_fox_site$grid <- unique(Fox$grid) 

#write.table(occ_fox_site, "Occ_fox_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_fox_site_long <- data.frame(melt(occ_fox_site, id.vars="grid"))
names(occ_fox_site_long) <- c("Grid","Year","Occupancy")
occ_fox_site_long$Year <- as.numeric(as.character(occ_fox_site_long$Year))
occ_fox_site_long$Species <- "Red fox"

############################## Gray wolf #######################################

umfm_wolf <- unmarkedMultFrame(y = `Wolf`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_wolf) 

fm_wolf_global <- colext(~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri), umfm_wolf)

summary(fm_wolf_global)

# Backwards elimination
final_model_wolf <- iterate_backward_elimination(fm_wolf_global, threshold = 0.157)
summary(final_model_wolf)

final_model_wolf <- colext(
  psiformula = final_model_wolf@psiformula,
  gammaformula = final_model_wolf@gamformula,
  epsilonformula = final_model_wolf@epsformula,
  pformula = final_model_wolf@detformula,
  data = umfm_wolf)

summary(final_model_wolf)

variables_wolf <- data.frame(matrix(ncol = 10, nrow = 9))
colnames(variables_wolf) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_wolf$Species <- "Gray wolf"
variables_wolf$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","OakForest","PineForest","Shrub","Urban_agri")

variables_wolf[1,3] <- backTransform(final_model_wolf,type="psi")@estimate
variables_wolf[1,4] <- SE(backTransform(final_model_wolf,type="psi"))

variables_wolf[1,5] <- backTransform(linearComb(final_model_wolf,c(1,0),type="col"))@estimate
variables_wolf[1,6] <- SE(backTransform(linearComb(final_model_wolf,c(1,0),type="col")))
variables_wolf[7,5] <- backTransform(linearComb(final_model_wolf,c(1,1),type="col"))@estimate
variables_wolf[7,6] <- SE(backTransform(linearComb(final_model_wolf,c(1,1),type="col")))

variables_wolf[1,7] <- backTransform(final_model_wolf,type="ext")@estimate
variables_wolf[1,8] <- SE(backTransform(final_model_wolf,type="ext"))

variables_wolf[1,9] <- backTransform(linearComb(final_model_wolf,c(1,0,0,0),type="det"))@estimate
variables_wolf[1,10] <- SE(backTransform(linearComb(final_model_wolf,c(1,0,0,0),type="det")))
variables_wolf[2,9] <- backTransform(linearComb(final_model_wolf,c(1,1,0,0),type="det"))@estimate
variables_wolf[2,10] <- SE(backTransform(linearComb(final_model_wolf,c(1,1,0,0),type="det")))
variables_wolf[3,9] <- backTransform(linearComb(final_model_wolf,c(1,0,1,0),type="det"))@estimate
variables_wolf[3,10] <- SE(backTransform(linearComb(final_model_wolf,c(1,0,1,0),type="det")))
variables_wolf[4,9] <- backTransform(linearComb(final_model_wolf,c(1,0,0,1),type="det"))@estimate
variables_wolf[4,10] <- SE(backTransform(linearComb(final_model_wolf,c(1,0,0,1),type="det")))

# Get yearly occupancy estimates
m3 <- nonparboot(final_model_wolf,B = 1000)
occ_pred_wolf<- data.frame(year = c(2015:2022),
                           smoothed_occ = smoothed(final_model_wolf)[2,],
                           SE = m3@smoothed.mean.bsse[2,])
occ_pred_wolf$Species <- 'Gray wolf'

wolf_year <- lm(smoothed_occ ~ year, occ_pred_wolf)
summary(wolf_year)

# Get site-specific occupancy estimates obtained from model output
occ_wolf_site <- data.frame(t(final_model_wolf@smoothed[2, , ]))
colnames(occ_wolf_site) <- unique(yearlySiteCovs(umfm_wolf)$year)
occ_wolf_site$grid <- unique(Wolf$grid) 

#write.table(occ_wolf_site, "Occ_wolf_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_wolf_site_long <- data.frame(melt(occ_wolf_site, id.vars="grid"))
names(occ_wolf_site_long) <- c("Grid","Year","Occupancy")
occ_wolf_site_long$Year <- as.numeric(as.character(occ_wolf_site_long$Year))
occ_wolf_site_long$Species <- "Gray wolf"

############################ Iberian ibex ######################################

umfm_ibex <- unmarkedMultFrame(y = `Ibex`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_ibex) 

fm_ibex_global <- colext(~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(CostDist_road)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri),
                         ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(PineForest)+scale(Shrub)+scale(Urban_agri), umfm_ibex)

summary(fm_ibex_global)

# Backwards elimination
final_model_ibex <- iterate_backward_elimination(fm_ibex_global, threshold = 0.157)
summary(final_model_ibex)

final_model_ibex <- colext(
  psiformula = final_model_ibex@psiformula,
  gammaformula = final_model_ibex@gamformula,
  epsilonformula = final_model_ibex@epsformula,
  pformula = final_model_ibex@detformula,
  data = umfm_ibex)

summary(final_model_ibex)

variables_ibex <- data.frame(matrix(ncol = 10, nrow = 9))
colnames(variables_ibex) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_ibex$Species <- "Iberian ibex"
variables_ibex$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","OakForest","PineForest","Shrub","Urban_agri")

variables_ibex[1,3] <- backTransform(final_model_ibex,type="psi")@estimate
variables_ibex[1,4] <- SE(backTransform(final_model_ibex,type="psi"))

variables_ibex[1,5] <- backTransform(final_model_ibex,type="col")@estimate
variables_ibex[1,6] <- SE(backTransform(final_model_ibex,type="col"))

variables_ibex[1,7] <- backTransform(final_model_ibex,type="ext")@estimate
variables_ibex[1,8] <- SE(backTransform(final_model_ibex,type="ext"))

variables_ibex[1,9] <- backTransform(linearComb(final_model_ibex,c(1,0,0,0,0),type="det"))@estimate
variables_ibex[1,10] <- SE(backTransform(linearComb(final_model_ibex,c(1,0,0,0,0),type="det")))
variables_ibex[2,9] <- backTransform(linearComb(final_model_ibex,c(1,1,0,0,0),type="det"))@estimate
variables_ibex[2,10] <- SE(backTransform(linearComb(final_model_ibex,c(1,1,0,0,0),type="det")))
variables_ibex[3,9] <- backTransform(linearComb(final_model_ibex,c(1,0,1,0,0),type="det"))@estimate
variables_ibex[3,10] <- SE(backTransform(linearComb(final_model_ibex,c(1,0,1,0,0),type="det")))
variables_ibex[6,9] <- backTransform(linearComb(final_model_ibex,c(1,0,0,1,0),type="det"))@estimate
variables_ibex[6,10] <- SE(backTransform(linearComb(final_model_ibex,c(1,0,0,1,0),type="det")))
variables_ibex[8,9] <- backTransform(linearComb(final_model_ibex,c(1,0,0,0,1),type="det"))@estimate
variables_ibex[8,10] <- SE(backTransform(linearComb(final_model_ibex,c(1,0,0,0,1),type="det")))

# Get yearly occupancy estimates
m3 <- nonparboot(final_model_ibex,B = 1000)
occ_pred_ibex<- data.frame(year = c(2015:2022),
                           smoothed_occ = smoothed(final_model_ibex)[2,],
                           SE = m3@smoothed.mean.bsse[2,])
occ_pred_ibex$Species <- 'Iberian ibex'

ibex_year <- lm(smoothed_occ ~ year, occ_pred_ibex)
summary(ibex_year)

# Get site-specific occupancy estimates obtained from model output
occ_ibex_site <- data.frame(t(final_model_ibex@smoothed[2, , ]))
colnames(occ_ibex_site) <- unique(yearlySiteCovs(umfm_ibex)$year)
occ_ibex_site$grid <- unique(Ibex$grid) 

#write.table(occ_ibex_site, "Occ_ibex_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_ibex_site_long <- data.frame(melt(occ_ibex_site, id.vars="grid"))
names(occ_ibex_site_long) <- c("Grid","Year","Occupancy")
occ_ibex_site_long$Year <- as.numeric(as.character(occ_ibex_site_long$Year))
occ_ibex_site_long$Species <- "Iberian ibex"


###################### Create tables for all species ########################### 

setwd("H:/Occupancy/Publication_data/Results_bs")

# Occupancy estimates for each species per year
occ_pred_all <- rbind(occ_pred_cattle,occ_pred_horse,occ_pred_deer,occ_pred_boar,occ_pred_fox,occ_pred_wolf,occ_pred_ibex)
occ_pred_all$Species <- as.factor(occ_pred_all$Species)
occ_pred_all$Species <- factor(occ_pred_all$Species, c("Domestic cattle", "Domestic horse", "Roe deer","Wild boar","Red fox","Gray wolf","Iberian ibex"))
write.table(occ_pred_all, "occ_pred_all.csv",sep=";",dec=".",row.names=F,col.names=T)

# Long format
occ_all_site_long <- rbind(occ_cattle_site_long,occ_horse_site_long,occ_deer_site_long,occ_boar_site_long,occ_fox_site_long,occ_wolf_site_long,occ_ibex_site_long)
occ_all_site_long$Species <- as.factor(occ_all_site_long$Species)
occ_all_site_long$Species <- factor(occ_all_site_long$Species, c("Domestic cattle", "Domestic horse", "Roe deer","Wild boar","Red fox","Gray wolf","Iberian ibex"))
write.table(occ_all_site_long, "occ_all_site_long.csv",sep=";",dec=".",row.names=F,col.names=T)

# Occupancy estimates for each species per grid cell and year
occ_all_grid <- cbind(occ_cattle_site_long[-4],occ_horse_site_long[3],occ_deer_site_long[3],occ_boar_site_long[3],occ_fox_site_long[3],occ_wolf_site_long[3],occ_ibex_site_long[3])
names(occ_all_grid) <- c("Grid","Year","Cattle","Horse","RoeDeer","WildBoar","Fox","Wolf","Ibex")
write.table(occ_all_grid, "occ_all_grid.csv",sep=";",dec=".",row.names=F,col.names=T)

# Effect of variables for each species
var_effect <- rbind(variables_cattle,variables_horse,variables_deer,variables_boar,variables_fox,variables_wolf,variables_ibex)
write.table(var_effect, "var_effect.csv",sep=";",dec=".",row.names=F,col.names=T)

####################### Interaction with other species ######################### 
########## Test for possible interactions with a linear regression #############

cattle_int<- lm(Cattle ~ Horse+RoeDeer+WildBoar+Fox+Wolf+Ibex, occ_all_grid)
summary(cattle_int)

####

horse_int <- lm(Horse ~ Cattle+RoeDeer+WildBoar+Fox+Wolf+Ibex, occ_all_grid)
summary(horse_int)

####

deer_int <- lm(RoeDeer ~ Cattle+Horse+WildBoar+Fox+Wolf+Ibex, occ_all_grid)
summary(deer_int)

####

boar_int <- lm(WildBoar ~ Cattle+Horse+RoeDeer+Fox+Wolf+Ibex, occ_all_grid)
summary(boar_int)

####

fox_int <- lm(Fox ~ Cattle+Horse+RoeDeer+WildBoar+Wolf+Ibex, occ_all_grid)
summary(fox_int)

####

wolf_int <- lm(Wolf ~ Cattle+Horse+RoeDeer+WildBoar+Fox+Ibex, occ_all_grid)
summary(wolf_int)

####

ibex_int <- lm(Ibex ~ Cattle+Horse+RoeDeer+WildBoar+Fox+Wolf, occ_all_grid)
summary(ibex_int)


################################################################################
########################## Create figures ######################################
################################################################################

get_only_legend <- function(plot) { 
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  legend <- plot_table$grobs[[legend_plot]] 
  return(legend) 
} 

############### Figure 3. Effect sizes of the variables  ####################### 

## Occupancy

occ_effect <- var_effect[1:4]
occ_effect <- na.omit(occ_effect)

occ_effect$var <- as.character(interaction(occ_effect$Species,occ_effect$Variable))
occ_effect$var <- factor(occ_effect$var, levels = occ_effect$var)
occ_effect$Variable <- as.factor(occ_effect$Variable)
occ_effect$Variable <- factor(occ_effect$Variable, c("Intercept", "CostDist_road","OakForest","PineForest","Shrub","Urban_agri"))

point_shape <- rep(c(15,16,17,0,1,2),7)
point_col <- c(rep("#1B9E77",6),rep("#D95F02",6),rep("#7570B3",6),rep("#A6761D",6),rep("#E6AB02",6),rep("#E7298A",6),rep("#66A61E",6))

occupancy_plot <- ggplot(aes(y = var),data=occ_effect) + 
  geom_point(aes(x=Occupancy,shape=Variable), size=5,color=point_col) +
  scale_shape_manual(values = point_shape) +
  geom_linerange(aes(xmin=Occupancy-Occ_SE, xmax=Occupancy+Occ_SE),color=point_col) +
  geom_hline(yintercept = c(6.5,12.5,18.5,24.5,30.5,36.5)) +
  geom_segment(aes(x = 1, y = 42.5, xend = 1, yend = 36.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 1, y = 36.5, xend = 1, yend = 30.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 1, y = 30.5, xend = 1, yend = 24.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 0.891, y = 24.5, xend = 0.891, yend = 18.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 0.637, y = 18.5, xend = 0.637, yend = 12.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 1, y = 12.5, xend = 1, yend = 6.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.000, y = 6.5, xend = 0.000, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev, labels=rev(occ_effect$Variable)) +
  xlim(-0.1,1) +
  ggtitle("Initial occupancy") +
  theme_bw() +
  theme(legend.position = c(0.11,0.87),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=16),
        plot.title = element_text(size=24,hjust=0.5,face="bold"),
        axis.line = element_line( linetype = "solid"))

## Colonization

col_effect <- var_effect[c(1,2,5,6)]
col_effect <- col_effect[which(col_effect$Variable != "Sensitivity" & col_effect$Variable != "Rad_max" & col_effect$Variable != "Temp_min"),]

col_effect$var <- as.character(interaction(col_effect$Species,col_effect$Variable))
col_effect$var <- factor(col_effect$var, levels = col_effect$var)
col_effect$Variable <- as.factor(col_effect$Variable)
col_effect$Variable <- factor(col_effect$Variable, c("Intercept", "CostDist_road","OakForest","PineForest","Shrub","Urban_agri"))

colonization_plot <- ggplot(aes(y = var),data=col_effect) + 
  geom_point(aes(x=Colonization,shape=Variable), size=5,color=point_col) +
  scale_shape_manual(values = point_shape) +
  geom_linerange(aes(xmin=Colonization-Col_SE, xmax=Colonization+Col_SE),color=point_col) +
  geom_hline(yintercept = c(6.5,12.5,18.5,24.5,30.5,36.5)) +
  geom_segment(aes(x = 0.515, y = 42.5, xend = 0.515, yend = 36.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 1, y = 36.5, xend = 1, yend = 30.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 0.481, y = 30.5, xend = 0.481, yend = 24.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 1, y = 24.5, xend = 1, yend = 18.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 1, y = 18.5, xend = 1, yend = 12.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 0.003, y = 12.5, xend = 0.003, yend = 6.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.001, y = 6.5, xend = 0.001, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev, labels=rev(col_effect$Variable)) +
  xlim(-0.1,1) +
  ggtitle("Colonization probability") +
  theme_bw() +
  theme(legend.position = c(0.11,0.87),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=16),
        plot.title = element_text(size=24,hjust=0.5,face="bold"),
        axis.line = element_line( linetype = "solid"))

## Extinction

ext_effect <- var_effect[c(1,2,7,8)]
ext_effect <- ext_effect[which(ext_effect$Variable != "Sensitivity" & ext_effect$Variable != "Rad_max" & ext_effect$Variable != "Temp_min"),]

ext_effect$var <- as.character(interaction(ext_effect$Species,ext_effect$Variable))
ext_effect$var <- factor(ext_effect$var, levels = ext_effect$var)
ext_effect$Variable <- as.factor(ext_effect$Variable)
ext_effect$Variable <- factor(occ_effect$Variable, c("Intercept", "CostDist_road","OakForest","PineForest","Shrub","Urban_agri"))

extinction_plot <- ggplot(aes(y = var),data=ext_effect) + 
  geom_point(aes(x=Extinction,shape=Variable), size=5,color=point_col) +
  scale_shape_manual(values = point_shape) +
  geom_linerange(aes(xmin=Extinction-Ext_SE, xmax=Extinction+Ext_SE),color=point_col) +
  geom_hline(yintercept = c(6.5,12.5,18.5,24.5,30.5,36.5)) +
  geom_segment(aes(x = 0.309, y = 42.5, xend = 0.309, yend = 36.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 0.124, y = 36.5, xend = 0.124, yend = 30.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 0.000, y = 30.5, xend = 0.000, yend = 24.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 0.015, y = 24.5, xend = 0.015, yend = 18.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 0.412, y = 18.5, xend = 0.412, yend = 12.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 0.000, y = 12.5, xend = 0.000, yend = 6.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.001, y = 6.5, xend = 0.001, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev, labels=rev(ext_effect$Variable)) +
  xlim(-0.1,1) +
  ggtitle("Extinction probability") +
  theme_bw() +
  theme(legend.position = c(0.89,0.87),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=16),
        plot.title = element_text(size=24,hjust=0.5,face="bold"),
        axis.line = element_line( linetype = "solid"))

## Detection

det_effect <- var_effect[c(1,2,9,10)]
det_effect <- det_effect[which(det_effect$Variable != "CostDist_road"),]

det_effect$var <- as.character(interaction(det_effect$Species,det_effect$Variable))
det_effect$var <- factor(det_effect$var, levels = det_effect$var)
det_effect$Variable <- factor(det_effect$Variable, c("Intercept", "Sensitivity", "Rad_max","Temp_min","OakForest","PineForest","Shrub","Urban_agri"))

point_shape_det <- rep(c(15,16,17,18,0,1,2,5),7)
point_col_det <- c(rep("#1B9E77",8),rep("#D95F02",8),rep("#7570B3",8),rep("#A6761D",8),rep("#E6AB02",8),rep("#E7298A",8),rep("#66A61E",8))

detection_plot <- ggplot(aes(y = var),data=det_effect) + 
  geom_point(aes(x=Detection,shape=Variable), size=5,color=point_col_det) +
  scale_shape_manual(values = point_shape_det) +
  geom_linerange(aes(xmin=Detection-Det_SE, xmax=Detection+Det_SE),color=point_col_det) +
  geom_hline(yintercept = c(8.5,16.5,24.5,32.5,40.5,48.5)) +
  geom_segment(aes(x = 0.277, y = 54.5, xend = 0.277, yend = 48.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 0.292, y = 48.5, xend = 0.292, yend = 40.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 0.388, y = 40.5, xend = 0.388, yend = 32.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 0.275, y = 32.5, xend = 0.275, yend = 24.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 0.111, y = 24.5, xend = 0.111, yend = 16.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 0.027, y = 16.5, xend = 0.027, yend = 8.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.044, y = 8.5, xend = 0.044, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev) +
  xlim(-0.1,1) +
  ggtitle("Detection probability") +
  theme_bw() +
  theme(legend.position = c(0.9,0.84),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=16),
        plot.title = element_text(size=24,hjust=0.5,face="bold"),
        axis.line = element_line( linetype = "solid"))

var_effect$Species <- as.factor(var_effect$Species)
var_effect$Species <- factor(var_effect$Species, c("Domestic cattle", "Domestic horse", "Roe deer","Wild boar","Red fox","Gray wolf","Iberian ibex"))

plot_legend <-ggplot(var_effect, aes(y=Occupancy, fill=Species)) +
  geom_bar() +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=20,face="bold"),
        legend.text = element_text(size=20))

legend <- get_only_legend(plot_legend) 

Figure_3 <- arrangeGrob(occupancy_plot, colonization_plot, extinction_plot, detection_plot, legend, layout_matrix = rbind(c(1, 2),c(3, 4),c(5,5)),heights = c(1/2, 1/2, 1/9))
ggsave(file="Figure_3_Variable_effects.png",Figure_3, width=20, height=28,dpi=300)


#################### Figure 4. Occupancy trends ################################

occ_all_mod_legend <- ggplot(data=occ_pred_all, aes(x=year, y=smoothed_occ)) +
  geom_ribbon(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE, fill=Species),alpha=0.3) +
  geom_line(aes(x=year, y=smoothed_occ, col=Species),size=1) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  theme(legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=24,face="bold"),
        legend.text = element_text(size=20))

legend <- get_only_legend(occ_all_mod_legend)

occ_all_mod <- ggplot(data=occ_pred_all, aes(x=year, y=smoothed_occ)) +
  geom_ribbon(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE, fill=Species),alpha=0.3) +
  geom_line(aes(x=year, y=smoothed_occ, col=Species,linetype=Species),size=1.5) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  scale_linetype_manual(values=c("dashed","twodash","dashed","dashed","dashed","dashed","solid")) +
  labs(title="Occupancy",x="Year", y="Smoothed occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=18), 
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=26,face="bold"),
        axis.title.x= element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1.2, 'cm'),
        legend.title = element_text(size=20,face="bold"),
        legend.position="none")

occ_all_lm <- ggplot(occ_all_site_long, aes(Year, Occupancy)) +
  geom_point(aes(x=Year, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Year, y=Occupancy, col=Species,fill=Species,linetype=Species),alpha=0.3,size=1.5,method="lm") +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  scale_linetype_manual(values=c("dashed","twodash","dashed","dashed","dashed","dashed","solid")) +
  labs(title="",x="Year", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=18), 
        axis.title=element_text(size=20,face="bold"),
        plot.title = element_blank(),
        legend.position="none")

Figure_4 <- arrangeGrob(occ_all_mod, occ_all_lm, legend, layout_matrix = rbind(c(1,3),c(2, 3)),widths=c(1/2,1/8))
ggsave(file="Figure_4_Occupancy_trends.png",Figure_4, width=16, height=20,dpi=300)


############### Figure 5. Spatial use correlations  ############################ 

cattle <- melt(occ_all_grid, id.vars=c("Grid", "Year","Cattle"))
names(cattle) <- c("Grid","Year","Cattle","Species","Occupancy")

cattle_int <- ggplot(cattle, aes(Cattle, Occupancy)) +
  geom_point(aes(x=Cattle, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Cattle, y=Occupancy, col=Species,fill=Species,linetype=Species),alpha=0.3,size=1,method="lm") +
  scale_fill_manual(values=c("#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  scale_linetype_manual(values=c("solid","dashed","dashed","dashed","solid","dashed")) +
  labs(title="Occupancy estimates",x="Cattle occupancy", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        plot.title = element_blank(),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=20,face="bold"),
        legend.position="none")

horse <- melt(occ_all_grid, id.vars=c("Grid", "Year","Horse"))
names(horse) <- c("Grid","Year","Horse","Species","Occupancy")

horse_int <- ggplot(horse, aes(Horse, Occupancy)) +
  geom_point(aes(x=Horse, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Horse, y=Occupancy, col=Species,fill=Species,linetype=Species),alpha=0.3,size=1,method="lm") +
  scale_fill_manual(values=c("#1B9E77","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  scale_linetype_manual(values=c("solid","dashed","dashed","solid","dashed","solid")) +
  labs(title="Occupancy estimates",x="Horse occupancy", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        plot.title = element_blank(),
        axis.title.y= element_blank(),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=20,face="bold"),
        legend.position="none")

wolf <- melt(occ_all_grid, id.vars=c("Grid", "Year","Wolf"))
names(wolf) <- c("Grid","Year","Wolf","Species","Occupancy")

wolf_int <- ggplot(wolf, aes(Wolf, Occupancy)) +
  geom_point(aes(x=Wolf, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Wolf, y=Occupancy, col=Species,fill=Species,linetype=Species),alpha=0.3,size=1,method="lm") +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#66A61E"))  +
  scale_linetype_manual(values=c("solid","dashed","twodash","dashed","dashed","solid")) +
  labs(title="Occupancy estimates",x="Wolf occupancy", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        plot.title = element_blank(),
        axis.title.y= element_blank(),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=20,face="bold"),
        legend.position="none")

legend <- get_only_legend(occ_all_mod_legend)

Figure_5 <- arrangeGrob(cattle_int, horse_int, wolf_int, legend, nrow=1,widths=c(1/3,1/3,1/3,1/5))
ggsave(file="Figure_5_Interaction_plot.png",Figure_5, width=20, height=7,dpi=300)


################################################################################
################################################################################

