####################################################################################################################
###### Ecological Dynamics and Coexistence Patterns of Wild and Domestic Mammals in an Abandoned Landscape #########
################################### Zuleger, Perino & Pereira (2024) ###############################################
####################################################################################################################

############################################# Load required data ###################################################

set.seed(123)

library(tidyr)
library(psych)
library(unmarked)
library(AICcmodavg)
library(ggplot2)
library(reshape)
library(MASS)
library(gridExtra)

#### Load required data from GitHub

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

################################## Check for correlation between covariates #########################################

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

corPlot(site_covs,method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Site Covariates",xlas=2)
# removed CostDist_urban, verbatimElevation, Agriculture and Rock
corPlot(site_covs[, c(2,4,6,8:10,13)],method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Site Covariates",xlas=2)

corPlot(obs_covs[1:7],method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Observation Covariates",xlas=2)
# keep Temp_min, Rad_max and Rain_cum
corPlot(obs_covs[, c(1,4,6)],method="spearman",use = "complete.obs",scale=F, main="Correlation Plot Observation Covariates",xlas=2)

################################## Manual backwards elimination of variables ####################################### 

### Only shown for one species to save space, but this was done for every study species and all variables that remained significant (p < 0.1) for one of the species/parameters were passed on to the automated model selection. 

#### Domestic cattle #####

umfm_cattle <- unmarkedMultFrame(y = `Cattle`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_cattle) 

### Detection probability 

fm_cattle_p_global <- colext(~1, ~1, ~1, ~scale(Sensitivity)+scale(Temp_min)+scale(Rad_max)+scale(Rain_sum)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_cattle)
summary(fm_cattle_p_global)

fm_cattle_p_1 <- colext(~1, ~1, ~1, ~scale(Sensitivity)+scale(Temp_min)+scale(Rad_max)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_cattle)
summary(fm_cattle_p_1)

fm_cattle_p_2 <- colext(~1, ~1, ~1, ~scale(Sensitivity)+scale(Temp_min)+scale(Rad_max)+scale(HighShrub)+scale(LowShrub), umfm_cattle)
summary(fm_cattle_p_2)

fm_cattle_p_3 <- colext(~1, ~1, ~1, ~scale(Sensitivity)+scale(Rad_max)+scale(HighShrub)+scale(LowShrub), umfm_cattle)
summary(fm_cattle_p_3)

# Sensitivity, Radiation, High Shrub & Low shrub

### Occupancy

fm_cattle_psi_global <- colext(~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(Urban)+scale(LowShrub)+scale(HighShrub), ~1, ~1, ~1, umfm_cattle)
summary(fm_cattle_psi_global) 

fm_cattle_psi_1 <- colext(~scale(CostDist_road)+scale(OakForest)+scale(Urban)+scale(LowShrub)+scale(HighShrub), ~1, ~1, ~1, umfm_cattle)
summary(fm_cattle_psi_1) 

fm_cattle_psi_2 <- colext(~scale(OakForest)+scale(Urban)+scale(LowShrub)+scale(HighShrub), ~1, ~1, ~1, umfm_cattle)
summary(fm_cattle_psi_2) 

fm_cattle_psi_3 <- colext(~scale(OakForest)+scale(LowShrub)+scale(HighShrub), ~1, ~1, ~1, umfm_cattle)
summary(fm_cattle_psi_3) 

fm_cattle_psi_4 <- colext(~scale(OakForest)+scale(LowShrub), ~1, ~1, ~1, umfm_cattle)
summary(fm_cattle_psi_4)

fm_cattle_psi_5 <-  colext(~scale(LowShrub), ~1, ~1, ~1, umfm_cattle)
summary(fm_cattle_psi_5)

fm_cattle_psi_6 <- colext(~1, ~1, ~1, ~1, umfm_cattle)
summary(fm_cattle_psi_6)

### Colonization 

fm_cattle_col_global <-  colext(~1, ~scale(year)+scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(Urban)+scale(HighShrub)+scale(LowShrub), ~1, ~1, umfm_cattle)
summary(fm_cattle_col_global) 

fm_cattle_col_1 <- colext(~1, ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(Urban)+scale(HighShrub)+scale(LowShrub), ~1, ~1, umfm_cattle)
summary(fm_cattle_col_1) 

fm_cattle_col_2 <- colext(~1, ~scale(CostDist_water)+scale(OakForest)+scale(Urban)+scale(HighShrub)+scale(LowShrub), ~1, ~1, umfm_cattle)
summary(fm_cattle_col_2) 

fm_cattle_col_3 <- colext(~1, ~scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub), ~1, ~1, umfm_cattle)
summary(fm_cattle_col_3) 

fm_cattle_col_4 <- colext(~1, ~scale(CostDist_water)+scale(OakForest)+scale(LowShrub), ~1, ~1, umfm_cattle)
summary(fm_cattle_col_4) 

fm_cattle_col_5 <- colext(~1, ~scale(OakForest)+scale(LowShrub), ~1, ~1, umfm_cattle)
summary(fm_cattle_col_5) 

fm_cattle_col_6 <- colext(~1, ~scale(OakForest), ~1, ~1, umfm_cattle)
summary(fm_cattle_col_6)

# OakForest

### Extinction 

fm_cattle_ext_global <- colext(~1, ~1, ~scale(year)+scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(Urban)+scale(HighShrub)+scale(LowShrub), ~1, umfm_cattle) 
summary(fm_cattle_ext_global) 

fm_cattle_ext_1 <- colext(~1, ~1, ~scale(year)+scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub), ~1, umfm_cattle) 
summary(fm_cattle_ext_1) 

fm_cattle_ext_2 <- colext(~1, ~1, ~scale(year)+scale(CostDist_road)+scale(OakForest)+scale(HighShrub)+scale(LowShrub), ~1, umfm_cattle) 
summary(fm_cattle_ext_2) 

fm_cattle_ext_3 <- colext(~1, ~1, ~scale(CostDist_road)+scale(OakForest)+scale(HighShrub)+scale(LowShrub), ~1, umfm_cattle) 
summary(fm_cattle_ext_3) 

# CostDist_road, OakForest, HighShrub & LowShrub

################################## Automated model selection #####################################################

########### Create all possible model combinations  ############

psi_vars <- c("scale(HighShrub)")
gamma_vars <- c("scale(OakForest)")
epsilon_vars <- c("scale(CostDist_road)","scale(CostDist_water)","scale(OakForest)","scale(HighShrub)","scale(LowShrub)","scale(Urban)")
p_vars <- c("scale(Sensitivity)","scale(Rad_max)","scale(Temp_min)","scale(HighShrub)","scale(LowShrub)","scale(OakForest)")

psi_formulas <- list()
for (i in seq_along(psi_vars)) {
  tmp <- combn(psi_vars, i)
  tmp <- apply(tmp, 2, paste, collapse="+")
  tmp <- paste0("~", tmp)
  psi_formulas[[i]] <- tmp
}
psi_formulas <- unlist(psi_formulas)
psi_formulas[length(psi_formulas)+1] <- "~1"

gamma_formulas <- list()
for (i in seq_along(gamma_vars)) {
  tmp <- combn(gamma_vars, i)
  tmp <- apply(tmp, 2, paste, collapse="+")
  tmp <- paste0("~", tmp)
  gamma_formulas[[i]] <- tmp
}
gamma_formulas <- unlist(gamma_formulas)
gamma_formulas[length(gamma_formulas)+1] <- "~1"

epsilon_formulas <- list()
for (i in seq_along(epsilon_vars)) {
  tmp <- combn(epsilon_vars, i)
  tmp <- apply(tmp, 2, paste, collapse="+")
  tmp <- paste0("~", tmp)
  epsilon_formulas[[i]] <- tmp
}
epsilon_formulas <- unlist(epsilon_formulas)
epsilon_formulas[length(epsilon_formulas)+1] <- "~1"

p_formulas <- list()
for (i in seq_along(p_vars)) {
  tmp <- combn(p_vars, i)
  tmp <- apply(tmp, 2, paste, collapse="+")
  tmp <- paste0("~", tmp)
  p_formulas[[i]] <- tmp
}
p_formulas <- unlist(p_formulas)
p_formulas[length(p_formulas)+1] <- "~1"

models <- data.frame(psi_formula = 1:length(p_formulas), gamma_formula = 1:length(p_formulas),epsilon_formula = 1:length(p_formulas), p_formula = 1:length(p_formulas))

for(psi in 1:length(psi_formulas)) {
  psi_formula <- psi_formulas[psi]
  for(gamma in 1:length(gamma_formulas)) {
    gamma_formula <- gamma_formulas[gamma]
    for(epsilon in 1:length(epsilon_formulas)) {
      epsilon_formula <- epsilon_formulas[epsilon]
      for(p in 1:length(p_formulas)) {
        p_formula <- p_formulas[p]
        models[p,1] <- psi_formulas[psi]
        models[p,2] <- gamma_formulas[gamma]
        models[p,3] <- epsilon_formulas[epsilon]
        models[p,4] <- p_formulas[p]
      }
      if (epsilon==1) {models_epsilon <- models} else {models_epsilon <- rbind(models_epsilon, models)}
    }
    if (gamma==1) {models_gamma <- models_epsilon} else {models_gamma <- rbind(models_gamma,models_epsilon)}
  }
  if (psi==1) {models_all <- models_gamma} else {models_all <- rbind(models_all, models_gamma)}
}

########### Automated model selection #############

setwd(dir = file.path(t,"Dynamics_Coexistence_Peneda-main/Results"))

# Run models parallel on several cores because otherwise it will take forever #

library(parallel)
library(doParallel)

detectCores() # number of cores available
workers=makeCluster(50) # set number of cores used for analysis
registerDoParallel(workers)
getDoParWorkers()

###### Cattle ###### 

umfm_cattle <- unmarkedMultFrame(y = `Cattle`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_cattle) 

#### Create global model with variables obtained from backwards elimination ####

fm_cattle_global <- colext(~scale(HighShrub),
                           ~scale(OakForest),
                           ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban),
                           ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_cattle)

summary(fm_cattle_global)

#### Test goodness of fit of the global model to correct for overdispersion if necessary ####

gof_cattle_global <- mb.gof.test(fm_cattle_global, plot.hist=FALSE,nsim=1000)
gof_cattle_global
chat <- gof_cattle_global$c.hat.est
n <- 16

#### Run every possible model combination and estimate AIC / QAIC (based on c-hat value from goodness of fit test) ####

# Code is commented as without several cores it will take too long
# If several cores are not available skip to line 298

# models_cattle <- foreach(x = 1:nrow(models_all), .combine='rbind') %dopar% {
#    fm_cattle <- unmarked::colext(psiformula= as.formula(models_all$psi_formula[x]), gammaformula = as.formula(models_all$gamma_formula[x]), epsilonformula = as.formula(models_all$epsilon_formula[x]), pformula = as.formula(models_all$p_formula[x]), umfm_cattle)
#    k <- AICcmodavg::AICc(fm_cattle,return.K=TRUE)
#    if(chat > 1) {k <- k +1}
#    LogLike <- unmarked::logLik(fm_cattle)
#    c(psi_formula = models_all$psi_formula[x], gamma_formula = models_all$gamma_formula[x],epsilon_formula = models_all$epsilon_formula[x], p_formula = models_all$p_formula[x],AIC = fm_cattle@AIC, QAIC=(-2*LogLike/chat)+(2*k))
#  }
# 
# models_cattle <- as.data.frame(models_cattle,stringsAsFactors = FALSE) 
# names(models_cattle) <- c("psi_formula","gamma_formula","epsilon_formula","p_formula","AIC","QAIC")
# if(chat<1) {models_cattle$QAIC <- models_cattle$AIC}
# models_cattle$AIC <- as.numeric(models_cattle$AIC)
# models_cattle$QAIC <- as.numeric(models_cattle$QAIC)
#  
# write.table(models_cattle, "models_cattle.csv",sep=";",dec=".",row.names=F,col.names=T)
# 
# #### Create best model based on AIC / QAIC ####
# 
# best_model_QAIC_cattle <- models_cattle[models_cattle$QAIC == min(models_cattle$QAIC),]
#  
# fm_cattle_best_QAIC <- colext(psiformula = as.formula(best_model_QAIC_cattle$psi_formula), gammaformula = as.formula(best_model_QAIC_cattle$gamma_formula), epsilonformula = as.formula(best_model_QAIC_cattle$epsilon_formula), pformula = as.formula(best_model_QAIC_cattle$p_formula), umfm_cattle)
# summary(fm_cattle_best_QAIC)
#  
fm_cattle_best_QAIC <- colext(~1, ~scale(OakForest), ~scale(CostDist_road)+scale(OakForest)+scale(HighShrub)+scale(LowShrub),~scale(Sensitivity)+scale(Rad_max)+scale(HighShrub)+scale(LowShrub), umfm_cattle)
summary(fm_cattle_best_QAIC)

#### Backtransform model estimates from logit scale (keeping Intercept at 1) ####

variables_cattle <- data.frame(matrix(ncol = 10, nrow = 10))
colnames(variables_cattle) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_cattle$Species <- "Domestic cattle"
variables_cattle$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","CostDist_water","OakForest","HighShrub","LowShrub","Urban")

variables_cattle[1,3] <- backTransform(fm_cattle_best_QAIC,type="psi")@estimate
variables_cattle[1,4] <- SE(backTransform(fm_cattle_best_QAIC,type="psi"))

variables_cattle[1,5] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0),type="col"))@estimate
variables_cattle[1,6] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0),type="col")))
variables_cattle[7,5] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,1),type="col"))@estimate
variables_cattle[7,6] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,1),type="col")))

variables_cattle[1,7] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,0),type="ext"))@estimate
variables_cattle[1,8] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,0),type="ext")))
variables_cattle[5,7] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,1,0,0,0),type="ext"))@estimate
variables_cattle[5,8] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,1,0,0,0),type="ext")))
variables_cattle[7,7] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,1,0,0),type="ext"))@estimate
variables_cattle[7,8] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,1,0,0),type="ext")))
variables_cattle[8,7] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,1,0),type="ext"))@estimate
variables_cattle[8,8] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,1,0),type="ext")))
variables_cattle[9,7] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,1),type="ext"))@estimate
variables_cattle[9,8] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,1),type="ext")))

variables_cattle[1,9] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,0),type="det"))@estimate
variables_cattle[1,10] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,0),type="det")))
variables_cattle[2,9] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,1,0,0,0),type="det"))@estimate
variables_cattle[2,10] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,1,0,0,0),type="det")))
variables_cattle[3,9] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,1,0,0),type="det"))@estimate
variables_cattle[3,10] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,1,0,0),type="det")))
variables_cattle[8,9] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,1,0),type="det"))@estimate
variables_cattle[8,10] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,1,0),type="det")))
variables_cattle[9,9] <- backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,1),type="det"))@estimate
variables_cattle[9,10] <- SE(backTransform(linearComb(fm_cattle_best_QAIC,c(1,0,0,0,1),type="det")))

#### Get occupancy estimates for the whole study area with non-parametric bootstrapping ####

m3 <- nonparboot(fm_cattle_best_QAIC,B = 1000)
occ_pred_cattle_best_QAIC <- data.frame(year = c(2015:2022),
                                        smoothed_occ = smoothed(fm_cattle_best_QAIC)[2,],
                                        SE = m3@smoothed.mean.bsse[2,])
occ_pred_cattle_best_QAIC$Species <- 'Domestic cattle'

cattle_year <- lm(smoothed_occ ~ year, occ_pred_cattle_best_QAIC)
summary(cattle_year)

#### Get site-specific occupancy estimates obtained from model output ####

occ_cattle_site <- data.frame(t(fm_cattle_best_QAIC@smoothed[2, , ]))
colnames(occ_cattle_site) <- unique(yearlySiteCovs(umfm_cattle)$year)
occ_cattle_site$grid <- unique(Cattle$grid) 

#write.table(occ_cattle_site, "Occ_cattle_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_cattle_site_long <- data.frame(melt(occ_cattle_site, id.vars="grid"))
names(occ_cattle_site_long) <- c("Grid","Year","Occupancy")
occ_cattle_site_long$Year <- as.numeric(as.character(occ_cattle_site_long$Year))
occ_cattle_site_long$Species <- "Domestic cattle"

###### Horse ###### 

umfm_horse <- unmarkedMultFrame(y = `Horse`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_horse) 

#### Create global model with variables obtained from backwards elimination ####

fm_horse_global <- colext(~scale(HighShrub),
                           ~scale(OakForest),
                           ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban),
                           ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_horse)

summary(fm_horse_global)

#### Test goodness of fit of the global model to correct for overdispersion if necessary ####

gof_horse_global <- mb.gof.test(fm_horse_global, plot.hist=FALSE,nsim=1000)
gof_horse_global
chat <- gof_horse_global$c.hat.est

#### Run every possible model combination and estimate AIC / QAIC (based on c-hat value from goodness of fit test) ####

# models_horse <- foreach(x = 1:nrow(models_all), .combine='rbind') %dopar% {
#    fm_horse <- unmarked::colext(psiformula= as.formula(models_all$psi_formula[x]), gammaformula = as.formula(models_all$gamma_formula[x]), epsilonformula = as.formula(models_all$epsilon_formula[x]), pformula = as.formula(models_all$p_formula[x]), umfm_horse)
#    k <- AICcmodavg::AICc(fm_horse,return.K=TRUE)
#    if(chat > 1) {k <- k +1}
#    LogLike <- unmarked::logLik(fm_horse)
#    c(psi_formula = models_all$psi_formula[x], gamma_formula = models_all$gamma_formula[x],epsilon_formula = models_all$epsilon_formula[x], p_formula = models_all$p_formula[x],AIC = fm_horse@AIC, QAIC=(-2*LogLike/chat)+(2*k))
# }
# 
# models_horse <- as.data.frame(models_horse,stringsAsFactors = FALSE)
# names(models_horse) <- c("psi_formula","gamma_formula","epsilon_formula","p_formula","AIC","QAIC")
# if(chat<1) {models_horse$QAIC <- models_horse$AIC}
# models_horse$AIC <- as.numeric(models_horse$AIC)
# models_horse$QAIC <- as.numeric(models_horse$QAIC)
# 
# #write.table(models_horse, "models_horse.csv",sep=";",dec=".",row.names=F,col.names=T)
# 
# #### Create best model based on AIC / QAIC ####
# 
# best_model_QAIC_horse <- models_horse[models_horse$QAIC == min(models_horse$QAIC),]
# 
# fm_horse_best_QAIC <- colext(psiformula = as.formula(best_model_QAIC_horse$psi_formula), gammaformula = as.formula(best_model_QAIC_horse$gamma_formula), epsilonformula = as.formula(best_model_QAIC_horse$epsilon_formula), pformula = as.formula(best_model_QAIC_horse$p_formula), umfm_horse)
# summary(fm_horse_best_QAIC)

fm_horse_best_QAIC <- colext(~1, ~scale(OakForest), ~scale(LowShrub)+scale(Urban), ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(OakForest)+scale(HighShrub)+scale(LowShrub), umfm_horse)
summary(fm_horse_best_AIC)

#### Backtransform model estimates from logit scale (keeping Intercept at 1) ####

variables_horse <- data.frame(matrix(ncol = 10, nrow = 10))
colnames(variables_horse) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_horse$Species <- "Domestic horse"
variables_horse$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","CostDist_water","OakForest","HighShrub","LowShrub","Urban")

variables_horse[1,3] <- backTransform(fm_horse_best_QAIC,type="psi")@estimate
variables_horse[1,4] <- SE(backTransform(fm_horse_best_QAIC,type="psi"))

variables_horse[1,5] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0),type="col"))@estimate
variables_horse[1,6] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0),type="col")))
variables_horse[7,5] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,1),type="col"))@estimate
variables_horse[7,6] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,1),type="col")))

variables_horse[1,7] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0),type="ext"))@estimate
variables_horse[1,8] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0),type="ext")))
variables_horse[9,7] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,1,0),type="ext"))@estimate
variables_horse[9,8] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,1,0),type="ext")))
variables_horse[10,7] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,1),type="ext"))@estimate
variables_horse[10,8] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,1),type="ext")))

variables_horse[1,9] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,0,0,0),type="det"))@estimate
variables_horse[1,10] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,0,0,0),type="det")))
variables_horse[2,9] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,1,0,0,0,0,0),type="det"))@estimate
variables_horse[2,10] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,1,0,0,0,0,0),type="det")))
variables_horse[3,9] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,1,0,0,0,0),type="det"))@estimate
variables_horse[3,10] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,1,0,0,0,0),type="det")))
variables_horse[4,9] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,1,0,0,0),type="det"))@estimate
variables_horse[4,10] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,1,0,0,0),type="det")))
variables_horse[7,9] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,1,0,0),type="det"))@estimate
variables_horse[7,10] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,1,0,0),type="det")))
variables_horse[8,9] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,0,1,0),type="det"))@estimate
variables_horse[8,10] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,0,1,0),type="det")))
variables_horse[9,9] <- backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,0,0,1),type="det"))@estimate
variables_horse[9,10] <- SE(backTransform(linearComb(fm_horse_best_QAIC,c(1,0,0,0,0,0,1),type="det")))

#### Get occupancy estimates for the whole study area with non-parametric bootstrapping ####

m3 <- nonparboot(fm_horse_best_QAIC,B = 1000)
occ_pred_horse_best_QAIC <- data.frame(year = c(2015:2022),
                                        smoothed_occ = smoothed(fm_horse_best_QAIC)[2,],
                                        SE = m3@smoothed.mean.bsse[2,])
occ_pred_horse_best_QAIC$Species <- 'Domestic horse'

horse_year <- lm(smoothed_occ ~ year, occ_pred_horse_best_QAIC)
summary(horse_year)

#### Get site-specific occupancy estimates obtained from model output ####

occ_horse_site <- data.frame(t(fm_horse_best_QAIC@smoothed[2, , ]))
colnames(occ_horse_site) <- unique(yearlySiteCovs(umfm_horse)$year)
occ_horse_site$grid <- unique(Horse$grid) 

#write.table(occ_horse_site, "Occ_horse_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_horse_site_long <- data.frame(melt(occ_horse_site, id.vars="grid"))
names(occ_horse_site_long) <- c("Grid","Year","Occupancy")
occ_horse_site_long$Year <- as.numeric(as.character(occ_horse_site_long$Year))
occ_horse_site_long$Species <- "Domestic horse"

###### Roe deer ###### 

umfm_deer <- unmarkedMultFrame(y = `Deer`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_deer) 

#### Create global model with variables obtained from backwards elimination ####

fm_deer_global <- colext(~scale(HighShrub),
                          ~scale(OakForest),
                          ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban),
                          ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_deer)

#### Test goodness of fit of the global model to correct for overdispersion if necessary ####

gof_deer_global <- mb.gof.test(fm_deer_global, plot.hist=FALSE,nsim=1000)
gof_deer_global
chat <- gof_deer_global$c.hat.est

#### Run every possible model combination and estimate AIC / QAIC (based on c-hat value from goodness of fit test) ####

# models_deer <- foreach(x = 1:nrow(models_all), .combine='rbind') %dopar% {
#  fm_deer <- unmarked::colext(psiformula= as.formula(models_all$psi_formula[x]), gammaformula = as.formula(models_all$gamma_formula[x]), epsilonformula = as.formula(models_all$epsilon_formula[x]), pformula = as.formula(models_all$p_formula[x]), umfm_deer)
#  k <- AICcmodavg::AICc(fm_deer,return.K=TRUE)
#  if(chat > 1) {k <- k +1}
#  LogLike <- unmarked::logLik(fm_deer)
#  c(psi_formula = models_all$psi_formula[x], gamma_formula = models_all$gamma_formula[x],epsilon_formula = models_all$epsilon_formula[x], p_formula = models_all$p_formula[x],AIC = fm_deer@AIC, QAIC=(-2*LogLike/chat)+(2*k))
# }
# 
# models_deer <- as.data.frame(models_deer,stringsAsFactors = FALSE) 
# names(models_deer) <- c("psi_formula","gamma_formula","epsilon_formula","p_formula","AIC","QAIC")
# if(chat<1) {models_deer$QAIC <- models_deer$AIC}
# models_deer$AIC <- as.numeric(models_deer$AIC)
# models_deer$QAIC <- as.numeric(models_deer$QAIC)
#  
# #write.table(models_deer, "models_deer.csv",sep=";",dec=".",row.names=F,col.names=T)
# 
# #### Create best model based on AIC / QAIC ####
# 
# best_model_QAIC_deer <- models_deer[models_deer$QAIC == min(models_deer$QAIC),]
#  
# fm_deer_best_QAIC <- colext(psiformula = as.formula(best_model_QAIC_deer$psi_formula), gammaformula = as.formula(best_model_QAIC_deer$gamma_formula), epsilonformula = as.formula(best_model_QAIC_deer$epsilon_formula), pformula = as.formula(best_model_QAIC_deer$p_formula), umfm_deer)
# summary(fm_deer_best_QAIC)

fm_deer_best_QAIC <- colext(~1, ~1, ~scale(HighShrub), ~scale(Sensitivity)+scale(Rad_max)+scale(OakForest)+scale(LowShrub), umfm_deer)
summary(fm_deer_best_QAIC)

#### Backtransform model estimates from logit scale (keeping Intercept at 1) ####

variables_deer <- data.frame(matrix(ncol = 10, nrow = 10))
colnames(variables_deer) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_deer$Species <- "European roe deer"
variables_deer$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","CostDist_water","OakForest","HighShrub","LowShrub","Urban")

variables_deer[1,3] <- backTransform(fm_deer_best_QAIC,type="psi")@estimate
variables_deer[1,4] <- SE(backTransform(fm_deer_best_QAIC,type="psi"))

variables_deer[1,5] <- backTransform(fm_deer_best_QAIC,type="col")@estimate
variables_deer[1,6] <- SE(backTransform(fm_deer_best_QAIC,type="col"))

variables_deer[1,7] <- backTransform(linearComb(fm_deer_best_QAIC,c(1,0),type="ext"))@estimate
variables_deer[1,8] <- SE(backTransform(linearComb(fm_deer_best_QAIC,c(1,0),type="ext")))
variables_deer[8,7] <- backTransform(linearComb(fm_deer_best_QAIC,c(1,1),type="ext"))@estimate
variables_deer[8,8] <- SE(backTransform(linearComb(fm_deer_best_QAIC,c(1,1),type="ext")))

variables_deer[1,9] <- backTransform(linearComb(fm_deer_best_QAIC,c(1,0,0,0,0),type="det"))@estimate
variables_deer[1,10] <- SE(backTransform(linearComb(fm_deer_best_QAIC,c(1,0,0,0,0),type="det")))
variables_deer[2,9] <- backTransform(linearComb(fm_deer_best_QAIC,c(1,1,0,0,0),type="det"))@estimate
variables_deer[2,10] <- SE(backTransform(linearComb(fm_deer_best_QAIC,c(1,1,0,0,0),type="det")))
variables_deer[3,9] <- backTransform(linearComb(fm_deer_best_QAIC,c(1,0,1,0,0),type="det"))@estimate
variables_deer[3,10] <- SE(backTransform(linearComb(fm_deer_best_QAIC,c(1,0,1,0,0),type="det")))
variables_deer[7,9] <- backTransform(linearComb(fm_deer_best_QAIC,c(1,0,0,1,0),type="det"))@estimate
variables_deer[7,10] <- SE(backTransform(linearComb(fm_deer_best_QAIC,c(1,0,0,1,0),type="det")))
variables_deer[9,9] <- backTransform(linearComb(fm_deer_best_QAIC,c(1,0,0,0,1),type="det"))@estimate
variables_deer[9,10] <- SE(backTransform(linearComb(fm_deer_best_QAIC,c(1,0,0,0,1),type="det")))

#### Get occupancy estimates for the whole study area with non-parametric bootstrapping ####

m3 <- nonparboot(fm_deer_best_QAIC,B = 1000)
occ_pred_deer_best_QAIC <- data.frame(year = c(2015:2022),
                                       smoothed_occ = smoothed(fm_deer_best_QAIC)[2,],
                                       SE = m3@smoothed.mean.bsse[2,])
occ_pred_deer_best_QAIC$Species <- 'European roe deer'

deer_year <- lm(smoothed_occ ~ year, occ_pred_deer_best_QAIC)
summary(deer_year)

#### Get site-specific occupancy estimates obtained from model output ####

occ_deer_site <- data.frame(t(fm_deer_best_QAIC@smoothed[2, , ]))
colnames(occ_deer_site) <- unique(yearlySiteCovs(umfm_deer)$year)
occ_deer_site$grid <- unique(Deer$grid) 

#write.table(occ_deer_site, "Occ_deer_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_deer_site_long <- data.frame(melt(occ_deer_site, id.vars="grid"))
names(occ_deer_site_long) <- c("Grid","Year","Occupancy")
occ_deer_site_long$Year <- as.numeric(as.character(occ_deer_site_long$Year))
occ_deer_site_long$Species <- "European roe deer"

###### Wild boar ###### 

umfm_boar <- unmarkedMultFrame(y = `Boar`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_boar) 

#### Create global model with variables obtained from backwards elimination ####

fm_boar_global <- colext(~scale(HighShrub),
                         ~scale(OakForest),
                         ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban),
                         ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_boar)

#### Test goodness of fit of the global model to correct for overdispersion if necessary ####

gof_boar_global <- mb.gof.test(fm_boar_global, plot.hist=FALSE,nsim=1000)
gof_boar_global
chat <- gof_boar_global$c.hat.est

#### Run every possible model combination and estimate AIC / QAIC (based on c-hat value from goodness of fit test) ####

# models_boar <- foreach(x = 1:nrow(models_all), .combine='rbind') %dopar% {
#  fm_boar <- unmarked::colext(psiformula= as.formula(models_all$psi_formula[x]), gammaformula = as.formula(models_all$gamma_formula[x]), epsilonformula = as.formula(models_all$epsilon_formula[x]), pformula = as.formula(models_all$p_formula[x]), umfm_boar)
#  k <- AICcmodavg::AICc(fm_boar,return.K=TRUE)
#  if(chat > 1) {k <- k +1}
#  LogLike <- unmarked::logLik(fm_boar)
#  c(psi_formula = models_all$psi_formula[x], gamma_formula = models_all$gamma_formula[x],epsilon_formula = models_all$epsilon_formula[x], p_formula = models_all$p_formula[x],AIC = fm_boar@AIC, QAIC=(-2*LogLike/chat)+(2*k))
# }
# 
# models_boar <- as.data.frame(models_boar,stringsAsFactors = FALSE) 
# names(models_boar) <- c("psi_formula","gamma_formula","epsilon_formula","p_formula","AIC","QAIC")
# if(chat<1) {models_boar$QAIC <- models_boar$AIC}
# models_boar$AIC <- as.numeric(models_boar$AIC)
# models_boar$QAIC <- as.numeric(models_boar$QAIC)
#  
# #write.table(models_boar, "models_boar.csv",sep=";",dec=".",row.names=F,col.names=T)
# 
# #### Create best model based on AIC / QAIC ####
# 
# best_model_QAIC_boar <- models_boar[models_boar$QAIC == min(models_boar$QAIC),]
#  
# fm_boar_best_QAIC <- colext(psiformula = as.formula(best_model_QAIC_boar$psi_formula), gammaformula = as.formula(best_model_QAIC_boar$gamma_formula), epsilonformula = as.formula(best_model_QAIC_boar$epsilon_formula), pformula = as.formula(best_model_QAIC_boar$p_formula), umfm_boar)
# summary(fm_boar_best_QAIC)

fm_boar_best_QAIC <- colext(~scale(HighShrub), ~1, ~1, ~scale(Sensitivity)+scale(OakForest), umfm_boar)
summary(fm_boar_best_QAIC)

#### Backtransform model estimates from logit scale (keeping Intercept at 1) ####

variables_boar <- data.frame(matrix(ncol = 10, nrow = 10))
colnames(variables_boar) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_boar$Species <- "Wild boar"
variables_boar$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","CostDist_water","OakForest","HighShrub","LowShrub","Urban")

variables_boar[1,3] <- backTransform(linearComb(fm_boar_best_QAIC,c(1,0),type="psi"))@estimate
variables_boar[1,4] <- SE(backTransform(linearComb(fm_boar_best_QAIC,c(1,0),type="psi")))
variables_boar[8,3] <- backTransform(linearComb(fm_boar_best_QAIC,c(1,1),type="psi"))@estimate
variables_boar[8,4] <- SE(backTransform(linearComb(fm_boar_best_QAIC,c(1,1),type="psi")))

variables_boar[1,5] <- backTransform(fm_boar_best_QAIC,type="col")@estimate
variables_boar[1,6] <- SE(backTransform(fm_boar_best_QAIC,type="col"))

variables_boar[1,7] <- backTransform(fm_boar_best_QAIC,type="ext")@estimate
variables_boar[1,8] <- SE(backTransform(fm_boar_best_QAIC,type="ext"))

variables_boar[1,9] <- backTransform(linearComb(fm_boar_best_QAIC,c(1,0,0),type="det"))@estimate
variables_boar[1,10] <- SE(backTransform(linearComb(fm_boar_best_QAIC,c(1,0,0),type="det")))
variables_boar[2,9] <- backTransform(linearComb(fm_boar_best_QAIC,c(1,1,0),type="det"))@estimate
variables_boar[2,10] <- SE(backTransform(linearComb(fm_boar_best_QAIC,c(1,1,0),type="det")))
variables_boar[7,9] <- backTransform(linearComb(fm_boar_best_QAIC,c(1,0,1),type="det"))@estimate
variables_boar[7,10] <- SE(backTransform(linearComb(fm_boar_best_QAIC,c(1,0,1),type="det")))


#### Get occupancy estimates for the whole study area with non-parametric bootstrapping ####

m3 <- nonparboot(fm_boar_best_QAIC,B = 1000)
occ_pred_boar_best_QAIC <- data.frame(year = c(2015:2022),
                                      smoothed_occ = smoothed(fm_boar_best_QAIC)[2,],
                                      SE = m3@smoothed.mean.bsse[2,])
occ_pred_boar_best_QAIC$Species <- 'Wild boar'

boar_year <- lm(smoothed_occ ~ year, occ_pred_boar_best_QAIC)
summary(boar_year)

#### Get site-specific occupancy estimates obtained from model output ####

occ_boar_site <- data.frame(t(fm_boar_best_QAIC@smoothed[2, , ]))
colnames(occ_boar_site) <- unique(yearlySiteCovs(umfm_boar)$year)
occ_boar_site$grid <- unique(Boar$grid) 

#write.table(occ_boar_site, "Occ_boar_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_boar_site_long <- data.frame(melt(occ_boar_site, id.vars="grid"))
names(occ_boar_site_long) <- c("Grid","Year","Occupancy")
occ_boar_site_long$Year <- as.numeric(as.character(occ_boar_site_long$Year))
occ_boar_site_long$Species <- "Wild boar"

###### Red fox ###### 

umfm_fox <- unmarkedMultFrame(y = `Fox`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_fox) 

#### Create global model with variables obtained from backwards elimination ####

fm_fox_global <- colext(~scale(HighShrub),
                         ~scale(OakForest),
                         ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban),
                         ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_fox)

#### Test goodness of fit of the global model to correct for overdispersion if necessary ####

gof_fox_global <- mb.gof.test(fm_fox_global, plot.hist=FALSE,nsim=1000)
gof_fox_global
chat <- gof_fox_global$c.hat.est

#### Run every possible model combination and estimate AIC / QAIC (based on c-hat value from goodness of fit test) ####

# models_fox <- foreach(x = 1:nrow(models_all), .combine='rbind') %dopar% {
#  fm_fox <- unmarked::colext(psiformula= as.formula(models_all$psi_formula[x]), gammaformula = as.formula(models_all$gamma_formula[x]), epsilonformula = as.formula(models_all$epsilon_formula[x]), pformula = as.formula(models_all$p_formula[x]), umfm_fox)
#  k <- AICcmodavg::AICc(fm_fox,return.K=TRUE)
#  if(chat > 1) {k <- k +1}
#  LogLike <- unmarked::logLik(fm_fox)
#  c(psi_formula = models_all$psi_formula[x], gamma_formula = models_all$gamma_formula[x],epsilon_formula = models_all$epsilon_formula[x], p_formula = models_all$p_formula[x],AIC = fm_fox@AIC, QAIC=(-2*LogLike/chat)+(2*k))
# }
# 
# models_fox <- as.data.frame(models_fox,stringsAsFactors = FALSE) 
# names(models_fox) <- c("psi_formula","gamma_formula","epsilon_formula","p_formula","AIC","QAIC")
# if(chat<1) {models_fox$QAIC <- models_fox$AIC}
# models_fox$AIC <- as.numeric(models_fox$AIC)
# models_fox$QAIC <- as.numeric(models_fox$QAIC)
#  
# #write.table(models_fox, "models_fox.csv",sep=";",dec=".",row.names=F,col.names=T)
# 
# #### Create best model based on AIC / QAIC ####
# 
# best_model_QAIC_fox <- models_fox[models_fox$QAIC == min(models_fox$QAIC),]
#  
# fm_fox_best_QAIC <- colext(psiformula = as.formula(best_model_QAIC_fox$psi_formula), gammaformula = as.formula(best_model_QAIC_fox$gamma_formula), epsilonformula = as.formula(best_model_QAIC_fox$epsilon_formula), pformula = as.formula(best_model_QAIC_fox$p_formula), umfm_fox)
# summary(fm_fox_best_QAIC)

fm_fox_best_QAIC <- colext(~1, ~1, ~1, ~scale(Sensitivity)+scale(OakForest)+scale(LowShrub), umfm_fox)
summary(fm_fox_best_QAIC)
# NaNs produced --> potential convergence issues --> simplify model
# OakForest has lowest effect -> remove OakForest
k <- AICcmodavg::AICc(fm_fox_best_QAIC,return.K=TRUE)
(-2*logLik(fm_fox_best_QAIC)/chat)+(2*k)

fm_fox_best_QAIC <- colext(~1, ~1, ~1, ~scale(Sensitivity)+scale(LowShrub), umfm_fox)
summary(fm_fox_best_QAIC)
# no more convergence issues
k <- AICcmodavg::AICc(fm_fox_best_QAIC,return.K=TRUE)
(-2*logLik(fm_fox_best_QAIC)/chat)+(2*k)
# QAIC is very close to QAIC from original model


#### Backtransform model estimates from logit scale (keeping Intercept at 1) ####

variables_fox <- data.frame(matrix(ncol = 10, nrow = 10))
colnames(variables_fox) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_fox$Species <- "Red fox"
variables_fox$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","CostDist_water","OakForest","HighShrub","LowShrub","Urban")

variables_fox[1,3] <- backTransform(fm_fox_best_QAIC,type="psi")@estimate
variables_fox[1,4] <- SE(backTransform(fm_fox_best_QAIC,type="psi"))

variables_fox[1,5] <- backTransform(fm_fox_best_QAIC,type="col")@estimate
variables_fox[1,6] <- SE(backTransform(fm_fox_best_QAIC,type="col"))

variables_fox[1,7] <- backTransform(fm_fox_best_QAIC,type="ext")@estimate
variables_fox[1,8] <- SE(backTransform(fm_fox_best_QAIC,type="ext"))

variables_fox[1,9] <- backTransform(linearComb(fm_fox_best_QAIC,c(1,0,0),type="det"))@estimate
variables_fox[1,10] <- SE(backTransform(linearComb(fm_fox_best_QAIC,c(1,0,0),type="det")))
variables_fox[2,9] <- backTransform(linearComb(fm_fox_best_QAIC,c(1,1,0),type="det"))@estimate
variables_fox[2,10] <- SE(backTransform(linearComb(fm_fox_best_QAIC,c(1,1,0),type="det")))
variables_fox[9,9] <- backTransform(linearComb(fm_fox_best_QAIC,c(1,0,1),type="det"))@estimate
variables_fox[9,10] <- SE(backTransform(linearComb(fm_fox_best_QAIC,c(1,0,1),type="det")))


#### Get occupancy estimates for the whole study area with non-parametric bootstrapping ####

m3 <- nonparboot(fm_fox_best_QAIC,B = 1000)
occ_pred_fox_best_QAIC <- data.frame(year = c(2015:2022),
                                      smoothed_occ = smoothed(fm_fox_best_QAIC)[2,],
                                      SE = m3@smoothed.mean.bsse[2,])
occ_pred_fox_best_QAIC$Species <- 'Red fox'

fox_year <- lm(smoothed_occ ~ year, occ_pred_fox_best_QAIC)
summary(fox_year)

#### Get site-specific occupancy estimates obtained from model output ####

occ_fox_site <- data.frame(t(fm_fox_best_QAIC@smoothed[2, , ]))
colnames(occ_fox_site) <- unique(yearlySiteCovs(umfm_fox)$year)
occ_fox_site$grid <- unique(Fox$grid) 

#write.table(occ_fox_site, "Occ_fox_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_fox_site_long <- data.frame(melt(occ_fox_site, id.vars="grid"))
names(occ_fox_site_long) <- c("Grid","Year","Occupancy")
occ_fox_site_long$Year <- as.numeric(as.character(occ_fox_site_long$Year))
occ_fox_site_long$Species <- "Red fox"

###### Gray wolf ###### 

umfm_wolf <- unmarkedMultFrame(y = `Wolf`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_wolf) 

#### Create global model with variables obtained from backwards elimination ####

fm_wolf_global <- colext(~scale(HighShrub),
                        ~scale(OakForest),
                        ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban),
                        ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_wolf)

#### Test goodness of fit of the global model to correct for overdispersion if necessary ####

gof_wolf_global <- mb.gof.test(fm_wolf_global, plot.hist=FALSE,nsim=1000)
gof_wolf_global
chat <- gof_wolf_global$c.hat.est

#### Run every possible model combination and estimate AIC / QAIC (based on c-hat value from goodness of fit test) ####

# models_wolf <- foreach(x = 1:nrow(models_all), .combine='rbind') %dopar% {
#  fm_wolf <- unmarked::colext(psiformula= as.formula(models_all$psi_formula[x]), gammaformula = as.formula(models_all$gamma_formula[x]), epsilonformula = as.formula(models_all$epsilon_formula[x]), pformula = as.formula(models_all$p_formula[x]), umfm_wolf)
#  k <- AICcmodavg::AICc(fm_wolf,return.K=TRUE)
#  if(chat > 1) {k <- k +1}
#  LogLike <- unmarked::logLik(fm_wolf)
#  c(psi_formula = models_all$psi_formula[x], gamma_formula = models_all$gamma_formula[x],epsilon_formula = models_all$epsilon_formula[x], p_formula = models_all$p_formula[x],AIC = fm_wolf@AIC, QAIC=(-2*LogLike/chat)+(2*k))
# }
# 
# models_wolf <- as.data.frame(models_wolf,stringsAsFactors = FALSE) 
# names(models_wolf) <- c("psi_formula","gamma_formula","epsilon_formula","p_formula","AIC","QAIC")
# if(chat<1) {models_wolf$QAIC <- models_wolf$AIC}
# models_wolf$AIC <- as.numeric(models_wolf$AIC)
# models_wolf$QAIC <- as.numeric(models_wolf$QAIC)
#  
# #write.table(models_wolf, "models_wolf.csv",sep=";",dec=".",row.names=F,col.names=T)
# 
# #### Create best model based on AIC / QAIC ####
# 
# best_model_QAIC_wolf <- models_wolf[models_wolf$AIC == min(models_wolf$QAIC),]
#  
# fm_wolf_best_QAIC <- colext(psiformula = as.formula(best_model_QAIC_wolf$psi_formula), gammaformula = as.formula(best_model_QAIC_wolf$gamma_formula), epsilonformula = as.formula(best_model_QAIC_wolf$epsilon_formula), pformula = as.formula(best_model_QAIC_wolf$p_formula), umfm_wolf)
# summary(fm_wolf_best_QAIC)

fm_wolf_best_QAIC <- colext(~1, ~1, ~scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban), ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub), umfm_wolf)
summary(fm_wolf_best_QAIC)
# Convergence issues --> unreasonable estimates for extinction (especially urban)
k <- AICcmodavg::AICc(fm_wolf_best_QAIC,return.K=TRUE)
(-2*logLik(fm_wolf_best_QAIC)/chat)+(2*k)

fm_wolf_best_QAIC <- colext(~1, ~1, ~scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub), ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub), umfm_wolf)
summary(fm_wolf_best_QAIC)
# no convergence issues
k <- AICcmodavg::AICc(fm_wolf_best_QAIC,return.K=TRUE)
(-2*logLik(fm_wolf_best_QAIC)/chat)+(2*k)

#### Backtransform model estimates from logit scale (keeping Intercept at 1) ####

variables_wolf <- data.frame(matrix(ncol = 10, nrow = 10))
colnames(variables_wolf) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_wolf$Species <- "Gray wolf"
variables_wolf$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","CostDist_water","OakForest","HighShrub","LowShrub","Urban")

variables_wolf[1,3] <- backTransform(fm_wolf_best_QAIC,type="psi")@estimate
variables_wolf[1,4] <- SE(backTransform(fm_wolf_best_QAIC,type="psi"))

variables_wolf[1,5] <- backTransform(fm_wolf_best_QAIC,type="col")@estimate
variables_wolf[1,6] <- SE(backTransform(fm_wolf_best_QAIC,type="col"))

variables_wolf[1,7] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,0),type="ext"))@estimate
variables_wolf[1,8] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,0),type="ext")))
variables_wolf[6,7] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,1,0,0,0),type="ext"))@estimate
variables_wolf[6,8] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,1,0,0,0),type="ext")))
variables_wolf[7,7] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,1,0,0),type="ext"))@estimate
variables_wolf[7,8] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,1,0,0),type="ext")))
variables_wolf[8,7] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,1,0),type="ext"))@estimate
variables_wolf[8,8] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,1,0),type="ext")))
variables_wolf[9,7] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,1),type="ext"))@estimate
variables_wolf[9,8] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,1),type="ext")))

variables_wolf[1,9] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,0,0),type="det"))@estimate
variables_wolf[1,10] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,0,0),type="det")))
variables_wolf[2,9] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,1,0,0,0,0),type="det"))@estimate
variables_wolf[2,10] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,1,0,0,0,0),type="det")))
variables_wolf[3,9] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,1,0,0,0),type="det"))@estimate
variables_wolf[3,10] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,1,0,0,0),type="det")))
variables_wolf[4,9] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,1,0,0),type="det"))@estimate
variables_wolf[4,10] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,1,0,0),type="det")))
variables_wolf[8,9] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,1,0),type="det"))@estimate
variables_wolf[8,10] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,1,0),type="det")))
variables_wolf[9,9] <- backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,0,1),type="det"))@estimate
variables_wolf[9,10] <- SE(backTransform(linearComb(fm_wolf_best_QAIC,c(1,0,0,0,0,1),type="det")))

m3 <- nonparboot(fm_wolf_best_QAIC,B = 1000)
occ_pred_wolf_best_QAIC <- data.frame(year = c(2015:2022),
                                     smoothed_occ = smoothed(fm_wolf_best_QAIC)[2,],
                                     SE = m3@smoothed.mean.bsse[2,])
occ_pred_wolf_best_QAIC$Species <- 'Gray wolf'

wolf_year <- lm(smoothed_occ ~ year, occ_pred_wolf_best_QAIC)
summary(wolf_year)

#### Get site-specific occupancy estimates obtained from model output ####

occ_wolf_site <- data.frame(t(fm_wolf_best_QAIC@smoothed[2, , ]))
colnames(occ_wolf_site) <- unique(yearlySiteCovs(umfm_wolf)$year)
occ_wolf_site$grid <- unique(Wolf$grid) 

#write.table(occ_wolf_site, "Occ_wolf_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_wolf_site_long <- data.frame(melt(occ_wolf_site, id.vars="grid"))
names(occ_wolf_site_long) <- c("Grid","Year","Occupancy")
occ_wolf_site_long$Year <- as.numeric(as.character(occ_wolf_site_long$Year))
occ_wolf_site_long$Species <- "Gray wolf"

###### Iberian ibex ###### 

#### Create global model with variables obtained from backwards elimination ####

umfm_ibex <- unmarkedMultFrame(y = `Ibex`[-c(1,2)],siteCovs = SiteCovs_2015, obsCovs = ObsCovs[c(1,4,6)], numPrimary = 8,yearlySiteCovs = yearlySiteCovs)
summary(umfm_ibex) 

fm_ibex_global <- colext(~scale(HighShrub),
                         ~scale(OakForest),
                         ~scale(CostDist_road)+scale(CostDist_water)+scale(OakForest)+scale(HighShrub)+scale(LowShrub)+scale(Urban),
                         ~scale(Sensitivity)+scale(Rad_max)+scale(Temp_min)+scale(HighShrub)+scale(LowShrub)+scale(OakForest), umfm_ibex)

#### Test goodness of fit of the global model to correct for overdispersion if necessary ####

gof_ibex_global <- mb.gof.test(fm_ibex_global, plot.hist=FALSE,nsim=1000)
gof_ibex_global
chat <- gof_ibex_global$c.hat.est

#### Run every possible model combination and estimate AIC / QAIC (based on c-hat value from goodness of fit test) ####

# models_ibex <- foreach(x = 1:nrow(models_all), .combine='rbind') %dopar% {
#  fm_ibex <- unmarked::colext(psiformula= as.formula(models_all$psi_formula[x]), gammaformula = as.formula(models_all$gamma_formula[x]), epsilonformula = as.formula(models_all$epsilon_formula[x]), pformula = as.formula(models_all$p_formula[x]), umfm_ibex)
#  k <- AICcmodavg::AICc(fm_ibex,return.K=TRUE)
#  if(chat > 1) {k <- k +1}
#  LogLike <- unmarked::logLik(fm_ibex)
#  c(psi_formula = models_all$psi_formula[x], gamma_formula = models_all$gamma_formula[x],epsilon_formula = models_all$epsilon_formula[x], p_formula = models_all$p_formula[x],AIC = fm_ibex@AIC, QAIC=(-2*LogLike/chat)+(2*k))
# }
# 
# models_ibex <- as.data.frame(models_ibex,stringsAsFactors = FALSE) 
# names(models_ibex) <- c("psi_formula","gamma_formula","epsilon_formula","p_formula","AIC","QAIC")
# if(chat<1) {models_ibex$QAIC <- models_ibex$AIC}
# models_ibex$AIC <- as.numeric(models_ibex$AIC)
# models_ibex$QAIC <- as.numeric(models_ibex$QAIC)
#  
# #write.table(models_ibex, "models_ibex.csv",sep=";",dec=".",row.names=F,col.names=T)
#  
# #### Create best model based on AIC / QAIC ####
# 
# best_model_QAIC_ibex <- models_ibex[models_ibex$QAIC == min(models_ibex$QAIC),]
#  
# fm_ibex_best_QAIC <- colext(psiformula = as.formula(best_model_QAIC_ibex$psi_formula), gammaformula = as.formula(best_model_QAIC_ibex$gamma_formula), epsilonformula = as.formula(best_model_QAIC_ibex$epsilon_formula), pformula = as.formula(best_model_QAIC_ibex$p_formula), umfm_ibex)
# summary(fm_ibex_best_QAIC)

fm_ibex_best_QAIC <- colext(~1, ~scale(OakForest), ~1, ~scale(Sensitivity)+scale(Rad_max)+scale(LowShrub)+scale(OakForest), umfm_ibex)
summary(fm_ibex_best_QAIC)

#### Backtransform model estimates from logit scale (keeping Intercept at 1) ####

variables_ibex <- data.frame(matrix(ncol = 10, nrow = 10))
colnames(variables_ibex) <- c("Species","Variable","Occupancy","Occ_SE","Colonization","Col_SE","Extinction","Ext_SE","Detection","Det_SE")
variables_ibex$Species <- "Iberian ibex"
variables_ibex$Variable <- c("Intercept", "Sensitivity","Rad_max","Temp_min","CostDist_road","CostDist_water","OakForest","HighShrub","LowShrub","Urban")

variables_ibex[1,3] <- backTransform(fm_ibex_best_QAIC,type="psi")@estimate
variables_ibex[1,4] <- SE(backTransform(fm_ibex_best_QAIC,type="psi"))

variables_ibex[1,5] <- backTransform(linearComb(fm_ibex_best_QAIC,c(1,0),type="col"))@estimate
variables_ibex[1,6] <- SE(backTransform(linearComb(fm_ibex_best_QAIC,c(1,0),type="col")))
variables_ibex[7,5] <- backTransform(linearComb(fm_ibex_best_QAIC,c(1,1),type="col"))@estimate
variables_ibex[7,6] <- SE(backTransform(linearComb(fm_ibex_best_QAIC,c(1,1),type="col")))

variables_ibex[1,7] <- backTransform(fm_ibex_best_QAIC,type="ext")@estimate
variables_ibex[1,8] <- SE(backTransform(fm_ibex_best_QAIC,type="ext"))

variables_ibex[1,9] <- backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,0,0,0),type="det"))@estimate
variables_ibex[1,10] <- SE(backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,0,0,0),type="det")))
variables_ibex[2,9] <- backTransform(linearComb(fm_ibex_best_QAIC,c(1,1,0,0,0),type="det"))@estimate
variables_ibex[2,10] <- SE(backTransform(linearComb(fm_ibex_best_QAIC,c(1,1,0,0,0),type="det")))
variables_ibex[3,9] <- backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,1,0,0),type="det"))@estimate
variables_ibex[3,10] <- SE(backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,1,0,0),type="det")))
variables_ibex[9,9] <- backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,0,1,0),type="det"))@estimate
variables_ibex[9,10] <- SE(backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,0,1,0),type="det")))
variables_ibex[7,9] <- backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,0,0,1),type="det"))@estimate
variables_ibex[7,10] <- SE(backTransform(linearComb(fm_ibex_best_QAIC,c(1,0,0,0,1),type="det")))

m3 <- nonparboot(fm_ibex_best_QAIC,B = 1000)
occ_pred_ibex_best_QAIC <- data.frame(year = c(2015:2022),
                                     smoothed_occ = smoothed(fm_ibex_best_QAIC)[2,],
                                     SE = m3@smoothed.mean.bsse[2,])
occ_pred_ibex_best_QAIC$Species <- 'Iberian ibex'

ibex_year <- lm(smoothed_occ ~ year, occ_pred_ibex_best_QAIC)
summary(ibex_year)

#### Get site-specific occupancy estimates obtained from model output ####

occ_ibex_site <- data.frame(t(fm_ibex_best_QAIC@smoothed[2, , ]))
colnames(occ_ibex_site) <- unique(yearlySiteCovs(umfm_ibex)$year)
occ_ibex_site$grid <- unique(Ibex$grid) 

#write.table(occ_ibex_site, "Occ_ibex_site.csv",sep=";",dec=".",row.names=F,col.names=T)

occ_ibex_site_long <- data.frame(melt(occ_ibex_site, id.vars="grid"))
names(occ_ibex_site_long) <- c("Grid","Year","Occupancy")
occ_ibex_site_long$Year <- as.numeric(as.character(occ_ibex_site_long$Year))
occ_ibex_site_long$Species <- "Iberian ibex"

######## Create table for all species ####### 

## Occupancy estimates for each species per year
occ_pred_all <- rbind(occ_pred_cattle_best_QAIC,occ_pred_horse_best_QAIC,occ_pred_deer_best_QAIC,occ_pred_boar_best_QAIC,occ_pred_fox_best_QAIC,occ_pred_wolf_best_QAIC,occ_pred_ibex_best_QAIC)
occ_pred_all$Species <- as.factor(occ_pred_all$Species)
occ_pred_all$Species <- factor(occ_pred_all$Species, c("Domestic cattle", "Domestic horse", "European roe deer","Wild boar","Red fox","Gray wolf","Iberian ibex"))
write.table(occ_pred_all, "occ_pred_all.csv",sep=";",dec=".",row.names=F,col.names=T)

## Long format
occ_all_site_long <- rbind(occ_cattle_site_long,occ_horse_site_long,occ_deer_site_long,occ_boar_site_long,occ_fox_site_long,occ_wolf_site_long,occ_ibex_site_long)
occ_all_site_long$Species <- as.factor(occ_all_site_long$Species)
occ_all_site_long$Species <- factor(occ_all_site_long$Species, c("Domestic cattle", "Domestic horse", "European roe deer","Wild boar","Red fox","Gray wolf","Iberian ibex"))
write.table(occ_all_site_long, "occ_all_site_long.csv",sep=";",dec=".",row.names=F,col.names=T)

## Occupancy estimates for each species per grid cell and year
occ_all_grid <- cbind(occ_cattle_site_long[-4],occ_horse_site_long[3],occ_deer_site_long[3],occ_boar_site_long[3],occ_fox_site_long[3],occ_wolf_site_long[3],occ_ibex_site_long[3])
names(occ_all_grid) <- c("Grid","Year","Cattle","Horse","RoeDeer","WildBoar","Fox","Wolf","Ibex")
write.table(occ_all_grid, "occ_all_grid.csv",sep=";",dec=".",row.names=F,col.names=T)

## Effect of variables for each species
var_effect <- rbind(variables_cattle,variables_horse,variables_deer,variables_boar,variables_fox,variables_wolf,variables_ibex)
write.table(var_effect, "var_effect.csv",sep=";",dec=".",row.names=F,col.names=T)

###### Interaction with other species ####### 

#### Test for possible interactions with a linear regression ####

cattle_int<- lm(Cattle ~ Horse+RoeDeer+WildBoar+Fox+Wolf+Ibex, occ_all_grid)
summary(cattle_int)

cattle_int_final <- stepAIC(cattle_int)
summary(cattle_int_final)

####

horse_int <- lm(Horse ~ Cattle+RoeDeer+WildBoar+Fox+Wolf+Ibex, occ_all_grid)
summary(horse_int)

horse_int_final <- stepAIC(horse_int)
summary(horse_int_final)

####

deer_int <- lm(RoeDeer ~ Cattle+Horse+WildBoar+Fox+Wolf+Ibex, occ_all_grid)
summary(deer_int)

deer_int_final <- stepAIC(deer_int)
summary(deer_int_final)

####

boar_int <- lm(WildBoar ~ Cattle+Horse+RoeDeer+Fox+Wolf+Ibex, occ_all_grid)
summary(boar_int)

boar_int_final <- stepAIC(boar_int)
summary(boar_int_final)

####

fox_int <- lm(Fox ~ Cattle+Horse+RoeDeer+WildBoar+Wolf+Ibex, occ_all_grid)
summary(fox_int)

fox_int_final <- stepAIC(fox_int)
summary(fox_int_final)

####

wolf_int <- lm(Wolf ~ Cattle+Horse+RoeDeer+WildBoar+Fox+Ibex, occ_all_grid)
summary(wolf_int)

wolf_int_final <- stepAIC(wolf_int)
summary(wolf_int_final)

####

ibex_int <- lm(Ibex ~ Cattle+Horse+RoeDeer+WildBoar+Fox+Wolf, occ_all_grid)
summary(ibex_int)

ibex_int_final <- stepAIC(ibex_int)
summary(ibex_int_final)

############################################## Create figures ########################################################

setwd(dir = file.path(t,"Dynamics_Coexistence_Peneda-main/Figures"))

###### Figure 3. Occupancy trends ######

occ_all_mod_legend <- ggplot(data=occ_pred_all, aes(x=year, y=smoothed_occ)) +
  geom_ribbon(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE, fill=Species),alpha=0.3) +
  geom_line(aes(x=year, y=smoothed_occ, col=Species),size=1) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=12))

get_only_legend <- function(plot) { 
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  legend <- plot_table$grobs[[legend_plot]] 
  return(legend) 
} 

legend <- get_only_legend(occ_all_mod_legend)

occ_all_mod <- ggplot(data=occ_pred_all, aes(x=year, y=smoothed_occ)) +
  geom_ribbon(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE, fill=Species),alpha=0.3) +
  geom_line(aes(x=year, y=smoothed_occ, col=Species),size=1) +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  labs(title="Occupancy",x="Year", y="Smoothed occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), 
        plot.title=element_text(size=16,face="bold"),
        axis.title.x= element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(1.2, 'cm'),
        legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=12),
        legend.position="none")

occ_all_lm <- ggplot(occ_all_site_long, aes(Year, Occupancy)) +
  geom_point(aes(x=Year, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Year, y=Occupancy, col=Species,fill=Species),alpha=0.3,size=1,method="lm") +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  labs(title="",x="Year", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        plot.title = element_blank(),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.position="none")

png(file="Figure_3_Occupancy_trends.png", width=800, height=1000)
grid.arrange(occ_all_mod, occ_all_lm, legend, layout_matrix = rbind(c(1, 3),c(2, 3)),widths=c(1/2,1/6))
dev.off()

###### Figure 4. Effect sizes of the variables  ####### 

## Occupancy

occ_effect <- var_effect[1:4]
occ_effect <- occ_effect[which(occ_effect$Variable != "Sensitivity" & occ_effect$Variable != "Rad_max" & occ_effect$Variable != "Temp_min"),]

occ_effect$var <- as.character(interaction(occ_effect$Species,occ_effect$Variable))
occ_effect$var <- factor(occ_effect$var, levels = occ_effect$var)
occ_effect$Variable <- as.factor(occ_effect$Variable)
occ_effect$Variable <- factor(occ_effect$Variable, c("Intercept", "CostDist_road", "CostDist_water","OakForest","HighShrub","LowShrub","Urban"))

point_shape <- rep(c(15,16,17,18,0,1,2),7)
point_col <- c(rep("#1B9E77",7),rep("#D95F02",7),rep("#7570B3",7),rep("#A6761D",7),rep("#E6AB02",7),rep("#E7298A",7),rep("#66A61E",7))

occupancy_plot <- ggplot(aes(y = var),data=occ_effect) + 
  geom_point(aes(x=Occupancy,shape=Variable), size=5,color=point_col) +
  scale_shape_manual(values = point_shape) +
  geom_linerange(aes(xmin=Occupancy-Occ_SE, xmax=Occupancy+Occ_SE),color=point_col) +
  geom_hline(yintercept = c(7.5,14.5,21.5,28.5,35.5,42.5)) +
  geom_segment(aes(x = 0.845, y = 49.5, xend = 0.845, yend = 42.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 0.879, y = 42.5, xend = 0.879, yend = 35.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 1, y = 35.5, xend = 1, yend = 28.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 0.914, y = 28.5, xend = 0.914, yend = 21.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 0.467, y = 21.5, xend = 0.467, yend = 14.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 0.537, y = 14.5, xend = 0.537, yend = 7.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.000714, y = 7.5, xend = 0.000714, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev, labels=rev(occ_effect$Variable)) +
  xlim(-0.1,1) +
  ggtitle("Initial occupancy") +
  theme_bw() +
  theme(legend.position = c(0.13,0.81),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=15),
        plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.line = element_line( linetype = "solid"))

## Colonization

col_effect <- var_effect[c(1,2,5,6)]
col_effect <- col_effect[which(col_effect$Variable != "Sensitivity" & col_effect$Variable != "Rad_max" & col_effect$Variable != "Temp_min"),]

col_effect$var <- as.character(interaction(col_effect$Species,col_effect$Variable))
col_effect$var <- factor(col_effect$var, levels = col_effect$var)
col_effect$Variable <- as.factor(col_effect$Variable)
col_effect$Variable <- factor(col_effect$Variable, c("Intercept", "CostDist_road", "CostDist_water","OakForest","HighShrub","LowShrub","Urban"))

colonization_plot <- ggplot(aes(y = var),data=col_effect) + 
  geom_point(aes(x=Colonization,shape=Variable), size=5,color=point_col) +
  scale_shape_manual(values = point_shape) +
  geom_linerange(aes(xmin=Colonization-Col_SE, xmax=Colonization+Col_SE),color=point_col) +
  geom_hline(yintercept = c(7.5,14.5,21.5,28.5,35.5,42.5)) +
  geom_segment(aes(x = 0.398, y = 49.5, xend = 0.398, yend = 42.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 0.578, y = 42.5, xend = 0.578, yend = 35.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 0.528, y = 35.5, xend = 0.528, yend = 28.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 0.837, y = 28.5, xend = 0.837, yend = 21.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 0.876, y = 21.5, xend = 0.876, yend = 14.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 0.33, y = 14.5, xend = 0.33, yend = 7.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.0941, y = 7.5, xend = 0.0941, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev, labels=rev(col_effect$Variable)) +
  xlim(-0.1,1) +
  ggtitle("Colonization probability") +
  theme_bw() +
  theme(legend.position = c(0.13,0.81),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=15),
        plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.line = element_line(linetype = "solid"))

## Extinction

ext_effect <- var_effect[c(1,2,7,8)]
ext_effect <- ext_effect[which(ext_effect$Variable != "Sensitivity" & ext_effect$Variable != "Rad_max" & ext_effect$Variable != "Temp_min"),]

ext_effect$var <- as.character(interaction(ext_effect$Species,ext_effect$Variable))
ext_effect$var <- factor(ext_effect$var, levels = ext_effect$var)
ext_effect$Variable <- as.factor(ext_effect$Variable)
ext_effect$Variable <- factor(ext_effect$Variable, c("Intercept", "CostDist_road", "CostDist_water","OakForest","HighShrub","LowShrub","Urban"))

extinction_plot <- ggplot(aes(y = var),data=ext_effect) + 
  geom_point(aes(x=Extinction,shape=Variable), size=5,color=point_col) +
  scale_shape_manual(values = point_shape) +
  geom_linerange(aes(xmin=Extinction-Ext_SE, xmax=Extinction+Ext_SE),color=point_col) +
  geom_hline(yintercept = c(7.5,14.5,21.5,28.5,35.5,42.5)) +
  geom_segment(aes(x = 0.329, y = 49.5, xend = 0.329, yend = 42.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 0.115, y = 42.5, xend = 0.115, yend = 35.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 0.0157, y = 35.5, xend = 0.0157, yend = 28.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 0.102, y = 28.5, xend = 0.102, yend = 21.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 0.418, y = 21.5, xend = 0.418, yend = 14.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 0.0864, y = 14.5, xend = 0.0864, yend = 7.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.000277, y = 7.5, xend = 0.000277, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev, labels=rev(ext_effect$Variable)) +
  xlim(-0.1,1) +
  ggtitle("Extinction probability") +
  theme_bw() +
  theme(legend.position = c(0.875,0.81),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=15),
        plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.line = element_line(linetype = "solid"))

## Detection

det_effect <- var_effect[c(1,2,9,10)]
det_effect <- det_effect[which(det_effect$Variable != "CostDist_road" & det_effect$Variable != "CostDist_water" & det_effect$Variable != "Urban"),]

det_effect$var <- as.character(interaction(det_effect$Species,det_effect$Variable))
det_effect$var <- factor(det_effect$var, levels = det_effect$var)
det_effect$Variable <- factor(det_effect$Variable, c("Intercept", "Sensitivity", "Rad_max","Temp_min","OakForest","HighShrub","LowShrub"))

detection_plot <- ggplot(aes(y = var),data=det_effect) + 
  geom_point(aes(x=Detection,shape=Variable), size=5,color=point_col) +
  scale_shape_manual(values = point_shape) +
  geom_linerange(aes(xmin=Detection-Det_SE, xmax=Detection+Det_SE),color=point_col) +
  geom_hline(yintercept = c(7.5,14.5,21.5,28.5,35.5,42.5)) +
  geom_segment(aes(x = 0.282, y = 49.5, xend = 0.282, yend = 42.5),col="#1B9E77",linetype="dashed") +
  geom_segment(aes(x = 0.304, y = 42.5, xend = 0.304, yend = 35.5),col="#D95F02",linetype="dashed") +
  geom_segment(aes(x = 0.386, y = 35.5, xend = 0.386, yend = 28.5),col="#7570B3",linetype="dashed") +
  geom_segment(aes(x = 0.289, y = 28.5, xend = 0.289, yend = 21.5),col="#A6761D",linetype="dashed") +
  geom_segment(aes(x = 0.114, y = 21.5, xend = 0.114, yend = 14.5),col="#E6AB02",linetype="dashed") +
  geom_segment(aes(x = 0.0483, y = 14.5, xend = 0.0483, yend = 7.5),col="#E7298A",linetype="dashed") +
  geom_segment(aes(x = 0.0167, y = 7.5, xend = 0.0167, yend = 0.5),col="#66A61E",linetype="dashed") +
  scale_y_discrete(limits=rev) +
  xlim(-0.1,1) +
  ggtitle("Detection probability") +
  theme_bw() +
  theme(legend.position = c(0.9,0.81),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text.x = element_text(size=15),
        plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.line = element_line(linetype = "solid"))

var_effect$Species <- as.factor(var_effect$Species)
var_effect$Species <- factor(var_effect$Species, c("Domestic cattle", "Domestic horse", "European roe deer","Wild boar","Red fox","Gray wolf","Iberian ibex"))

plot_legend <-ggplot(var_effect, aes(y=Occupancy, fill=Species)) +
  geom_bar() +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=18,face="bold"),
        legend.text = element_text(size=15))

legend <- get_only_legend(plot_legend) 

png(file="Figure_4_Variable_effects.png", width=1100, height=1500)
grid.arrange(occupancy_plot, colonization_plot, extinction_plot, detection_plot, legend, layout_matrix = rbind(c(1, 2),c(3, 4),c(5,5)),heights = c(1/2, 1/2, 1/9))
dev.off()

###### Figure 5. Effect of the occurrence of A) domestic cattle, B) domestic horses and C) wolf  ####### 

cattle <- melt(occ_all_grid, id.vars=c("Grid", "Year","Cattle"))
names(cattle) <- c("Grid","Year","Cattle","Species","Occupancy")

cattle_int <- ggplot(cattle, aes(Cattle, Occupancy)) +
  geom_point(aes(x=Cattle, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Cattle, y=Occupancy, col=Species,fill=Species),alpha=0.3,size=0.7,method="lm") +
  scale_fill_manual(values=c("#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#D95F02","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  labs(title="Occupancy estimates",x="Cattle occupancy", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        plot.title = element_blank(),
        axis.text=element_text(size=15), 
        axis.title=element_text(size=16),
        legend.position="none")

horse <- melt(occ_all_grid, id.vars=c("Grid", "Year","Horse"))
names(horse) <- c("Grid","Year","Horse","Species","Occupancy")

horse_int <- ggplot(horse, aes(Horse, Occupancy)) +
  geom_point(aes(x=Horse, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Horse, y=Occupancy, col=Species,fill=Species),alpha=0.3,size=0.7,method="lm") +
  scale_fill_manual(values=c("#1B9E77","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#7570B3","#A6761D","#E6AB02","#E7298A","#66A61E"))  +
  labs(title="Occupancy estimates",x="Horse occupancy", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        plot.title = element_blank(),
        axis.title.y= element_blank(),
        axis.text=element_text(size=15), 
        axis.title=element_text(size=16),
        legend.position="none")

wolf <- melt(occ_all_grid, id.vars=c("Grid", "Year","Wolf"))
names(wolf) <- c("Grid","Year","Wolf","Species","Occupancy")

wolf_int <- ggplot(wolf, aes(Wolf, Occupancy)) +
  geom_point(aes(x=Wolf, y=Occupancy, col=Species)) +
  geom_smooth(aes(x=Wolf, y=Occupancy, col=Species,fill=Species),alpha=0.3,size=0.7,method="lm") +
  scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#66A61E")) +
  scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3","#A6761D","#E6AB02","#66A61E"))  +
  labs(title="Occupancy estimates",x="Wolf occupancy", y="Occupancy") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        plot.title = element_blank(),
        axis.title.y= element_blank(),
        axis.text=element_text(size=15), 
        axis.title=element_text(size=16),
        legend.position="none")

legend <- get_only_legend(occ_all_mod_legend)

png(file="Figure_5_Interaction_plot.png", width=1600, height=600)
grid.arrange(cattle_int, horse_int, wolf_int, legend, nrow=1,widths=c(1/3,1/3,1/3,1/7))
dev.off()


