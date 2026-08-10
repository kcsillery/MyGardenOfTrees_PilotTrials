library(asreml)
##asreml.license.activate()
## HCEE-JCDA-JBFD-BJFA
library(lubridate)
library(RColorBrewer)

library(ggplot2)
library(patchwork)
setwd('/home/kati/Dropbox/projects/MyGardenOfTrees_trials/Manuscripts/MGOT_PilotTrials/analysis/gardens/')
load("dat2022_germ.RData") ## see 2_lasso_varselect.r
load(file="best_lasso_variables.RData")

datA$Prcp_2022_spring <- apply(datA[,c("Prcp_2022_3", "Prcp_2022_4", "Prcp_2022_5")], 1, sum)
datF$Prcp_2022_spring <- apply(datF[,c("Prcp_2022_3", "Prcp_2022_4", "Prcp_2022_5")], 1, sum)
datA$T_mean_2022_spring <- apply(datA[,c("T_mean_2022_3", "T_mean_2022_4", "T_mean_2022_5")], 1, mean)
datF$T_mean_2022_spring <- apply(datF[,c("T_mean_2022_3", "T_mean_2022_4", "T_mean_2022_5")], 1, mean)

## keep only lasso selected climate variables + seed variables
seed_vars <- c("Seed_weight_g", "Moisture_content_percent")
env_vars <- c("Prcp_2022_spring", "T_mean_2022_spring", "clay_15.30cm_garden", "silt_15.30cm_garden", "nitrogen_15.30cm_garden", "bio12_1901_1980_seed_prov")
modelvars <- c("Genus", "ID", "Garden_ID", "GardenBlock", "Garden_block_spot", "Date", "Date_int", "DOY", "germ", "germ.max.all")
datA <- datA[, match(c(modelvars, seed_vars, env_vars), names(datA))]
datF <- datF[, match(c(modelvars, seed_vars, env_vars), names(datF))]

## scale variables
datA[, c(seed_vars, env_vars)] <- apply(datA[, c(seed_vars, env_vars)], 2, scale)
datF[, c(seed_vars, env_vars)] <- apply(datF[, c(seed_vars, env_vars)], 2, scale)

datA$GardenBlock <- as.factor(datA$GardenBlock)
datF$GardenBlock <- as.factor(datF$GardenBlock)

## fixed effects
myfixed_prcp <- germ.max.all ~ ID + Prcp_2022_spring + T_mean_2022_spring + clay_15.30cm_garden + silt_15.30cm_garden + nitrogen_15.30cm_garden + ID:Prcp_2022_spring
myfixed_temp <- germ.max.all ~ ID + Prcp_2022_spring + T_mean_2022_spring + clay_15.30cm_garden + silt_15.30cm_garden + nitrogen_15.30cm_garden + ID:T_mean_2022_spring
myfixed_clay <- germ.max.all ~ ID + Prcp_2022_spring + T_mean_2022_spring + clay_15.30cm_garden + silt_15.30cm_garden + nitrogen_15.30cm_garden + ID:clay_15.30cm_garden
myfixed_silt <- germ.max.all ~ ID + Prcp_2022_spring + T_mean_2022_spring + clay_15.30cm_garden + silt_15.30cm_garden + nitrogen_15.30cm_garden + ID:silt_15.30cm_garden
myfixed_N <- germ.max.all ~ ID + Prcp_2022_spring + T_mean_2022_spring + clay_15.30cm_garden + silt_15.30cm_garden + nitrogen_15.30cm_garden + ID:nitrogen_15.30cm_garden


## ##################################
##             Abies
## ##################################

## precipitation
## ###############
modA.prcp.fixedslopes <- asreml(
    fixed = myfixed_prcp,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datA)

summary(modA.prcp.fixedslopes)$varcomp
wald.asreml(modA.prcp.fixedslopes, ssType = "conditional")


## temperature
## ###############
modA.temp.fixedslopes <- asreml(
    fixed = myfixed_temp,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datA)

summary(modA.temp.fixedslopes)$varcomp
wald.asreml(modA.temp.fixedslopes, ssType = "conditional")


## clay
## ###############
modA.clay.fixedslopes <- asreml(
    fixed = myfixed_clay,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datA)

summary(modA.clay.fixedslopes)$varcomp
wald.asreml(modA.clay.fixedslopes, ssType = "conditional")


## silt
## ###############
modA.silt.fixedslopes <- asreml(
    fixed = myfixed_silt,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datA)

summary(modA.silt.fixedslopes)$varcomp
wald.asreml(modA.silt.fixedslopes, ssType = "conditional")


## nitrogen
## ###############
modA.N.fixedslopes <- asreml(
    fixed = myfixed_N,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datA)

summary(modA.N.fixedslopes)$varcomp
wald.asreml(modA.N.fixedslopes, ssType = "conditional")





# ##################################
##             Fagus
## ##################################

## precipitation
## ###############
modF.prcp.fixedslopes <- asreml(
    fixed = myfixed_prcp,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datF)

summary(modF.prcp.fixedslopes)$varcomp
wald.asreml(modF.prcp.fixedslopes, ssType = "conditional")


## temperature
## ###############
modF.temp.fixedslopes <- asreml(
    fixed = myfixed_temp,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datF)

summary(modF.temp.fixedslopes)$varcomp
wald.asreml(modF.temp.fixedslopes, ssType = "conditional")


## clay
## ###############
modF.clay.fixedslopes <- asreml(
    fixed = myfixed_clay,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datF)

summary(modF.clay.fixedslopes)$varcomp
wald.asreml(modF.clay.fixedslopes, ssType = "conditional")


## silt
## ###############
modF.silt.fixedslopes <- asreml(
    fixed = myfixed_silt,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datF)

summary(modF.silt.fixedslopes)$varcomp
wald.asreml(modF.silt.fixedslopes, ssType = "conditional")


## nitrogen
## ###############
modF.N.fixedslopes <- asreml(
    fixed = myfixed_N,
    random = ~ Garden_ID + GardenBlock + ID:Garden_ID + ID:GardenBlock,
    family = asr_negative.binomial(dispersion = NA),
    na.action = na.method(x = ""),
    maxiter = 200,
    ai.sing = TRUE,
    data = datF)

summary(modF.N.fixedslopes)$varcomp
wald.asreml(modF.N.fixedslopes, ssType = "conditional")



