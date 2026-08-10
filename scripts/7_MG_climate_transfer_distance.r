library(dplyr)
library(ggsci)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(mgcv)


setwd('/home/kati/Dropbox/projects/MyGardenOfTrees_trials/Manuscripts/MGOT_PilotTrials/analysis/gardens/')

mle.a <- read.table("Abies_MLE.txt", header = T)
mle.f <- read.table("Fagus_MLE.txt", head=T)

## combine results with seed data
mle.a$Genus = "Abies"
mle.f$Genus = "Fagus"
mle.dat = rbind(mle.a, mle.f)

load("dat2022_germ.RData")
prov.a.keep = unique(mle.a$provenance)[unique(mle.a$provenance) %in% unique(datA$Provenance)]
mle.a = subset(mle.a, provenance %in% prov.a.keep)

## MLE was done for the two iranian prov combined: add same MLE estimate for both
tmp = mle.f[mle.f$provenance == "Iran Alborz Mountains", ]
tmp$provenance = "Iran Alborz Mountains High"
mle.f[mle.f$provenance == "Iran Alborz Mountains", "provenance"] = "Iran Alborz Mountains Middle"
mle.f = rbind(mle.f, tmp)

prov.f.keep = unique(mle.f$provenance)[unique(mle.f$provenance) %in% unique(datF$Provenance)]
mle.f = subset(mle.f, provenance %in% prov.f.keep)

mle.a$garden = gsub("21_21", "21_", mle.a$garden)
mle.a$garden = gsub("0", "", mle.a$garden)
mle.f$garden = gsub("21_21", "21_", mle.f$garden)
mle.f$garden = gsub("0", "", mle.f$garden)

## merge MLE and germ + env datasets
datF$garden  = as.character(datF$Garden_ID)
datA$garden  = as.character(datA$Garden_ID)
datF$provenance = datF$Provenance
datA$provenance = datA$Provenance
datA = merge(mle.a, datA, by=c("provenance", "garden"), all=T)
datF = merge(mle.f, datF, by=c("provenance", "garden"), all=T)
dat = rbind(datA, datF)

## exclude species
datA.aa = subset(datA, provenance != "Georgia Ambrolauri")
datF.fs = subset(datF, !provenance %in% c("Iran Alborz Mountains High", "Iran Alborz Mountains Middle"))


## env distances
## ###########################
datA$nitrogen_dist = datA$nitrogen_15.30cm - datA$nitrogen_15.30cm_garden
datA$clay_dist = datA$clay_15.30cm - datA$clay_15.30cm_garden

datF$nitrogen_dist = datF$nitrogen_15.30cm - datF$nitrogen_15.30cm_garden
datF$clay_dist = datF$clay_15.30cm - datF$clay_15.30cm_garden

datA$Prcp_2022_spring <- apply(datA[,c("Prcp_2022_3", "Prcp_2022_4", "Prcp_2022_5")], 1, sum)
datF$Prcp_2022_spring <- apply(datF[,c("Prcp_2022_3", "Prcp_2022_4", "Prcp_2022_5")], 1, sum)
datA$Prcp_spring_prov <- apply(datA[,c("Prcp_seed_prov_Mar", "Prcp_seed_prov_Apr", "Prcp_seed_prov_May")], 1, sum) 
datF$Prcp_spring_prov <- apply(datF[,c("Prcp_seed_prov_Mar", "Prcp_seed_prov_Apr", "Prcp_seed_prov_May")], 1, sum) 
datA$Prcp_dist <- datA$Prcp_spring_prov - datA$Prcp_2022_spring
datF$Prcp_dist <- datF$Prcp_spring_prov - datF$Prcp_2022_spring

datA$T_2022_spring <- apply(datA[,c("T_mean_2022_3", "T_mean_2022_4", "T_mean_2022_5")], 1, sum)
datF$T_2022_spring <- apply(datF[,c("T_mean_2022_3", "T_mean_2022_4", "T_mean_2022_5")], 1, sum)
datA$T_spring_prov <- apply(datA[,c("T_mean_seed_prov_Mar", "T_mean_seed_prov_Apr", "T_mean_seed_prov_May")], 1, sum) 
datF$T_spring_prov <- apply(datF[,c("T_mean_seed_prov_Mar", "T_mean_seed_prov_Apr", "T_mean_seed_prov_May")], 1, sum) 
datA$T_dist <- datA$T_spring_prov - datA$T_2022_spring
datF$T_dist <- datF$T_spring_prov - datF$T_2022_spring

##Prcp_distance_April_May
##T_mean_distance_April_May

modA <-gam(
  g_without_death ~
      s(Longitude_distance, k = 4) +
      s(Latitude_distance, k = 4) +
      s(T_dist, k = 4) +
      s(Prcp_dist, k = 4) +
      s(clay_dist, k = 4) +
      s(nitrogen_dist, k = 4) +
      s(ID, bs = "re") +
      s(Garden_ID, bs = "re"),
  data = datA,
  method = "REML"
)


modF <- gam(
  g_without_death ~
      s(Longitude_distance, k = 4) +
      s(Latitude_distance, k = 4) +
      s(T_dist, k = 4) +
      s(Prcp_dist, k = 4) +
      s(clay_dist, k = 4) +
      s(nitrogen_dist, k = 4) +
      s(ID, bs = "re") + ## each garden has its own intercept
      s(Garden_ID, bs = "re"), ## each prov has its own intercept
  data = datF,
  method = "REML"
)


get_p_label <- function(mod, term_label) {
  s <- summary(mod)$s.table
  p <- s[term_label, "p-value"]
  
  if (is.na(p)) return("")
  if (p >= 0.05) return("NS")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", signif(p, 2))
}

plot_fullgam_effect <- function(dat, mod, xvar, yvar = "g_without_death",
                                main = "", xlab = "", ylab = "Germination rate",
                                term_label = NULL, add_legend = FALSE,
                                k = 200) {
  
  dat <- dat[is.finite(dat[[xvar]]) & is.finite(dat[[yvar]]), ]
  
  plot(dat[[xvar]], dat[[yvar]],
       col = adjustcolor(dat$col, alpha.f = 0.75),
       pch = 20, cex = 0.8,
       xlab = xlab, ylab = ylab, main = main)
  
  nd <- data.frame(
    Longitude_distance = 0,
    Latitude_distance  = 0,
    T_dist             = 0,
    Prcp_dist          = 0,
    clay_dist          = 0,
    nitrogen_dist      = 0,
    ID                 = dat$ID[1],
    Garden_ID          = dat$Garden_ID[1]
  )
  
  nd <- nd[rep(1, k), ]
  nd[[xvar]] <- seq(min(dat[[xvar]], na.rm = TRUE),
                    max(dat[[xvar]], na.rm = TRUE),
                    length.out = k)
  
  pr <- predict(mod, newdata = nd, se.fit = TRUE,
                exclude = c("s(ID)", "s(Garden_ID)"))
  
  lines(nd[[xvar]], pr$fit, lwd = 2)
  lines(nd[[xvar]], pr$fit + 2 * pr$se.fit, lty = 2)
  lines(nd[[xvar]], pr$fit - 2 * pr$se.fit, lty = 2)
  
  if (!is.null(term_label)) {
    lab <- get_p_label(mod, term_label)
    usr <- par("usr")
    text(usr[1] + 0.03 * diff(usr[1:2]),
         usr[4] - 0.08 * diff(usr[3:4]),
         lab, adj = c(0, 1), cex = 0.9)
  }
  
  if (add_legend) {
    leg <- dat[!duplicated(dat$ID), c("ID", "col")]
    legend("topright",
           legend = leg$ID,
           col = leg$col,
           pch = 20,
           cex = 0.65,
           bty = "o",
           box.col = NA,
           bg = adjustcolor("white", alpha.f = 0.7))
  }
}

pdf("SFig_Transfer_distance.pdf", 5, 10)

par(mfcol = c(6, 2),
    mar = c(3, 3, 1, .5),
    mgp = c(2, 0.7, 0))

plot_fullgam_effect(datA, modA, "Latitude_distance",  main = "Abies", xlab = "Latitude transfer distance", term_label = "s(Latitude_distance)", add_legend = TRUE)
plot_fullgam_effect(datA, modA, "Longitude_distance", main = "", xlab = "Longitude transfer distance", term_label = "s(Longitude_distance)")
plot_fullgam_effect(datA, modA, "T_dist",             main = "", xlab = "Spring temperature transfer distance", term_label = "s(T_dist)")
plot_fullgam_effect(datA, modA, "Prcp_dist",          main = "", xlab = "Spring precipitation transfer distance", term_label = "s(Prcp_dist)")
plot_fullgam_effect(datA, modA, "nitrogen_dist",      main = "", xlab = "Soil nitrogen transfer distance", term_label = "s(nitrogen_dist)")
plot_fullgam_effect(datA, modA, "clay_dist",          main = "", xlab = "Soil clay transfer distance", term_label = "s(clay_dist)")

plot_fullgam_effect(datF, modF, "Latitude_distance",  main = "Fagus", xlab = "Latitude transfer distance", ylab = "", term_label = "s(Latitude_distance)", add_legend = TRUE)
plot_fullgam_effect(datF, modF, "Longitude_distance", main = "", xlab = "Longitude transfer distance", ylab = "", term_label = "s(Longitude_distance)")
plot_fullgam_effect(datF, modF, "T_dist",             main = "", xlab = "Spring temperature transfer distance", ylab = "", term_label = "s(T_dist)")
plot_fullgam_effect(datF, modF, "Prcp_dist",          main = "", xlab = "Spring precipitation transfer distance", ylab = "", term_label = "s(Prcp_dist)")
plot_fullgam_effect(datF, modF, "nitrogen_dist",      main = "", xlab = "Soil nitrogen transfer distance", ylab = "", term_label = "s(nitrogen_dist)")
plot_fullgam_effect(datF, modF, "clay_dist",          main = "", xlab = "Soil clay transfer distance", ylab = "", term_label = "s(clay_dist)")

dev.off()
