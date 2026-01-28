## ##########################################################################

## Katalin Csillery

## 28 Jan 2026

## For the manuscript:
## Gene–environment interactions govern early regeneration in fir and
## beech: evidence from participatory provenance trials across Europe
## by Katalin Csilléry, Justine Charlet de Sauvage, Madleina Caduff,
## Johannes Alt, Marjorie Bison, Mert Celik, Nicole Ponta, Daniel
## Wegmann

## ###########################################################################


library(lubridate)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(asreml)

load(file="dat.RData")

## ######################################################
##                     Abies
## ######################################################

datA <- subset(dat.all, Genus=="Abies")
datA <- subset(datA, ! ID %in% c("AA_HR1","AA_IT1"))
##datA <- subset(datA, ! Garden_ID %in% c(16, 5, 3))

## make factors
datA$Garden_ID <- as.factor(datA$Garden_ID)
datA$Block <- as.factor(datA$Block)
datA$Garden_block_spot <- as.factor(datA$Garden_block_spot)
datA$ID <- as.factor(datA$ID)
datA$Year <- format(datA$Date, "%Y")

## order by Garden_block_spot, then Date
datA <- datA[order(datA$Garden_ID, datA$Garden_block_spot, datA$Date), ]

## change 2025 entries to 2024 last date
last.date.2024 <- max(datA[datA$Year==2024, "Date"])
datA$Date[datA$Year==2025] <- last.date.2024
datA <- subset(datA, Year<2025)

## cohort aggregate
datAsum <- aggregate(Numind ~ Garden_ID + ID + Date,
                      data = datA, FUN = sum)
datAsum$Year <- format(datAsum$Date, "%Y")

## Define N0 = max aggregated count in 2022 (max germination rate)
datAsum$N0 <- NA_real_
uid <- unique(paste(datAsum$Garden_ID, datAsum$ID, sep="|"))
for (u in uid) {
  i <- paste(datAsum$Garden_ID, datAsum$ID, sep="|") == u
  j <- i & datAsum$Year == "2022"
  if (any(j)) datAsum$N0[i] <- max(datAsum$Numind[j], na.rm=TRUE)
}
datAsum <- datAsum[datAsum$N0 > 0, ] ## drop cohorts never germinated

## order
datAsum <- datAsum[order(datAsum$Garden_ID,
                           datAsum$ID,
                           datAsum$Date), ]

## keys + ordering
datAsum$Date <- as.Date(datAsum$Date)
datAsum <- datAsum[order(datAsum$Garden_ID, datAsum$ID, datAsum$Date), ]
key <- paste(datAsum$Garden_ID, datAsum$ID, sep="|")
uid <- unique(key)

## define PeakDate = first date reaching N0 (survival starts here)
datAsum$PeakDate <- as.Date(NA)
keep <- rep(TRUE, nrow(datAsum))

for (u in uid) {
  i <- key == u
  r <- which(i & datAsum$Numind >= datAsum$N0)   # first reach of N0
  if (!length(r)) { keep[i] <- FALSE; next }       # never reaches N0 -> drop cohort
  d <- datAsum$Date[r[1]]
  datAsum$PeakDate[i] <- d
  keep[i] <- datAsum$Date[i] >= d                 # drop germination phase
}

datAsum <- datAsum[keep, ]
key <- paste(datAsum$Garden_ID, datAsum$ID, sep="|")  # refresh key

## no new seedlings after peak: remove Numind > N0 AFTER peak
bad <- datAsum$Numind > datAsum$N0
sum(bad)
datAsum <- datAsum[!bad, ]

## define survival time (AgeYears = 0 at peak)
datAsum$AgeDays <- as.numeric(difftime(datAsum$Date, datAsum$PeakDate, units="days"))

## 30 days window for the peak dat
medPeak <- median(datAsum$PeakDate, na.rm = TRUE)
shift_days <- as.numeric(datAsum$PeakDate - medPeak)
sd(shift_days)


## final binomial validity check (must be TRUE)
stopifnot(all(datAsum$Numind >= 0),
          all(datAsum$N0 > 0),
          all(datAsum$Numind <= datAsum$N0))

## ##################
##     model
## ##################

datAsum$Obs <- factor(seq_len(nrow(datAsum)))
datAsum$PeakDOY <- as.numeric(format(datAsum$PeakDate, "%j"))  # day-of-year, dropped because not significant

## create time grid
## Choose grid points for plotting
age_seq <- seq(0, as.numeric(quantile(datAsum$AgeDays, 0.95)), length.out = 200)
age_tsf <- log(age_seq + 1)
age2_tsf <- age_tsf^2
t_levels <- paste0("t", seq_along(age_seq))
## Map each observed AgeDaysTsf to nearest grid value
nearest_idx <- function(x, grid) {
  vapply(x, function(z) which.min(abs(grid - z)), integer(1))
}

datAsum$tgrid <- factor(t_levels[nearest_idx(datAsum$AgeDays, age_tsf)],
                        levels = t_levels)

seasons1 <- c(-30, 0, 30, 90, 180, 270)
seasons2 <- c(-30, 0, 60, 90, 180)
datAsum$AgeBin <- cut(datAsum$AgeDays,
  breaks = c(seasons1, 365 + seasons2, Inf),
  include.lowest = TRUE
  )
levels(datAsum$AgeBin)
tab <- table(datAsum$AgeBin)
tab
tab[tab < 20]   # bins with little data


## run the model
## #########################
m_surv <- asreml(
  fixed  = Numind ~ AgeBin,
  random = ~ Garden_ID * ID + Obs,
  ##  family = asr_binomial(total = N0, dispersion=NA), all bounded
  family = asr_binomial(total = N0),
  data   = datAsum,
  maxit  = 50
)
m_surv <- update.asreml(m_surv)
summary(m_surv)
wald.asreml(m_surv)

## ######################
## prediction
## ######################
grid_g <- expand.grid(
  AgeBin    = levels(datAsum$AgeBin),
  Garden_ID = levels(datAsum$Garden_ID)
)

pred_g <- predict(m_surv,
  classify = "AgeBin:Garden_ID",
  newdata  = grid_g,
  type     = "link"
)

pred_g_df <- data.frame(
  AgeBin    = pred_g$pvals$AgeBin,
  Garden_ID = pred_g$pvals$Garden_ID,
  eta       = pred_g$pvals$predicted.value
)

ilogit <- function(x) exp(x)/(1+exp(x))
pred_g_df$fit <- ilogit(pred_g_df$eta)

## population mean and CI
## #############################

grid_pop <- data.frame(AgeBin = levels(datAsum$AgeBin))

pred_pop <- predict(m_surv,
  classify="AgeBin",
  newdata=grid_pop,
  type="link",
  sed=TRUE
)

pred_pop_df <- na.omit(data.frame(
  AgeBin = pred_pop$pvals$AgeBin,
  eta    = pred_pop$pvals$predicted.value,
  se     = pred_pop$pvals$std.error
))

pred_pop_df$fit <- ilogit(pred_pop_df$eta)
pred_pop_df$lo  <- ilogit(pred_pop_df$eta - 1.96*pred_pop_df$se)
pred_pop_df$hi  <- ilogit(pred_pop_df$eta + 1.96*pred_pop_df$se)

## mid-points for bins
bin_mid <- function(lbl) {
  # lbl like "[-30,0]" or "(0,30]" or "(545,Inf]"
  a <- sub("^[\\[(]([^,]+),.*$", "\\1", lbl)
  b <- sub("^.*,(.+)[\\])]$", "\\1", lbl)
  a <- as.numeric(a)
  b <- suppressWarnings(as.numeric(b))
  if (is.na(b)) b <- a
  (a + b) / 2
}

pred_g_df$AgeMid <- sapply(as.character(pred_g_df$AgeBin), bin_mid)
pred_pop_df$AgeMid <- sapply(as.character(pred_pop_df$AgeBin), bin_mid)


## ######################################################
##                     Fagus
## ######################################################

datF <- subset(dat.all, Genus=="Fagus")

## make factors
datF$Garden_ID <- as.factor(datF$Garden_ID)
datF$Block <- as.factor(datF$Block)
datF$Garden_block_spot <- as.factor(datF$Garden_block_spot)
datF$ID <- as.factor(datF$ID)
datF$Year <- format(datF$Date, "%Y")

## order by Garden_block_spot, then Date
datF <- datF[order(datF$Garden_ID, datF$Garden_block_spot, datF$Date), ]

## change 2025 entries to 2024 last date
last.date.2024 <- max(datF[datF$Year==2024, "Date"])
datF$Date[datF$Year==2025] <- last.date.2024
datF <- subset(datF, Year<2025)

## cohort aggregate
datFsum <- aggregate(Numind ~ Garden_ID + ID + Date,
                      data = datF, FUN = sum)
datFsum$Year <- format(datFsum$Date, "%Y")

## Define N0 = max aggregated count in 2022 (max germination rate)
datFsum$N0 <- NA_real_
uid <- unique(paste(datFsum$Garden_ID, datFsum$ID, sep="|"))
for (u in uid) {
  i <- paste(datFsum$Garden_ID, datFsum$ID, sep="|") == u
  j <- i & datFsum$Year == "2022"
  if (any(j)) datFsum$N0[i] <- max(datFsum$Numind[j], na.rm=TRUE)
}
datFsum <- datFsum[datFsum$N0 > 0, ] ## drop cohorts never germinated

## order
datFsum <- datFsum[order(datFsum$Garden_ID,
                           datFsum$ID,
                           datFsum$Date), ]

## keys + ordering
datFsum$Date <- as.Date(datFsum$Date)
datFsum <- datFsum[order(datFsum$Garden_ID, datFsum$ID, datFsum$Date), ]
key <- paste(datFsum$Garden_ID, datFsum$ID, sep="|")
uid <- unique(key)

## define PeakDate = first date reaching N0 (survival starts here)
datFsum$PeakDate <- as.Date(NA)
keep <- rep(TRUE, nrow(datFsum))

for (u in uid) {
  i <- key == u
  r <- which(i & datFsum$Numind >= datFsum$N0)   # first reach of N0
  if (!length(r)) { keep[i] <- FALSE; next }       # never reaches N0 -> drop cohort
  d <- datFsum$Date[r[1]]
  datFsum$PeakDate[i] <- d
  keep[i] <- datFsum$Date[i] >= d                 # drop germination phase
}

datFsum <- datFsum[keep, ]
key <- paste(datFsum$Garden_ID, datFsum$ID, sep="|")  # refresh key

## no new seedlings after peak: remove Numind > N0 AFTER peak
bad <- datFsum$Numind > datFsum$N0
sum(bad)
datFsum <- datFsum[!bad, ]

## define survival time (AgeYears = 0 at peak)
datFsum$AgeDays <- as.numeric(difftime(datFsum$Date, datFsum$PeakDate, units="days"))

## 30 days window for the peak dat
medPeak <- median(datFsum$PeakDate, na.rm = TRUE)
shift_days <- as.numeric(datFsum$PeakDate - medPeak)
sd(shift_days)


## final binomial validity check (must be TRUE)
stopifnot(all(datFsum$Numind >= 0),
          all(datFsum$N0 > 0),
          all(datFsum$Numind <= datFsum$N0))

## ##################
##     model
## ##################

datFsum$Obs <- factor(seq_len(nrow(datFsum)))
datFsum$PeakDOY <- as.numeric(format(datFsum$PeakDate, "%j"))  # day-of-year, dropped because not significant

## create time grid
## Choose grid points for plotting
age_seq <- seq(0, as.numeric(quantile(datFsum$AgeDays, 0.95)), length.out = 200)
age_tsf <- log(age_seq + 1)
age2_tsf <- age_tsf^2
t_levels <- paste0("t", seq_along(age_seq))
## Map each observed AgeDaysTsf to nearest grid value
nearest_idx <- function(x, grid) {
  vapply(x, function(z) which.min(abs(grid - z)), integer(1))
}

datFsum$tgrid <- factor(t_levels[nearest_idx(datFsum$AgeDays, age_tsf)],
                        levels = t_levels)

seasons1 <- c(-30, 0, 30, 90, 180, 270)
seasons2 <- c(-30, 0, 60, 90, 180)
datFsum$AgeBin <- cut(datFsum$AgeDays,
  breaks = c(seasons1, 365 + seasons2, Inf),
  include.lowest = TRUE
  )
levels(datFsum$AgeBin)
tab <- table(datFsum$AgeBin)
tab
tab[tab < 20]   # bins with little data


## run the model
## #########################
m_survF <- asreml(
  fixed  = Numind ~ AgeBin,
  random = ~ Garden_ID * ID + Obs,
  family = asr_binomial(total = N0),
  data   = datFsum,
  maxit  = 50
)
m_survF <- update.asreml(m_survF)
summary(m_survF)
wald.asreml(m_survF)

## version with all components
m_survF_dp <- asreml(
  fixed  = Numind ~ AgeBin,
  random = ~ Garden_ID * ID + Obs,
  family = asr_binomial(total = N0, dispersion=NA),
  data   = datFsum,
  maxit  = 50
)
## not stable....

## ######################
## prediction
## ######################
grid_g_F <- expand.grid(
  AgeBin    = levels(datFsum$AgeBin),
  Garden_ID = levels(datFsum$Garden_ID)
)

pred_g_F <- predict(m_survF,
  classify = "AgeBin:Garden_ID",
  newdata  = grid_g_F,
  type     = "link"
)

pred_g_F_df <- data.frame(
  AgeBin    = pred_g_F$pvals$AgeBin,
  Garden_ID = pred_g_F$pvals$Garden_ID,
  eta       = pred_g_F$pvals$predicted.value
)

ilogit <- function(x) exp(x)/(1+exp(x))
pred_g_F_df$fit <- ilogit(pred_g_F_df$eta)

## population mean and CI
## #############################

grid_pop <- data.frame(AgeBin = levels(datFsum$AgeBin))

pred_pop_F <- predict(m_survF,
  classify="AgeBin",
  newdata=grid_pop,
  type="link",
  sed=TRUE
)

pred_pop_F_df <- na.omit(data.frame(
  AgeBin = pred_pop_F$pvals$AgeBin,
  eta    = pred_pop_F$pvals$predicted.value,
  se     = pred_pop_F$pvals$std.error
))

pred_pop_F_df$fit <- ilogit(pred_pop_F_df$eta)
pred_pop_F_df$lo  <- ilogit(pred_pop_F_df$eta - 1.96*pred_pop_F_df$se)
pred_pop_F_df$hi  <- ilogit(pred_pop_F_df$eta + 1.96*pred_pop_F_df$se)

pred_g_F_df$AgeMid <- sapply(as.character(pred_g_F_df$AgeBin), bin_mid)
pred_pop_F_df$AgeMid <- sapply(as.character(pred_pop_F_df$AgeBin), bin_mid)


## ============================
## Combine Abies + Fagus preds
## ============================

## 0) Tag species
pred_g_df$Species     <- "Abies"
pred_pop_df$Species   <- "Abies"
pred_g_F_df$Species   <- "Fagus"
pred_pop_F_df$Species <- "Fagus"

## 1) Make sure AgeBin is character before rbind (avoids factor level mismatches)
pred_g_df$AgeBin      <- as.character(pred_g_df$AgeBin)
pred_pop_df$AgeBin    <- as.character(pred_pop_df$AgeBin)
pred_g_F_df$AgeBin    <- as.character(pred_g_F_df$AgeBin)
pred_pop_F_df$AgeBin  <- as.character(pred_pop_F_df$AgeBin)

## 2) Row-bind
pred_g_all   <- rbind(pred_g_df,   pred_g_F_df)
pred_pop_all <- rbind(pred_pop_df, pred_pop_F_df)

## 3) Enforce common AgeBin ordering (use Abies bin levels)
age_levels <- levels(datAsum$AgeBin)

pred_g_all$AgeBin   <- factor(pred_g_all$AgeBin,   levels = age_levels)
pred_pop_all$AgeBin <- factor(pred_pop_all$AgeBin, levels = age_levels)

## Drop any rows whose bins don't exist in the chosen ordering
pred_g_all   <- pred_g_all[!is.na(pred_g_all$AgeBin), ]
pred_pop_all <- pred_pop_all[!is.na(pred_pop_all$AgeBin), ]

age_mid <- sapply(age_levels, bin_mid)
mid_map <- data.frame(
  AgeBin  = factor(age_levels, levels = age_levels),
  AgeMid  = as.numeric(age_mid),
  stringsAsFactors = FALSE
)

## 5) Add AgeMid to both prediction tables
pred_g_all  <- merge(pred_g_all,  mid_map, by = "AgeBin", all.x = TRUE, sort = FALSE)
pred_pop_all <- merge(pred_pop_all, mid_map, by = "AgeBin", all.x = TRUE, sort = FALSE)

pred_g_all$AgeMid <- pred_g_all$AgeMid.x
pred_pop_all$AgeMid <- pred_pop_all$AgeMid.x

## Ensure proper ordering for plotting
pred_g_all  <- pred_g_all[order(pred_g_all$Species, pred_g_all$Garden_ID, pred_g_all$AgeMid), ]
pred_pop_all <- pred_pop_all[order(pred_pop_all$Species, pred_pop_all$AgeMid), ]

## ============================
## Season shading + labels
## ============================

## season colors
seasons_df <- data.frame(
  xmin = c(-30,   0,   90, 180, 270, 335, 365, 455),
  xmax = c(  0,  90,  180, 270, 335, 365, 455, 545),
  season = c(
    "spring",  ## pre-peak / emergence
    "summer",
    "autumn",
    "winter",
    "spring",
    "spring",  ## late spring transition
    "summer",
    "autumn"))

## strong but not overwhelming greys
season_greys <- c(
  spring = "grey90",
  summer = "grey80",
  autumn = "grey70",
  winter = "grey60"
)

species_cols <- c(
  "Abies" = "#E18600",
  "Fagus" = "#5B84B1"
)

## Build label positions (top of panel)
season_labels <- data.frame(
  xmid  = (seasons_df$xmin + seasons_df$xmax)/2,
  label = tools::toTitleCase(as.character(seasons_df$season)),
  stringsAsFactors = FALSE
)

season_labels$label[6] <- ""


## ============================
## Garden labels at end (>= 0.2)
## ============================

thr <- 0.15

## End point per Species x Garden_ID = last AgeMid row
key <- paste(pred_g_all$Species, pred_g_all$Garden_ID, sep="|")
end_df <- do.call(rbind, lapply(split(pred_g_all, key), function(d) {
  d[which.max(d$AgeMid), c("Species","Garden_ID","AgeMid","fit")]
}))
rownames(end_df) <- NULL

lab_df <- end_df[end_df$fit >= thr, ]

## Simple vertical spreading to reduce overlap
spread_y <- function(y, min_sep = 0.03, ymin = 0.02, ymax = 0.98) {
  y <- pmin(pmax(y, ymin), ymax)
  o <- order(y)
  y2 <- y[o]
  if (length(y2) > 1) {
    for (i in 2:length(y2)) if ((y2[i] - y2[i-1]) < min_sep) y2[i] <- y2[i-1] + min_sep
    for (i in (length(y2)-1):1) if ((y2[i+1] - y2[i]) < min_sep) y2[i] <- y2[i+1] - min_sep
  }
  y2 <- pmin(pmax(y2, ymin), ymax)
  out <- y; out[o] <- y2
  out
}
lab_df$ylab <- spread_y(lab_df$fit, min_sep = 0.035)

## Margin x positions (no x-axis extension; we draw outside panel)
x_end <- max(pred_g_all$AgeMid, na.rm = TRUE)
x_min <- min(pred_g_all$AgeMid, na.rm = TRUE)
x_rng <- x_end - x_min

lab_df$x0   <- lab_df$AgeMid
lab_df$x1   <- x_end + 0.015 * x_rng
lab_df$xtxt <- x_end + 0.020 * x_rng

## ============================
## Plot + PDF export
## ============================

species_cols <- c("Abies"="#E18600","Fagus"="#5B84B1")

## x-axis ticks at bin mids, labels as days
x_breaks <- mid_map$AgeMid
x_labels <- as.character(round(mid_map$AgeMid))

##pdf("Figure_survival.pdf", width = 10.5, height = 7.5, onefile = TRUE)

p <- ggplot() +

  ## seasonal background (grey shades; no legend)
  geom_rect(
    data = seasons_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = season_greys[as.character(seasons_df$season)],
    alpha = 0.35,
    inherit.aes = FALSE
  ) +

  ## season labels above the panel
  geom_text(
    data = season_labels,
    aes(x = xmid, y = 1.04, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.5,
    colour = "grey25"
  ) +

  ## garden-specific curves (species color; NO legend)
  geom_line(
    data = pred_g_all,
    aes(x = AgeMid, y = fit, group = interaction(Species, Garden_ID), colour = Species),
    alpha = 0.30, linewidth = 0.9,
    show.legend = FALSE
  ) +

  ## CI ribbons (species fill)
  geom_ribbon(
    data = pred_pop_all,
    aes(x = AgeMid, ymin = lo, ymax = hi, fill = Species),
    alpha = 0.22
  ) +

  ## mean curves (species colour)
  geom_line(
    data = pred_pop_all,
    aes(x = AgeMid, y = fit, colour = Species),
    linewidth = 2.1
  ) +

  ## tick leader lines (no arrows; NO legend)
  geom_segment(
    data = lab_df,
    aes(x = x0, y = fit, xend = x1, yend = ylab, colour = Species),
    inherit.aes = FALSE,
    linewidth = 0.6,
    show.legend = FALSE
  ) +

  ## garden ID labels in right margin (NO legend)
  geom_text(
    data = lab_df,
    aes(x = xtxt, y = ylab, label = Garden_ID, colour = Species),
    inherit.aes = FALSE,
    hjust = 0, size = 5.0, fontface = "bold",
    show.legend = FALSE
  ) +

  scale_colour_manual(values = species_cols) +
  scale_fill_manual(values = species_cols) +

  ##scale_x_continuous(breaks = x_breaks, labels = x_labels) +

  labs(
    x = "Days since peak cohort abundance",
    y = "Survival probability",
    colour = "Species",
    fill   = "Species"
  ) +

  coord_cartesian(ylim = c(0, 1), clip = "off") +

  theme_classic(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14),
    plot.margin  = margin(t = 22, r = 120, b = 10, l = 10)
  )  

print(p)
##dev.off()

## other panel with variance components
## ######################################

## helper: extract + clean varcomp table
vc_extract <- function(mod, species) {
  vc <- mod$varcomp
  vc$Term <- rownames(vc)
  vc$Species <- species

  ## keep only the entries we want to show
  keep <- vc$Term %in% c("Garden_ID", "Garden_ID:ID", "ID", "Obs", "units!R")
  vc <- vc[keep, c("Species","Term","component","std.error","z.ratio")]

  ## numeric safety
  vc$component <- as.numeric(vc$component)
  vc$std.error <- as.numeric(vc$std.error)
  vc$z.ratio   <- as.numeric(vc$z.ratio)

  ## rename to nice labels
  lab <- c(
      "units!R"      = "Residual",
      "Obs"          = "OLRE",
      "ID"           = "Provenance",
      "Garden_ID:ID" = "Garden × Provenance",
    "Garden_ID"    = "Garden")
  vc$Component <- lab[vc$Term]

  ## drop NAs and exact zeros if you prefer (optional)
  vc <- vc[is.finite(vc$component), ]

  ## compute proportions (latent scale)
  tot <- sum(vc$component, na.rm = TRUE)
  vc$PropVar <- vc$component / tot
  vc$PropSE  <- vc$std.error / tot  # simple approx (good enough for plotting)

  ## order for plotting (top to bottom)
  vc$Component <- factor(vc$Component,
                         levels = c("Residual","OLRE","Provenance","Garden × Provenance","Garden"))

  vc
}

vc_abies <- vc_extract(summary(m_surv), "Abies")
vc_abies <- vc_abies[nrow(vc_abies):1, ]
vc_fagus <- vc_extract(summary(m_survF), "Fagus")
vc_fagus <- vc_fagus[nrow(vc_fagus):1, ]

var_df <- rbind(vc_abies, vc_fagus)

### combined plot
## #####################

p_var <- ggplot(var_df, aes(x = PropVar, y = Component, fill = Species)) +

  geom_col(position = position_dodge(width = 0.7),
           width = 0.65) +

  geom_errorbarh(aes(xmin = pmax(0, PropVar - 1.96*PropSE),
                     xmax = PropVar + 1.96*PropSE),
                 position = position_dodge(width = 0.7),
                 height = 0.25, linewidth = 0.6) +

  scale_fill_manual(values = species_cols) +

  labs(x = "Proportion of latent-scale variance", y = NULL, fill = "Species") +

  theme_classic(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14),
    plot.margin  = margin(t = 22, r = 10, b = 10, l = 10)
  )


final_fig <- p + p_var +
    plot_layout(widths = c(1.7, 1)) +
    plot_annotation(
        tag_levels = "a",                     # lowercase labels: a, b, c
        tag_prefix = "",                      # no parentheses
        tag_suffix = "") &
    theme(
        legend.position = "top",
        legend.justification = "center",
        legend.direction = "horizontal",
        panel.spacing = unit(0, "pt"),
        plot.margin = margin(10, 10, 10, 10),
        plot.tag.position = c(0, 1)
  )

ggsave("survival_with_varcomp.pdf", final_fig,
       width = 15, height = 7, device = cairo_pdf)

ggsave("/home/kati/Dropbox/Overleaf/MyGardenOfTrees_PilotPhase_MS/New_Figures/Figure_survival.pdf", final_fig,
       width = 15, height = 7, device = cairo_pdf)

