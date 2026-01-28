library(asreml)
##asreml.license.activate()
## HCEE-JCDA-JBFD-BJFA
library(lubridate)
library(RColorBrewer)

library(ggplot2)
library(patchwork)

## ======================================================================

## understanding which environmental factors that affect germination
## and the relative roles of garden vs env (variance components)

## ======================================================================

setwd('/home/kati/Dropbox/projects/MyGardenOfTrees_trials/Manuscripts/MGOT_PilotTrials/analysis/gardens/')
load("dat2022_germ.RData") ## see 2_lasso_varselect.r
load(file="best_lasso_variables.RData")
##unlist(group_list_fagus); unlist(group_list_abies)
intersect(unlist(group_list_fagus), unlist(group_list_abies))

datA$Prcp_2022_spring <- apply(datA[,c("Prcp_2022_3", "Prcp_2022_4", "Prcp_2022_5")], 1, sum)
datF$Prcp_2022_spring <- apply(datF[,c("Prcp_2022_3", "Prcp_2022_4", "Prcp_2022_5")], 1, sum)
datF$Prcp_seed_prov <- apply(datF[,c("Prcp_2022_3", "Prcp_2022_4", "Prcp_2022_5")], 1, sum)
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

## fixed effects
myfixed <- as.formula(paste("germ.max.all ~", paste(env_vars, collapse = " + ")))
## two missing provenances for each species!
myfixed_seed <- as.formula(paste("germ.max.all ~", paste(c(seed_vars, env_vars), collapse = " + ")))

datA$GardenBlock <- as.factor(datA$GardenBlock)
datF$GardenBlock <- as.factor(datF$GardenBlock)

## mixed models Abies 
## ##################
modA <- asreml(
  fixed  = myfixed_seed,
  random = ~ ID * (Garden_ID + GardenBlock),
  family = asr_negative.binomial(dispersion=NA),
  na.action = na.method(x = ""),
  maxiter = 80,
  ai.sing = TRUE,
  data   = datA
)
modA <- update.asreml(modA)
summary(modA)$varcomp
wald.asreml(modA)
pdf("asreml_fit_abies.pdf")
plot(modA)
dev.off()
## variance components
barplot(summary(modA)$varcomp[,1],
        names.arg=unlist(lapply(strsplit(rownames(summary(modA)$varcomp),"!"), function(a) a[1])))
## fixed effects
barplot(modA$coefficients$fixed[,1])

## mixed models Fagus
## ##################
modF <- asreml(
  fixed  = myfixed_seed,
  random = ~ ID * (Garden_ID + GardenBlock),
  family = asr_negative.binomial(dispersion=NA),
  na.action = na.method(x = ""),
  maxiter = 80,
  ai.sing = TRUE,
  data   = datF
)
modF <- update.asreml(modF)
summary(modF)$varcomp
wald.asreml(modF)
pdf("asreml_fit_fagus.pdf")
plot(modF)
dev.off()

### figure
## ####################

### ===================================================================
### 1. CLEAN LOOKUP TABLE
### ===================================================================

lookup <- data.frame(
    var = c(seed_vars, env_vars),
  nice = c("Seed weight", "Seed moisture content", "Garden precipitation spring", "Garden temperature spring", "Garden soil clay content", "Garden soil silt content", "Garden soil nitrogen content",  "Provenance total annual precipitation"),
  stringsAsFactors = FALSE
)


### ===================================================================
### 2. RANDOM EFFECTS (combined panel)
### ===================================================================

vcA <- as.data.frame(summary(modA)$varcomp)
vcF <- as.data.frame(summary(modF)$varcomp)

vcA$Effect <- rownames(vcA)
vcF$Effect <- rownames(vcF)

vcA$Species <- "Abies"
vcF$Species <- "Fagus"

vc_df <- rbind(vcA, vcF)
rownames(vc_df) <- NULL

names(vc_df)[1:3] <- c("component","std_error","z")

# Rename effects
vc_df$Effect_clean <- vc_df$Effect
vc_df$Effect_clean[vc_df$Effect=="ID"] <- "Provenance"
vc_df$Effect_clean[vc_df$Effect=="Garden_ID"] <- "Garden"
vc_df$Effect_clean[vc_df$Effect=="ID:Garden_ID"] <- "Garden × Provenance"
vc_df$Effect_clean[vc_df$Effect=="ID:GardenBlock"]  <- "Block × Provenance"
vc_df$Effect_clean[vc_df$Effect=="GardenBlock"] <- "Block"
vc_df$Effect_clean[vc_df$Effect=="units!R"] <- "Residual"

# Order
effect_order <- c("Block","Garden","Block × Provenance","Garden × Provenance","Provenance","Residual")
vc_df$Effect_clean <- factor(vc_df$Effect_clean, levels = effect_order)

# Proportion variance
# compute per species
Prop <- numeric(nrow(vc_df))
for(i in 1:nrow(vc_df)){
  sp <- vc_df$Species[i]
  tot <- sum(vc_df$component[vc_df$Species==sp])
  Prop[i] <- vc_df$component[i] / tot
}
vc_df$Prop <- Prop

# transparency for z < 2
vc_df$alpha <- 1##ifelse(vc_df$z < 2, 0.3, 1)

# Standard error of proportions (same scaling as component)
se_prop <- numeric(nrow(vc_df))
for(i in 1:nrow(vc_df)) {
  sp <- vc_df$Species[i]
  tot <- sum(vc_df$component[vc_df$Species == sp])
  se_prop[i] <- vc_df$std_error[i] / tot
}
vc_df$se_prop <- se_prop



### ===================================================================
### 3. FIXED EFFECTS (Abies)
### ===================================================================

fxA <- as.data.frame(modA$coefficients$fixed)
names(fxA) <- "Estimate"
fxA$Cov <- rownames(fxA)
waldA <- as.data.frame(wald.asreml(modA))[,3:4]
waldA <- head(waldA, nrow(fxA))
names(waldA) <- c("Wald","p")
fxA <- cbind(fxA, waldA)
fxA <- fxA[fxA$Cov != "(Intercept)", ]

## match lookup
fxA$Cov_clean <- fxA$Cov
for(i in 1:nrow(lookup)){
  fxA$Cov_clean[fxA$Cov == lookup$var[i]] <- lookup$nice[i]
}

# Force EXACT lookup order
fxA$Cov_clean <- factor(fxA$Cov_clean, levels = lookup$nice)

# significance stars
fxA$stars <- ""
fxA$stars[fxA$p < 0.05] <- "*"
fxA$stars[fxA$p < 0.01] <- "**"
fxA$stars[fxA$p < 0.001] <- "***"

fxA$alpha <- ifelse(fxA$p < 0.05, 1, 0.3)


### ===================================================================
### 4. FIXED EFFECTS (Fagus)
### ===================================================================

fxF <- as.data.frame(modF$coefficients$fixed)
names(fxF) <- "Estimate"
fxF$Cov <- rownames(fxF)
waldF <- as.data.frame(wald.asreml(modF))[,3:4]
waldF <- head(waldF, nrow(fxF))
names(waldF) <- c("Wald","p")
fxF <- cbind(fxF, waldF)
fxF <- fxF[fxF$Cov != "(Intercept)", ]

                                        # match lookup
# Assign clean names using lookup
fxF$Cov_clean <- fxF$Cov
for(i in 1:nrow(lookup)){
  fxF$Cov_clean[fxF$Cov == lookup$var[i]] <- lookup$nice[i]
}
fxF$Cov_clean <- factor(fxF$Cov_clean, levels = lookup$nice)

# stars
fxF$stars <- ""
fxF$stars[fxF$p < 0.05] <- "*"
fxF$stars[fxF$p < 0.01] <- "**"
fxF$stars[fxF$p < 0.001] <- "***"

fxF$alpha <- ifelse(fxF$p < 0.05, 1, 0.3)


### ===================================================================
### 5. COLORS
### ===================================================================

species_cols <- c("Abies"="#E18600","Fagus"="#5B84B1")


### ===================================================================
### 6. PLOTS
### ===================================================================

## RANDOM EFFECTS (one panel)
p_random <- ggplot(vc_df,
                   aes(x = Prop, y = Effect_clean,
                       fill = Species, alpha = alpha)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_errorbar(aes(xmin = Prop - se_prop,
                    xmax = Prop + se_prop),
                position = position_dodge(width = 0.7),
                width = 0.3) +
  scale_fill_manual(values = species_cols) +
  scale_alpha_identity(guide="none") +
  labs(title="",
       x="Proportion of variance", y="") +
    theme_classic(base_size = 13) +
  theme(
    legend.position="none",
    axis.text.y = element_text(size=11)
  )

## FIXED EFFECTS — Abies
p_fxA <- ggplot(fxA,
                aes(x=Estimate, y=Cov_clean,
                    fill="Abies", alpha=alpha)) +
  geom_col(width=0.6) +
  geom_text(aes(label=stars), hjust=-0.2, size=5) +
  scale_fill_manual(values=species_cols, guide="none") +
  scale_alpha_identity() +
  labs(title="",
       x="Effect size", y="") +
    theme_classic(base_size = 13) +
  theme(legend.position="none",
        axis.text.y=element_text(size=11),
        plot.margin = margin(5,0,5,5))


# FIXED EFFECTS — Fagus
p_fxF <- ggplot(fxF,
                aes(x=Estimate, y=Cov_clean,
                    fill="Fagus", alpha=alpha)) +
  geom_col(width=0.6) +
  geom_text(aes(label=stars), hjust=-0.2, size=5) +
  scale_fill_manual(values=species_cols, guide="none") +
  scale_alpha_identity() +
  labs(title="",
       x="Effect size", y=NULL) +
    theme_classic(base_size = 13) +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(5,5,5,0))


### ===================================================================
### 7. FINAL COMBINED FIGURE (ONE LEGEND)
### ===================================================================
final_plot <-
  (
    p_random /
    (p_fxA + p_fxF)
  ) +
  plot_layout(heights = c(1.5, 1), guides = "collect") +
  plot_annotation(
    tag_levels = "a",                     # lowercase labels: a, b, c
    tag_prefix = "",                      # no parentheses
    tag_suffix = "") &
  theme(
    legend.position = "top",
    panel.spacing   = unit(0, "pt"),
    plot.margin     = margin(10, 10, 10, 10),   # space around figure
    plot.tag.position = c(0, 1)              # x, y in npc units → OUTSIDE
  )

ggsave(
  filename = "Figure_mixed_model_mg.pdf",
  plot = final_plot,
  width = 8,      # adjust width (in inches)
  height = 9,     # adjust height
  units = "in",
  dpi = 300
)

ggsave(
  filename = "/home/kati/Dropbox/Overleaf/MyGardenOfTrees_PilotPhase_MS/New_Figures/Figure_mixed_model_mg.pdf",
  plot = final_plot,
  width = 8,      # adjust width (in inches)
  height = 10,     # adjust height
  units = "in",
  dpi = 300
)


