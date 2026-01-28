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
library(glmnet)
library(readxl)

##    all data
## ###############
load("dat2022.RData") 
dat.all <- subset(dat.all, Date <= "2022-07-01") ## germination process only - no summer or autumn survival here

##    Abies
## ###############

datA <- subset(dat.all, Genus=="Abies")
## remove provenances that are not tested everywhere
##datA <- subset(datA, ! ID %in% c("AA_HR1", "AA_PL1", "AA_SI1", "AA_IT1"))
datA <- subset(datA, ! ID %in% c("AA_HR1","AA_IT1"))
## remove garden where sample sizes are small
##datA <- subset(datA, ! Garden_ID %in% c(3, 5, 16))

## make factors
datA$Garden_ID <- as.factor(datA$Garden_ID)
datA$Block <- as.factor(datA$Block)
datA$Garden_block_spot <- as.factor(datA$Garden_block_spot)
datA$ID <- as.factor(datA$ID)

##date
datA$Date <- as.Date(datA$Date)
datA$Date <- factor(datA$Date, ordered=T)
datA$Date_int <- as.integer(as.factor(datA$Date))

## response variable: cummulative max germination rate
datA$germ.max.all <- ave(datA$germ, datA$Garden_ID, datA$Block, datA$Seedling_spot_true, datA$ID, FUN = cummax, na.rm=T)

## order by Garden_block_spot, then Date
datA <- datA[order(datA$Garden_block_spot, datA$Date), ]

## GET RID OF TIME: keep last row of each Garden_block_spot
datA <- datA[!duplicated(datA$Garden_block_spot, fromLast = TRUE), ]

## order
datA <- datA[order(datA$ID, datA$Garden_ID, datA$Block, datA$Garden_block_spot), ]


##    Fagus
## ###############

datF <- subset(dat.all, Genus=="Fagus")
## remove provenances that are not tested everywhere
datF <- subset(datF, ! ID %in% c("FS_AT1"))
datF <- subset(datF, ! Garden_ID %in% c(3, 16))

## make factors
datF$Garden_ID <- as.factor(datF$Garden_ID)
datF$Block <- as.factor(datF$Block)
datF$Garden_block_spot <- as.factor(datF$Garden_block_spot)
datF$ID <- as.factor(datF$ID)

##date
datF$Date <- as.Date(datF$Date)
datF$Date <- factor(datF$Date, ordered=T)
datF$Date_int <- as.integer(as.factor(datF$Date))

## response variable: cummulative max germination rate
datF$germ.max.all <- ave(datF$germ, datF$Garden_ID, datF$Block, datF$Seedling_spot_true, datF$ID, FUN = cummax, na.rm=T)

## order by Garden_block_spot, then Date
datF <- datF[order(datF$Garden_block_spot, datF$Date), ]

## GET RID OF TIME: keep last row of each Garden_block_spot
datF <- datF[!duplicated(datF$Garden_block_spot, fromLast = TRUE), ]

## order
datF <- datF[order(datF$ID, datF$Garden_ID, datF$Block, datF$Garden_block_spot), ]


save(datA, datF, file="dat2022_germ.RData")

## ############################################################
## lasso regression for the environmental variables - ABIES
## ##########################################################

env_vars <- read_xlsx(file="Table_S1_environmental_variables_LASSO.xlsx")
X <- datA[, match(env_vars$Variable_short_name, names(datA))]
X <- X[, apply(X, 2, function(a) length(na.omit(a))==length(a))] ## select columns without NA
y <- datA$germ.max.all
X_scaled <- scale(X)
X_mat <- model.matrix(~ . - 1, data = as.data.frame(X_scaled))

set.seed(123)
lasso_fit <- cv.glmnet(
  X_mat, y,
  family = "poisson",
  alpha = 1          
)

## Extract coefficients
coef_df <- as.data.frame(as.matrix(coef(lasso_fit, s = "lambda.1se")))
coef_df$variable <- rownames(coef_df)
colnames(coef_df)[1] <- "coef"

## join with type
coef_df <- merge(coef_df, env_vars)

## Filter out intercept and zero coefficients
coef_df <- coef_df[coef_df$variable != "(Intercept)" & coef_df$coef != 0, ]

## Order by absolute coefficient
coef_df <- coef_df[order(abs(coef_df$coef)), ]

## table
lasso_table <- data.frame(
  variable = coef_df$variable,
  coef     = coef_df$coef,
  abs_coef = abs(coef_df$coef),
  type     = coef_df$type,
  sign     = ifelse(coef_df$coef > 0, "+", "-"),
  row.names = NULL,
  stringsAsFactors = FALSE
)
## sort by absolute strength
lasso_table <- lasso_table[order(-lasso_table$abs_coef), ]
write.csv(lasso_table, "LASSO_coefficients_abies.csv", row.names = FALSE)

## reducing variables to a smaller set:
p <- 1   # top %
n_keep <- ceiling(nrow(lasso_table) * p)
selected <- head(lasso_table[order(-lasso_table$abs_coef), ], n_keep)
X_sub <- X[, selected$variable, drop = FALSE]
cors <- cor(X_sub, use = "pairwise.complete.obs")
d <- as.dist(1 - abs(cors))
hc <- hclust(d, method = "average")
plot(hc, cex = 0.7, main = "Clustering of LASSO-selected variables")
clusters <- cutree(hc, h = 0.7)  ## adjust h to control cluster size
group_list_abies <- split(names(clusters), clusters)
lasso_table_abies <- lasso_table

## Plot
## #############
oi <- c(
  "Provenance climate" = "#0072B2",  # blue
  "Provenance soil"    = "#56B4E9",  # lighter blue

  "Garden climate"     = "#009E73",  # green
  "Garden soil"        = "#5DC8A8"   # lighter green
)

pdf("LASSO_coefficients_plot_abies.pdf", width = 7, height = 8) 

par(
  mar = c(6, 12.5, 4, 12),
  xpd = NA,
  family = "sans",
  cex = 0.9
)

cols <- oi  

bar_cols <- cols[lasso_table$type]

ord <- order(lasso_table$coef)
lasso_table <- lasso_table[ord, ]
bar_cols <- bar_cols[ord]

bp <- barplot(
  lasso_table$coef,
  horiz = TRUE,
  col = bar_cols,
  border = NA,
  las = 1,
  xlab = "LASSO coefficient",
  main = "LASSO-selected Predictors",
  space = 0.7
)

text(
  x = min(lasso_table$coef) - 0.12 * diff(range(lasso_table$coef)),
  y = bp,
  labels = lasso_table$variable,
  pos = 2,
  cex = 0.85
)

abline(v = 0, lwd = 1.4, lty = 2, col = "grey40")
abline(v = pretty(lasso_table$coef), col = "#D0D0D050", lwd = 0.7)

legend(
  x = max(lasso_table$coef) * 1.3,
  y = max(bp),
  legend = names(cols),
  fill = cols,
  border = NA,
  cex = 0.9,
  bty = "n"
)

dev.off()

## ############################################################
## lasso regression for the environmental variables - FAGUS
## ########################################################## 

env_vars <- read.csv(file="env_vars_modified.csv")
X <- datF[, match(env_vars$Variable_short_name, names(datF))]
X <- X[, apply(X, 2, function(a) length(na.omit(a))==length(a))]
y <- datF$germ.max.all
X_scaled <- scale(X)
X_mat <- model.matrix(~ . - 1, data = as.data.frame(X_scaled))

set.seed(123)
lasso_fit <- cv.glmnet(
  X_mat, y,
  family = "poisson",
  alpha = 1          
)

## Extract coefficients
coef_df <- as.data.frame(as.matrix(coef(lasso_fit, s = "lambda.1se")))
coef_df$variable <- rownames(coef_df)
colnames(coef_df)[1] <- "coef"

## join with type
coef_df <- merge(coef_df, env_vars)

## Filter out intercept and zero coefficients
coef_df <- coef_df[coef_df$variable != "(Intercept)" & coef_df$coef != 0, ]

## Order by absolute coefficient
coef_df <- coef_df[order(abs(coef_df$coef)), ]

## table
lasso_table <- data.frame(
  variable = coef_df$variable,
  coef     = coef_df$coef,
  abs_coef = abs(coef_df$coef),
  type     = coef_df$type,
  sign     = ifelse(coef_df$coef > 0, "+", "-"),
  row.names = NULL,
  stringsAsFactors = FALSE
)
## sort by absolute strength
lasso_table <- lasso_table[order(-lasso_table$abs_coef), ]
write.csv(lasso_table, "LASSO_coefficients_fagus.csv", row.names = FALSE)

## reducing variables to a smaller set:
p <- 0.50   # top %
n_keep <- ceiling(nrow(lasso_table) * p)
selected <- head(lasso_table[order(-lasso_table$abs_coef), ], n_keep)
X_sub <- X[, selected$variable, drop = FALSE]
cors <- cor(X_sub, use = "pairwise.complete.obs")
d <- as.dist(1 - abs(cors))
hc <- hclust(d, method = "average")
plot(hc, cex = 0.7, main = "Clustering of LASSO-selected variables")
clusters <- cutree(hc, h = 0.6)  # adjust h to control cluster size
group_list_fagus <- split(names(clusters), clusters)

lasso_table_fagus <- lasso_table

## Plot
## #############

pdf("LASSO_coefficients_plot_fagus.pdf", width = 7, height = 8) 

par(
  mar = c(6, 12.5, 4, 12),
  xpd = NA,
  family = "sans",
  cex = 0.9
)

cols <- oi 

bar_cols <- cols[lasso_table$type]

ord <- order(lasso_table$coef)
lasso_table <- lasso_table[ord, ]
bar_cols <- bar_cols[ord]

bp <- barplot(
  lasso_table$coef,
  horiz = TRUE,
  col = bar_cols,
  border = NA,
  las = 1,
  xlab = "LASSO coefficient",
  main = "LASSO-selected Predictors",
  space = 0.7
)

text(
  x = min(lasso_table$coef) - 0.12 * diff(range(lasso_table$coef)),
  y = bp,
  labels = lasso_table$variable,
  pos = 2,
  cex = 0.85
)

abline(v = 0, lwd = 1.4, lty = 2, col = "grey40")
abline(v = pretty(lasso_table$coef), col = "#D0D0D050", lwd = 0.7)

legend(
  x = max(lasso_table$coef) * 1.3,
  y = max(bp),
  legend = names(cols),
  fill = cols,
  border = NA,
  cex = 0.9,
  bty = "n"
)

dev.off()


save(group_list_abies, group_list_fagus, lasso_table_abies,lasso_table_fagus, file="best_lasso_variables.RData")

intersect(unlist(group_list_fagus), unlist(group_list_abies))

