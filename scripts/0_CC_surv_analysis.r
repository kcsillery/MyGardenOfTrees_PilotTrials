## ##########################################################################

## Katalin Csillery, with input from Johannes Alt, Mert Celik,
## Marjorie Bison

## 28 Jan 2026

## For the manuscript:
## Gene–environment interactions govern early regeneration in fir and
## beech: evidence from participatory provenance trials across Europe
## by Katalin Csilléry, Justine Charlet de Sauvage, Madleina Caduff,
## Johannes Alt, Marjorie Bison, Mert Celik, Nicole Ponta, Daniel
## Wegmann

## ###########################################################################


library(Hmisc)
library("PerformanceAnalytics")
library(plyr)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(paletteer)
library(cowplot)
library(germinationmetrics)
library(forestmodel)

## raw data figures: time vs number of germinated seeds
## #####################################################
dat <- read.csv("climate_chamber_pheno_stages_seeds.csv")

## number of germinated seeds
dat$seed.sum <- apply(dat[, grep("Seed_", names(dat))], 1, function(a) sum(a>0))

## seed data
seed.dat <- read.csv("seed_provenance_details.csv")
dat <- merge(dat, seed.dat, by=c("Provenance", "Species"), all.x=T)
dat$Genus <- dat$Genus.x

dat.final <- subset(dat, Position_chamber=="finish")
dat.abs <- subset(dat, Position_chamber=="absolut")  ## check for roots that were not visible yet
dat <- subset(dat, ! Position_chamber %in% c("shelf"))

## removes "ROF1" "CHF1" with 0 final germinations
dat <- subset(dat, !ID %in% unique(dat.abs$ID[dat.abs$seed.sum==0]))

abies.plot =
    ggplot(NULL,
           aes(x = Days_since_setup, y = seed.sum)) + 
    geom_line(data=subset(dat, Species %in% c("Abies alba", "Abies nordmanniana") & Position_chamber!="absolut"),
              aes(color = ID, linetype=Moisture_treatment), linewidth = 1.2) +
    geom_jitter(data=subset(dat, Species %in% c("Abies alba", "Abies nordmanniana") & Position_chamber=="absolut"),
               aes(color = ID, shape=Moisture_treatment), size=3, position=position_jitter(width=.8, height=0)) +
    xlab("Time (Days)") + ylab("Cumulative germination count") +
    scale_x_continuous(breaks= seq(0,105, by=5)) +
    scale_y_continuous(breaks = seq (0,75,by = 5)) +
    scale_color_paletteer_d("ltc::crbhits") +
      ##ggtitle("Abies alba & Abies nordmanniana") +
    theme("light", text = element_text(size=rel(4)), legend.text = element_text(size=rel(4))) + theme_bw()
    

fagus.plot =
    ggplot(subset(dat, Species %in% c("Fagus sylvatica", "Fagus orientalis") & Position_chamber!="absolut"),
           aes(x = Days_since_setup, y = seed.sum)) + 
    geom_line(aes(color = ID, linetype=Moisture_treatment), linewidth = 1.2) +
    geom_jitter(data=subset(dat, Species %in% c("Fagus sylvatica", "Fagus orientalis") & Position_chamber=="absolut"),
               aes(color = ID, shape=Moisture_treatment), size=3, position=position_jitter(width=.8, height=0)) +
    xlab("Time (Days)") + ylab("Cumulative germination count") +
    scale_x_continuous(breaks= seq(0,105, by=5)) +
    scale_y_continuous(breaks = seq (0,80,by = 5)) +
   scale_color_paletteer_d("MoMAColors::Lupi") +
      ##ggtitle("Fagus sylvatica & Fagus caspica") +
    theme("light", text = element_text(size=rel(4)), legend.text = element_text(size=rel(4))) + theme_bw()

pdf("CC_rawdata.pdf", 10, 8)
plot_grid(abies.plot, fagus.plot, nrow = 2, align = "v", rel_heights = c(1/2, 1/2))
dev.off()

## Abies colors
abies.col = c("#CBC106FF", "#27993CFF", "#1C6838FF", "#8EBCB5FF", "#389CA7FF", "#4D83ABFF", "#CB7B26FF", "#BF565DFF", "#9E163CFF")
abies.id = levels(factor(dat$ID[dat$Genus=="Abies"]))

## Fagus colors
fagus.col = c("#61BEA4FF", "#B6E7E0FF", "#DAA5ACFF", "#98A54FFF") ## col for IT1: "#AA3F5DFF", not used below
fagus.col.all = c("#61BEA4FF", "#B6E7E0FF", "#AA3F5DFF", "#DAA5ACFF", "#98A54FFF") ## col for IT1: , not used below
fagus.id = levels(factor(dat$ID[dat$Genus=="Fagus"]))

## mycols = rbind(data.frame(Genus="Abies", ID=abies.id, col=abies.col),
##                data.frame(Genus="Fagus", ID=fagus.id, col=fagus.col.all))
## save(mycols, file="../mycols.RData")

## write.table(merge(mycols, seed.dat, by="ID", all=T), file="../mycols.txt")

## caluclate germination metrics
## ##############################

## exclude provenances with too little germination - ITF1
dat <- subset(dat, ! ID %in% c("FS_CH1", "FS_RO1", "FS_IT1"))

## seed.sum for the two treatments jointly
## number of germinated seeds for dry and moist together
seed.sum.cmb <- aggregate(dat$seed.sum, list(dat$ID, dat$Date), function(a) sum(a)/2)
names(seed.sum.cmb) <- c("ID", "Date", "seed.sum")
seed.sum.cmb$Moisture_treatment <- "joint"
dat.sub <- subset(dat, Moisture_treatment == "dry")
dat.sub$Moisture_treatment <- "joint"
dat.sub <- dat.sub[, -match("seed.sum", names(dat.sub))]
dat2 <- merge(dat.sub, seed.sum.cmb, by=c("ID", "Date", "Moisture_treatment"), all=T)
dat <- rbind(dat, dat2)
dat <- dat[, -grep("Seed_", names(dat))]

## germination statistics
germ.names <- c("GermPercent", "t50", "CVGermTime", "GermUncertainty", "GermSynchrony")
my.germ <- function(a, tot=100){
    int = a[, "Days_since_setup"]
    x = a[,"seed.sum"]-c(0, head(a[,"seed.sum"], -1))
    y = a[,"seed.sum"]
    out = c(round(GermPercent(germ.counts = x, total.seeds = tot)),
            round(t50(germ.counts = x, intervals = int, method = "coolbear")),
            round(CVGermTime(germ.counts = x, intervals = int),2),
            round(GermUncertainty(germ.counts = x, intervals = int),2),
            round(GermSynchrony(germ.counts = x, intervals = int),2))
    names(out) = germ.names
    out
}

long.prov = unique(dat$ID [dat$Days_since_setup > 73])
dat1 <- subset(dat, ID %in%long.prov  & Moisture_treatment %in% c("dry", "moist"))
res1 <- ddply(dat1, .(Genus, ID, Moisture_treatment), my.germ )
res1$cycle = "long"
dat2 <- subset(dat, Days_since_setup <73 & Moisture_treatment %in% c("dry", "moist"))
res2 <- ddply(dat2, .(Genus, ID, Moisture_treatment), my.germ)
res2$cycle = "short"
dat3 <- subset(dat, ID %in%long.prov & Moisture_treatment %in% c("joint"))
res3 <- ddply(dat3, .(Genus, ID, Moisture_treatment), my.germ, tot=100)
res3$cycle = "long"
dat4 <- subset(dat, Days_since_setup <73 & Moisture_treatment %in% c("joint"))
res4 <- ddply(dat4, .(Genus, ID, Moisture_treatment), my.germ, tot=100)
res4$cycle = "short"
res <- merge(res1, res2, all=T)
res <- merge(res, res3, all=T)
res <- merge(res, res4, all=T)

## combine with seed data
res <- merge(res, seed.dat, by=c("ID", "Genus"), all.x=T)
names(res)[16] = "Moisture_content"
seed.names = c("Latitude", "Longitude", "Altitude", "Seed_weight_g", "Moisture_content")
germ.vars = match(c(germ.names, seed.names), names(res))

## write out to a table
res.print =  res[, c("ID", "Genus", "Moisture_treatment", "cycle", "GermPercent", "t50", "GermUncertainty", "CVGermTime", "GermSynchrony")]
res.print = res.print[order(res$Genus, res$ID, res$Moisture_treatment, res$cycle), ]
##latex(res.print, rowlabel = NULL, rowname = NULL)

## combine with HMM data
load("mm_dat.RData")
all.res <- merge(res, mm.dat[c("ID", "g_MLE", "gamma_MLE", "delta_MLE", "alpha_MLE", "beta_MLE")], by="ID", all=T)
germ.vars = c(germ.names, "g_MLE", "delta_MLE", seed.names)
    
## plot correlations
pdf("CC_trait_correlations_abies.pdf", 10, 10)
chart.Correlation(subset(all.res, Genus == "Abies" & Moisture_treatment == "joint" & cycle == "short")[, germ.vars],
                  histogram=TRUE, pch=19, cex.lab=2)
dev.off()
pdf("CC_trait_correlations_fagus.pdf", 10, 10)
chart.Correlation(subset(all.res, Genus == "Fagus" & Moisture_treatment == "joint" & cycle == "short")[, germ.vars],
                  histogram=TRUE, pch=19)
dev.off()

pdf("delta_sw_lon.pdf",6,6)

par(mfrow=c(2,2), bty="l",mar=c(4,4,1,0))
myres = subset(all.res, Genus == "Abies" & Moisture_treatment == "joint" & cycle == "short")
plot(myres$g_MLE~myres$GermPercent, col=abies.col[as.numeric(factor(myres$ID))], cex=2, pch=20,
     ylab="Germination rate (g)", xlab="GermPercentage")
abline(lm(myres$g_MLE~myres$GermPercent), col="grey", lty=2)
mtext(paste("r =", round(cor.test(myres$g_MLE, myres$GermPercent)$estimate,2)), 3, line=-10, adj=1, cex=0.9)
mtext(paste("p-value =", round(cor.test(myres$g_MLE, myres$GermPercent)$p.value,3)), 3, line=-11, adj=1, cex=0.9)

myres = subset(all.res, Genus == "Fagus" & Moisture_treatment == "joint" & cycle == "short")
plot(myres$g_MLE~myres$GermPercent, col=fagus.col[as.numeric(factor(myres$ID))], cex=2, pch=20,
     ylab="", xlab="GermPercentage")
abline(lm(myres$g_MLE~myres$GermPercent), col="grey", lty=2)
mtext(paste("r =", round(cor.test(myres$g_MLE, myres$GermPercent)$estimate,2)), 3, line=-10, adj=1, cex=0.9)
mtext(paste("p-value =", round(cor.test(myres$g_MLE, myres$GermPercent)$p.value,3)), 3, line=-11, adj=1, cex=0.9)

myres = subset(all.res, Genus == "Abies" & Moisture_treatment == "joint" & cycle == "short")
plot(myres$delta_MLE~myres$Seed_weight_g, col=abies.col[as.numeric(factor(myres$ID))], cex=2, pch=20,
     ylab="Development speed (delta)", xlab="Seed weight")
abline(lm(myres$delta_MLE~myres$Seed_weight_g), col="grey", lty=2)
mtext(paste("r =", round(cor.test(myres$delta_MLE, myres$Seed_weight_g)$estimate,2)), 3, line=-1, adj=1, cex=0.9)
mtext(paste("p-value =", round(cor.test(myres$delta_MLE, myres$Seed_weight_g)$p.value,3)), 3, line=-2, adj=1, cex=0.9)

myres = subset(all.res, Genus == "Fagus" & Moisture_treatment == "joint" & cycle == "short")
plot(myres$delta_MLE~myres$Seed_weight_g, col=fagus.col[as.numeric(factor(myres$ID))], cex=2, pch=20,
     ylab="", xlab="Seed weight")
abline(lm(myres$delta_MLE~myres$Seed_weight_g), col="grey", lty=2)
mtext(paste("r =", round(cor.test(myres$delta_MLE, myres$Seed_weight_g)$estimate,2)), 3, line=-1, adj=1, cex=0.9)
mtext(paste("p-value =", round(cor.test(myres$delta_MLE, myres$Seed_weight_g)$p.value,3)), 3, line=-2, adj=1, cex=0.9)

dev.off()

## source("myplotcorr.r")
## library(ellipse)
## x = cor(subset(all.res, Genus == "Abies" & Moisture_treatment == "joint" & cycle == "short")[, germ.vars])
## my.plotcorr(x, col = terrain.colors(length(x[1,]))[5*x + 6], upper.panel="number", diag='ellipse')

all.res.sub = subset(all.res, Moisture_treatment %in% c("moist", "dry"))
all.res.sub$Moisture_treatment = factor(all.res.sub$Moisture_treatment, ordered=T, levels=c("moist", "dry"))

## plot effect of treatment
abies.g = ggplot(data=subset(all.res.sub, cycle == "short" & Genus=="Abies"),
       aes(y=GermPercent,x=Moisture_treatment, group = ID, colour = ID)) + xlab("") + ylab("") + 
    geom_line(linewidth = 1) + theme_bw()+ theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = abies.col) + ggtitle("Germination %")
abies.t50 = ggplot(data=subset(all.res.sub, cycle == "short" & Genus=="Abies"),
       aes(y=t50,x=Moisture_treatment, group = ID, colour = ID)) + xlab("") + ylab("") + 
    geom_line(linewidth = 1) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = abies.col) + ggtitle("t50")
abies.gu = ggplot(data=subset(all.res.sub, cycle == "short" & Genus=="Abies"),
       aes(y=GermUncertainty,x=Moisture_treatment, group = ID, colour = ID)) + xlab("") + ylab("") + 
    geom_line(linewidth = 1) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = abies.col) + ggtitle("Germination\nUncertainty")

fagus.g = ggplot(data=subset(all.res.sub, cycle == "short" & Genus=="Fagus"),
       aes(y=GermPercent,x=Moisture_treatment, group = ID, colour = ID)) + xlab("") + ylab("") + 
    geom_line(linewidth = 1) + theme_bw()+ theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = fagus.col) + ggtitle("Germination %")
fagus.t50 = ggplot(data=subset(all.res.sub, cycle == "short" & Genus=="Fagus"),
       aes(y=t50,x=Moisture_treatment, group = ID, colour = ID)) + xlab("") + ylab("") + 
    geom_line(linewidth = 1) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = fagus.col) + ggtitle("t50")
fagus.gu = ggplot(data=subset(all.res.sub, cycle == "short" & Genus=="Fagus"),
       aes(y=GermUncertainty,x=Moisture_treatment, group = ID, colour = ID)) + xlab("") + ylab("") + 
    geom_line(linewidth = 1) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = fagus.col) + ggtitle("Germination\nUncertainty")


## pdf("GermMetrics_Treatment.pdf",8,6)
## plot_grid(abies.g, abies.t50, abies.gu, fagus.g, fagus.t50, fagus.gu, nrow = 2, align = "hv", rel_heights = c(1/2, 1/2, 1/2))
## dev.off()

## ## plot effect of cycle
## all.res.sub = subset(all.res, Moisture_treatment %in% c("joint"))
## all.res.sub$Moisture_treatment = factor(all.res.sub$Moisture_treatment, ordered=T, levels=c("moist", "dry"))

## abies.g = ggplot(data=subset(all.res.sub, Genus=="Abies"),
##        aes(y=GermPercent,x=cycle, group = ID, colour = ID)) + xlab("") + ylab("") + 
##     geom_line(linewidth = 1) + theme_bw()+ theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = abies.col) + ggtitle("Germination %")
## abies.t50 = ggplot(data=subset(all.res.sub, Genus=="Abies"),
##        aes(y=t50,x=cycle, group = ID, colour = ID)) + xlab("") + ylab("") + 
##     geom_line(linewidth = 1) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = abies.col) + ggtitle("t50")
## abies.gu = ggplot(data=subset(all.res.sub, Genus=="Abies"),
##        aes(y=GermUncertainty,x=cycle, group = ID, colour = ID)) + xlab("") + ylab("") + 
##     geom_line(linewidth = 1) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = abies.col) + ggtitle("Germination Uncertainty")

## fagus.g = ggplot(data=subset(all.res.sub, Genus=="Fagus"),
##        aes(y=GermPercent,x=cycle, group = ID, colour = ID)) + ylab("") + 
##     geom_line(linewidth = 1) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = fagus.col)
## fagus.t50 = ggplot(data=subset(all.res.sub, Genus=="Fagus"),
##        aes(y=t50,x=cycle, group = ID, colour = ID)) + ylab("") + 
##     geom_line(linewidth = 1) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = fagus.col)
## fagus.gu = ggplot(data=subset(all.res.sub, Genus=="Fagus"),
##        aes(y=GermUncertainty,x=cycle, group = ID, colour = ID)) + ylab("") + 
##     geom_line(linewidth = 1) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = fagus.col)

## pdf("GermMetrics_Cycle.pdf",8,6)
## plot_grid(abies.g, abies.t50, abies.gu, fagus.g, fagus.t50, fagus.gu, nrow = 2, align = "hv", rel_heights = c(1/2, 1/2, 1/2))
## dev.off()

## KM survival plots
## ##################
datsurv = read.csv("dat_surv.csv")
datsurv$germinated.or.not = datsurv$germinated.or.not + 1

## merge with seed data
seed1 = read.csv("seed_provenance_details.csv")
datsurv = merge(datsurv, seed1, by="Provenance", all.x=T)

## exclude provenances with too little germination - ITF1
datsurv <- subset(datsurv, ! ID %in% c("FS_CH1", "FS_RO1", "FS_IT1"))

## exclude long cycle
##datsurv <- subset(datsurv, Germination.day < 73)

datsurv.abies = subset(datsurv, Genus == "Abies")
KM.abies <- survfit(Surv(Germination.day, germinated.or.not) ~ ID, data = datsurv.abies)

abies.surv.plot <- ggsurvplot(KM.abies, data = datsurv.abies, palette = abies.col, ylab="Probability of not germinating",
    xlab="Time (Days)", break.time.by = 10, ggtheme = theme_light(), linewidth= 1.5, ylim=c(0.3,1))

datsurv.fagus = subset(datsurv, Genus == "Fagus")
KM.fagus <- survfit(Surv(Germination.day, germinated.or.not) ~ ID, data = datsurv.fagus)

fagus.surv.plot <- ggsurvplot(KM.fagus, data = datsurv.fagus, palette = fagus.col, ylab="Probability of not germinating",
    xlab="Time (Days)", break.time.by = 10, ggtheme = theme_light(), linewidth= 1.5, ylim=c(0.3,1))

## pdf("CC_KMplot.pdf", 5,5)
## print(abies.surv.plot)
## print(fagus.surv.plot)
## ##plot_grid(abies.surv.plot, fagus.surv.plot, nrow = 2, align = "v", rel_heights = c(1/2, 1/2))
## dev.off()

## CPH model
## #############
datsurv = read.csv("dat_surv.csv")

## merge with seed data
seed1 = read.csv("seed_provenance_details.csv")
datsurv = merge(datsurv, seed1, by=c("Species", "Provenance"), all.x=T)

## exclude provenances with too little germination - ITF1
datsurv <- subset(datsurv, ! ID %in% c("FS_CH1", "FS_RO1", "FS_IT1"))

## exclude long cycle
##datsurv <- subset(datsurv, Germination.day < 73)

datsurv.abies = subset(datsurv, Genus == "Abies")
datsurv.fagus = subset(datsurv, Genus == "Fagus")

## order ID by germination %
abies.order = unique(subset(all.res.sub[order(all.res.sub$t50),], Genus=="Abies")$ID)
fagus.order = unique(subset(all.res.sub[order(all.res.sub$t50),], Genus=="Fagus")$ID)
datsurv.abies$ID = factor(datsurv.abies$ID, levels=unique(datsurv.abies$ID[order(datsurv.abies$Longitude)])) ## west to east order
datsurv.fagus$ID = factor(datsurv.fagus$ID, levels=unique(datsurv.fagus$ID[order(datsurv.fagus$Longitude)])) ## west to east order
datsurv.abies$Treatment = factor(datsurv.abies$Treatment, levels=c("moist", "dry"))
datsurv.fagus$Treatment = factor(datsurv.fagus$Treatment, levels=c("moist", "dry"))

## abies
cph.abies <- coxph(Surv(Germination.day, germinated.or.not, type = c("right")) ~ ID + Treatment, data = datsurv.abies)
modid.abies = forestmodel::forest_model(cph.abies)

cph.abies.cov <- coxph(Surv(Germination.day, germinated.or.not, type = c("right")) ~ Longitude + Moisture_content_percent + Seed_weight_g, data = datsurv.abies)
modcov.abies = forestmodel::forest_model(cph.abies.cov)

## fagus
cph.fagus <- coxph(Surv(Germination.day, germinated.or.not, type = c("right")) ~ ID + Treatment, data = datsurv.fagus)
modid.fagus = forestmodel::forest_model(cph.fagus)

cph.fagus.cov <- coxph(Surv(Germination.day, germinated.or.not, type = c("right")) ~ Longitude + Moisture_content_percent + Seed_weight_g, data = datsurv.fagus)
modcov.fagus = forestmodel::forest_model(cph.fagus.cov)

## pdf("CPH_all.pdf", 11, 15)
## plot_grid(modid.abies, modcov.abies, modid.fagus, modcov.fagus, nrow = 4, ncol=1, align="hv", rel_heights = c(1, .37, .6, .37))
## dev.off()


############################

## # Install the GitHub (development) version of forestmodel
## if (!requireNamespace("devtools", quietly = TRUE)) {
##   install.packages("devtools")
## }
## devtools::install_github("NikNakk/forestmodel")

library(forestmodel)
library(survival)
library(grid)
library(gridExtra)
library(gtable)

# formatting options
fmt <- forest_model_format_options(
  colour = "black",
  shape = 15,
  point_size = 3,
  text_size = 4,
  banded = FALSE
)

p_abies_id <- forest_model(
  cph.abies,
  exponentiate = TRUE,
  format_options = fmt,
  theme = theme_forest()
)

p_abies_cov <- forest_model(
  cph.abies.cov,
  exponentiate = TRUE,
  format_options = fmt,
  theme = theme_forest()
)

p_fagus_id <- forest_model(
  cph.fagus,
  exponentiate = TRUE,
  format_options = fmt,
  theme = theme_forest()
)

p_fagus_cov <- forest_model(
  cph.fagus.cov,
  exponentiate = TRUE,
  format_options = fmt,
  theme = theme_forest()
)

# 2. Function to align column widths between a set of forestmodel plots -----

align_widths <- function(plot_list) {
  gt <- lapply(plot_list, ggplotGrob)
  maxw <- do.call(grid::unit.pmax, lapply(gt, function(x) x$widths))
  lapply(gt, function(x) { x$widths <- maxw; x })
}

# 1. Align Abies + Fagus plots inside species columns
abies_aligned <- align_widths(list(p_abies_id, p_abies_cov))
fagus_aligned <- align_widths(list(p_fagus_id, p_fagus_cov))

# 2. Create species titles
title_abies <- textGrob("", gp = gpar(fontsize = 16, fontface = "bold"))
title_fagus <- textGrob("", gp = gpar(fontsize = 16, fontface = "bold"))

# 3. Two rows with different heights
row1 <- arrangeGrob(
  abies_aligned[[1]],
  fagus_aligned[[1]],
  ncol = 2
)

row2 <- arrangeGrob(
  abies_aligned[[2]],
  fagus_aligned[[2]],
  ncol = 2
)

# 4. Titles row
titles <- arrangeGrob(
  title_abies,
  title_fagus,
  ncol = 2
)

# 5. Now assemble with unequal heights
final <- arrangeGrob(
  titles,
  row1,
  row2,
  ncol = 1,
  heights = unit.c(
    unit(1, "lines"),  # small title area
    unit(2, "null"),   # row 1 gets 2/3 of height
    unit(1, "null")    # row 2 gets 1/3 of height
  )
)



## pdf("CPH_all.pdf", 14, 6)
## grid.newpage()
## grid.draw(final)
## dev.off()


### ===================================================================
### FINAL COMBINED FIGURE
### ===================================================================

## extract ggplot from survplot
abies.surv.plot <- ggsurvplot(
  KM.abies,
  data = datsurv.abies,
  palette = abies.col,
  ylab = "Probability of not germinating",
  xlab = "Time (Days)",
  break.time.by = 10,
  ggtheme = theme_light(),
  size = 1.5,
  ylim = c(0.3, 1)
)$plot

fagus.surv.plot <- ggsurvplot(
  KM.fagus,
  data = datsurv.fagus,
  palette = fagus.col,
  ylab = "Probability of not germinating",
  xlab = "Time (Days)",
  break.time.by = 10,
  ggtheme = theme_light(),
  size = 1.5,
  ylim = c(0.3, 1)
)$plot

abies.surv.plot <- abies.surv.plot +
  theme(legend.title = element_blank())

fagus.surv.plot <- fagus.surv.plot +
  theme(legend.title = element_blank())


## supress legend from treatment effect plots
abies.surv.plot    <- abies.surv.plot    + theme(legend.position = "none")
fagus.surv.plot    <- fagus.surv.plot    + theme(legend.position = "none")

## add letters
tag_plot <- function(p, letter) {
  p +
    labs(tag = letter) +
    theme(
      plot.tag = element_text(
        face = "bold",
        size = 15       # ⬅ increase from ~9
      ),
      plot.tag.position = c(0.01, 0.99)
    )
}

abies.surv.plot <- tag_plot(abies.surv.plot, "a")
abies.g    <- tag_plot(abies.g,    "b")
##abies.t50 <- tag_plot(abies.t50, "")
##abies.gu  <- tag_plot(abies.gu,  "")
p_abies_id  <- tag_plot(p_abies_id,  "c")

fagus.surv.plot <- tag_plot(fagus.surv.plot, "d")
fagus.g    <- tag_plot(fagus.g,    "e")
##fagus.t50 <- tag_plot(fagus.t50, "")
##fagus.gu  <- tag_plot(fagus.gu,  "")
p_fagus_id  <- tag_plot(p_fagus_id,  "f")

##p_abies_cov <- tag_plot(p_abies_cov, "e")
##p_fagus_cov <- tag_plot(p_fagus_cov, "l")

library(patchwork)

second_row_abies <-
  (abies.g | abies.t50 | abies.gu | guide_area()) +
  plot_layout(
    widths = c(1, 1, 1, 1),
    guides = "collect"
  )
second_row_fagus <-
  (fagus.g | fagus.t50 | fagus.gu | guide_area()) +
  plot_layout(
    widths = c(1, 1, 1, 1),
    guides = "collect"
  )

row_heights <- c(
  2.5,  # survival
  1.2,  # small 3-panel row
  2.4,  # forest plot (ID)
  .9   # forest plot (covariates)
)

abies_column <-
  (abies.surv.plot /
   second_row_abies /
   p_abies_id /
   p_abies_cov) +
  plot_layout(heights = row_heights)

fagus_column <-
  (fagus.surv.plot /
   second_row_fagus /
   p_fagus_id /
   p_fagus_cov) +
    plot_layout(heights = row_heights)

final_plot <- abies_column | fagus_column

ggsave(
  filename = "/home/kati/Dropbox/Overleaf/MyGardenOfTrees_PilotPhase_MS/New_Figures/Figure_results_cc.pdf",
  plot     = final_plot,
  device   = cairo_pdf,   # recommended for journals
  width    = 36,
  height   = 31,
  units    = "cm"
)


