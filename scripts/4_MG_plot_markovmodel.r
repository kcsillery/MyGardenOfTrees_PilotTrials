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
library(dplyr)
library(ggsci)
library(scales)
library(RColorBrewer)

source("00_plot.r")

mle.a <- read.table("Abies_MLE.txt", header = T)
mle.f <- read.table("Fagus_MLE.txt", head=T)

# plot as matrix
plot(mle.a[,c("g", "alpha", "beta", "delta", "g_without_death", "p_survive_until_x_1", "p_survive_until_x_5", "p_survive_until_x_7")],
     col=factor(mle.a$garden))
plot(mle.a[,c("g", "alpha", "beta", "delta", "g_without_death", "p_survive_until_x_1", "p_survive_until_x_3", "p_survive_until_x_7")],
     col=as.factor(mle.a$provenance))

plot(mle.f[,c("g", "alpha", "beta", "delta", "g_without_death", "p_survive_until_x_1", "p_survive_until_x_5", "p_survive_until_x_7")],
     col=factor(mle.f$garden))
plot(mle.f[,c("g", "alpha", "beta", "delta", "g_without_death", "p_survive_until_x_1", "p_survive_until_x_3", "p_survive_until_x_7")],
     col=as.factor(mle.f$provenance))

## g and delta are uncorrelated, but Px is correlated to both
## more clustering based on garden then provenance in abies, not much clustering either way for fagus

## combine results with seed data
mle.a$Genus = "Abies"
mle.f$Genus = "Fagus"
seed.dat = read.csv("seed_dat_cols.csv")
seed.dat$provenance = seed.dat$Provenance
dat1 = merge(mle.a, seed.dat, by=c("provenance", "Genus"), all.x=T)
dat2 = merge(mle.f, seed.dat, by=c("provenance", "Genus"), all.x=T)
dat = rbind(dat1, dat2)

## remove, not germinated
dat = subset(dat, !ID %in% c("FS_CH1", "FS_RO1", "AA_AT1"))

dat$Garden_ID <- as.character(as.numeric(unlist(lapply(strsplit(dat$garden, split="_"), function(a) a[1]))))
dat$Block <- gsub("21", "", unlist(lapply(strsplit(dat$garden, split="_"), function(a) a[2])))
dat$Block[is.na(dat$Block)] <- ""
dat$Garden_ID21 <- paste(dat$Garden_ID, dat$Block, sep="")

## combine with climate chamber
load("../climate_chamber/mm_dat.RData")
all.dat <- merge(dat, mm.dat[c("ID", "g_MLE", "gamma_MLE", "delta_MLE", "alpha_MLE", "beta_MLE")], by="ID", all=T)

## order provenances by germ
all.dat.abies = subset(all.dat, Genus == "Abies")
all.dat.fagus = subset(all.dat, Genus == "Fagus")

myord.p = aggregate(all.dat$g_without_death, list(all.dat$Genus, all.dat$ID), mean, na.rm=T)
abies.ord.p = myord.p$Group.2[myord.p$Group.1=="Abies"][order(myord.p$x[myord.p$Group.1=="Abies"])]
fagus.ord.p = myord.p$Group.2[myord.p$Group.1=="Fagus"][order(myord.p$x[myord.p$Group.1=="Fagus"])]

all.dat.abies$ID = factor(all.dat.abies$ID, ordered=T, levels=abies.ord.p)
all.dat.fagus$ID = factor(all.dat.fagus$ID, ordered=T, levels=fagus.ord.p)

## order gardens by garden
myord.g = aggregate(all.dat$g_without_death, list(all.dat$Genus, all.dat$Garden_ID21), mean, na.rm=T)
abies.ord = myord.g$Group.2[myord.g$Group.1=="Abies"][order(myord.g$x[myord.g$Group.1=="Abies"], decreasing=T)]
fagus.ord = myord.g$Group.2[myord.g$Group.1=="Fagus"][order(myord.g$x[myord.g$Group.1=="Fagus"], decreasing=T)]

all.dat.abies$Garden_ID21 = factor(all.dat.abies$Garden_ID21, ordered=T, levels=abies.ord)
all.dat.fagus$Garden_ID21 = factor(all.dat.fagus$Garden_ID21, ordered=T, levels=fagus.ord)

all.dat.abies = all.dat.abies[order(all.dat.abies$ID, all.dat.abies$Garden_ID21), ]
all.dat.fagus = all.dat.fagus[order(all.dat.fagus$ID, all.dat.fagus$Garden_ID21), ]


## print equation
print.eq <- function(fit, df){
    cf <- round(coef(fit), 2) 
    eq <- paste0("y = ", cf[1],
                 ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x")
    pval <- paste0("p-value = ", round(summary(fit)$coeff[2,4], 4))
    mtext(eq, 3, line=-2)
    mtext(pval, 3, line=-3)
}

################################################
pdf("Figure_results_mg.pdf", 18, 10)
par(mfrow=c(2,4), bty="l", cex.lab=1.5, cex.axis=1.4, mar=c(7,6,5,0))

## abies
par(xaxt="n")
tmp <- na.omit(all.dat.abies[, c("g_without_death", "g_MLE", "ID", "Garden_ID21","never_germinated", "col")])
plot(tmp$g_without_death ~ tmp$g_MLE, col=tmp$col, pch=20,
     xlab="", ylab="")
tmp2 <- unique(tmp[, c("ID", "g_MLE")])
tmp2 <- tmp2[order(tmp2$g_MLE),]
par(xaxt="s")
axis(side = 1, at =tmp2$g_MLE+c(0,0,0,0,0,-0.008,0.008,0,0), labels=tmp2$ID, las=2)
axis(side = 3)
mtext("Germination rate (climate chamber)", 3, 3, cex=1.2)
mtext("Germination rate (field trials)", 2, 3, cex=1.2)

fit <- lm(g_without_death ~ g_MLE, dat=tmp)
newx <- seq(min(tmp$g_MLE), max(tmp$g_MLE), length.out=nrow(tmp))
preds <- predict(fit, newdata = data.frame(g_MLE=newx), interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey90', border = NA)
abline(fit)
points(tmp$g_without_death ~ tmp$g_MLE, col=tmp$col, pch=20, cex=2.5)
print.eq(fit)
abline(0,1, lty=2)
mtext("a", line=2.3, adj=0, cex=1.8)

tmp <- all.dat.abies
data <- data.frame(garden = tmp$ID, 
                   provenance = tmp$Garden_ID21,
                   val_to_plot = tmp$g_without_death, 
                   never_germinated = tmp$never_germinated,
                   mycol = tmp$col)
makeImageWithMeans(data, stat = mean)
mtext("Garden ID", 2, 3.5)
mtext("b", line=2.3, adj=0, cex=1.8)

## fagus
par(xaxt="n")
tmp <- na.omit(all.dat.fagus[, c("g_without_death", "g_MLE", "ID", "Garden_ID21","never_germinated", "col")])
plot(tmp$g_without_death ~ tmp$g_MLE, col=tmp$col, pch=20, xlim=c(0,.7),
     xlab="", ylab="")
tmp2 <- unique(tmp[, c("ID", "g_MLE")])
tmp2 <- tmp2[order(tmp2$g_MLE),]
par(xaxt="s")
axis(side = 1, at =tmp2$g_MLE+c(0,-0.008,0.008,0,0), labels=tmp2$ID, las=2)
axis(side = 3)
mtext("Germination rate (climate chamber)", 3, 3, cex=1.2)
mtext("Germination rate (field trials)", 2, 3, cex=1.2)

fit <- lm(g_without_death ~ g_MLE, dat=tmp)
newx <- seq(min(tmp$g_MLE), max(tmp$g_MLE), length.out=nrow(tmp))
preds <- predict(fit, newdata = data.frame(g_MLE=newx), interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey90', border = NA)
abline(fit)
points(tmp$g_without_death ~ tmp$g_MLE, col=tmp$col, pch=20, cex=2.5)
print.eq(fit)
abline(0,1, lty=2)
mtext("c", line=2.3, adj=0, cex=1.8)

tmp <- all.dat.fagus
data <- data.frame(garden = tmp$ID, 
                   provenance = tmp$Garden_ID21,
                   val_to_plot = tmp$g_without_death, 
                   never_germinated = tmp$never_germinated,
                   mycol = tmp$col)
makeImageWithMeans(data, stat = mean)
mtext("Garden ID", 2, 3.5)
mtext("d", line=2.3, adj=0, cex=1.8)

## pdf("Figure_resultsFT-CC_delta.pdf", 9, 10)
## par(mfrow=c(2,2), bty="l", cex.lab=1.2, mar=c(6,5,4,0))

## abies
par(xaxt="n")
tmp <- na.omit(all.dat.abies[, c("delta", "delta_MLE", "ID", "Garden_ID21","never_germinated", "col")])
plot(tmp$delta ~ tmp$delta_MLE, col=tmp$col, pch=20,
     xlab="", ylab="")
tmp2 <- unique(tmp[, c("ID", "delta_MLE")])
tmp2 <- tmp2[order(tmp2$delta_MLE),]
par(xaxt="s")
axis(side = 1, at =tmp2$delta_MLE, labels=tmp2$ID, las=2)
axis(side = 3)
mtext("Development speed (climate chamber)", 3, 3, cex=1.2)
mtext("Development speed (field trials)", 2, 3, cex=1.2)


fit <- lm(delta ~ delta_MLE, dat=tmp)
newx <- seq(min(tmp$delta_MLE), max(tmp$delta_MLE), length.out=nrow(tmp))
preds <- predict(fit, newdata = data.frame(delta_MLE=newx), interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey90', border = NA)
abline(fit)
points(tmp$delta ~ tmp$delta_MLE, col=tmp$col, pch=20, cex=2.5)
print.eq(fit)
abline(0,1, lty=2)
mtext("e", line=2.3, adj=0, cex=1.8)

tmp <- all.dat.abies
tmp$col[is.na(tmp$delta)] = NA
data <- data.frame(garden = tmp$ID, 
                   provenance = tmp$Garden_ID21,
                   val_to_plot = tmp$delta, 
                   never_germinated = tmp$never_germinated,
                   mycol = tmp$col)
makeImageWithMeans(data, stat = mean)
mtext("Garden ID", 2, 3.5)
mtext("f", line=2.3, adj=0, cex=1.8)

## fagus
par(xaxt="n")
tmp <- na.omit(all.dat.fagus[, c("delta", "delta_MLE", "ID", "Garden_ID21","never_germinated", "col")])
plot(tmp$delta ~ tmp$delta_MLE, col=tmp$col, pch=20,
     xlab="", ylab="")
tmp2 <- unique(tmp[, c("ID", "delta_MLE")])
tmp2 <- tmp2[order(tmp2$delta_MLE),]
par(xaxt="s")
axis(side = 1, at =tmp2$delta_MLE, labels=tmp2$ID, las=2)
axis(side = 3)
mtext("Development speed (climate chamber)", 3, 3, cex=1.2)
mtext("Development speed (field trials)", 2, 3, cex=1.2)

fit <- lm(delta ~ delta_MLE, dat=tmp)
newx <- seq(min(tmp$delta_MLE), max(tmp$delta_MLE), length.out=nrow(tmp))
preds <- predict(fit, newdata = data.frame(delta_MLE=newx), interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey90', border = NA)
abline(fit)
points(tmp$delta ~ tmp$delta_MLE, col=tmp$col, pch=20, cex=2.5)
print.eq(fit)
abline(0,1, lty=2)
mtext("g", line=2.3, adj=0, cex=1.8)

tmp <- all.dat.fagus
tmp$col[is.na(tmp$delta)] = NA
data <- data.frame(garden = tmp$ID, 
                   provenance = tmp$Garden_ID21,
                   val_to_plot = tmp$delta, 
                   never_germinated = tmp$never_germinated,
                   mycol = tmp$col)
makeImageWithMeans(data, stat = mean)
mtext("Garden ID", 2, 3.5)
mtext("h", line=2.3, adj=0, cex=1.8)

dev.off()

