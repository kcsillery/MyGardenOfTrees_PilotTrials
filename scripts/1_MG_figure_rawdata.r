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
library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)

## read data
load(file="dat2022.RData")

dat.all <- subset(dat.all, Date <= "2022-09-30") 

mycols <- read.csv("seed_dat_cols.csv")

## Abies colors
col_Abies <- mycols$col[mycols$Genus == "Abies"]
names(col_Abies) <- mycols$ID[mycols$Genus == "Abies"]

# Fagus colors
col_Fagus <- mycols$col[mycols$Genus == "Fagus"]
names(col_Fagus) <- mycols$ID[mycols$Genus == "Fagus"]


## Order gardens by mean GDD so facets are meaningful
garden.order <- with(subset(dat.all, Month %in% c(1:6)), ## take only month till June
                            tapply(GDD, Garden_ID, mean, na.rm = TRUE))
garden.order <- names(sort(garden.order))

dat.all$Garden_ID <- factor(dat.all$Garden_ID, levels = garden.order)

## create garden label
garden.meta <- unique(dat.all[c("Garden_ID", "Country_garden")])
garden.meta$Garden_label <- paste(paste("Garden", garden.meta$Garden_ID),
                                  " (", garden.meta$Country_garden, ")", sep="")
dat.all <- merge(dat.all, garden.meta, by = "Garden_ID", all.x = TRUE)

lab.vec <- garden.meta$Garden_label
names(lab.vec) <- garden.meta$Garden_ID
lab.vec <- lab.vec[lab.vec != "Garden 13 (NA)"]

## separate data for species
datA <- subset(dat.all, Genus == "Abies")
datF <- subset(dat.all, Genus == "Fagus")

## monotonous smoothing function
mono_smooth <- function(df) {
  df <- df[order(df$GDD), ]
  iso <- isoreg(df$GDD, df$germ.cum)

  tibble(
    GDD = df$GDD,
    germ.cum.mono = iso$yf,
    ID = df$ID,
    Garden_ID = df$Garden_ID
  )
}

datA.mono <- datA %>%
  group_by(Garden_ID, ID) %>%
  group_modify(~ mono_smooth(.x)) %>%
  ungroup()

datF.mono <- datF %>%
  group_by(Garden_ID, ID) %>%
  group_modify(~ mono_smooth(.x)) %>%
  ungroup()

nA.dates <- datA %>%
  group_by(Garden_ID) %>%
  summarise(dates = n_distinct(Date)) %>%
  mutate(x = Inf, y = Inf)

infoA <- datA %>%
  group_by(Garden_ID) %>%
  summarise(
    mean.GDD = mean(GDD, na.rm = TRUE),
    first.date = min(Date),
    last.date  = max(Date)
  ) %>%
  mutate(x = Inf, y = Inf)

nF.dates <- datF %>%
  group_by(Garden_ID) %>%
  summarise(dates = n_distinct(Date)) %>%
  mutate(x = Inf, y = Inf)

infoF <- datF %>%
  group_by(Garden_ID) %>%
  summarise(
    mean.GDD = mean(GDD, na.rm = TRUE),
    first.date = min(Date),
    last.date  = max(Date)
  ) %>%
  mutate(x = Inf, y = Inf)


### ============================================================
### theme
### ============================================================

theme_pub <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "grey60"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey40", fill = NA),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")
  ) 

## plot Abies

plot.Abies <- ggplot(datA) +

  geom_point(
    aes(x = GDD, y = germ.cum, color = ID),
    size = 0.9, alpha = 0.35
  ) +

  geom_line(
    data = datA.mono,
    aes(x = GDD, y = germ.cum.mono, color = ID,
        group = interaction(Garden_ID, ID)),
    linewidth = 0.8
  ) +

  geom_vline(
    data = infoA,
    aes(xintercept = 500),
    color = "black",
    linewidth = 0.4,
    linetype = "dashed"
  ) +

  # annotations (stacked)
    geom_text(data = nA.dates,
            aes(x = -Inf, y = Inf, label = paste0("N observations = ", dates)),
            hjust = -0.1, vjust = 2, size = 3, inherit.aes = FALSE) +

  geom_text(data = infoA,
            aes(x = -Inf, y = Inf, label = paste0("first: ", first.date)),
            hjust = -0.1, vjust = 3.5, size = 3, inherit.aes = FALSE) +

  geom_text(data = infoA,
            aes(x = -Inf, y = Inf, label = paste0("last: ", last.date)),
            hjust = -0.1, vjust = 5, size = 3, inherit.aes = FALSE) +

  scale_color_manual(values = col_Abies) +
    
    facet_wrap(~ Garden_ID, ncol = 4, scales = "free_x",
               labeller = labeller(Garden_ID = lab.vec)) +

  labs(
    x = "Growing degree-days (°C days)",
    y = "Cumulative germination (%)",
    title = expression(paste("a. Abies")),
    color = "Provenance"
  ) +

  guides(color = guide_legend(nrow = 2)) +

    theme_pub +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )



## ####################################



plot.Fagus <- ggplot(datF) +

  geom_point(
    aes(x = GDD, y = germ.cum, color = ID),
    size = 0.9, alpha = 0.35
  ) +

  geom_line(
    data = datF.mono,
    aes(x = GDD, y = germ.cum.mono, color = ID,
        group = interaction(Garden_ID, ID)),
    linewidth = 0.8
  ) +

  geom_vline(
    data = infoF,
    aes(xintercept = 500),
    color = "black",
    linewidth = 0.4,
    linetype = "dashed"
  ) +


  geom_text(data = nF.dates,
            aes(x = -Inf, y = Inf, label = paste0("N observations = ", dates)),
            hjust = -0.1, vjust = 2, size = 3, inherit.aes = FALSE) +

  geom_text(data = infoF,
            aes(x = -Inf, y = Inf, label = paste0("first: ", first.date)),
            hjust = -0.1, vjust = 3.5, size = 3, inherit.aes = FALSE) +

  geom_text(data = infoF,
            aes(x = -Inf, y = Inf, label = paste0("last: ", last.date)),
            hjust = -0.1, vjust = 5, size = 3, inherit.aes = FALSE) +

  scale_color_manual(values = col_Fagus) +

    facet_wrap(~ Garden_ID, ncol = 4, scales = "free_x",
               labeller = labeller(Garden_ID = lab.vec)) +

  labs(
    x = "Growing degree-days (°C days)",
    y = "Cumulative germination (%)",
    title = expression(paste("b. Fagus")),
    color = "Provenance"
  ) +

  guides(color = guide_legend(nrow = 2)) +

  theme_pub +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )


ggsave("Abies_monotone_GDD.pdf", plot.Abies, width = 10, height = 12, dpi = 300)
ggsave("Fagus_monotone_GDD.pdf", plot.Fagus, width = 10, height = 12, dpi = 300)
