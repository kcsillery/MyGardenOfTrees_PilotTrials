## ##########################################################################

## Madleina Caduff, with slight modifications from Katalin Csillery

## 28 Jan 2026

## For the manuscript:
## Gene–environment interactions govern early regeneration in fir and
## beech: evidence from participatory provenance trials across Europe
## by Katalin Csilléry, Justine Charlet de Sauvage, Madleina Caduff,
## Johannes Alt, Marjorie Bison, Mert Celik, Nicole Ponta, Daniel
## Wegmann

## ###########################################################################



##--------------------------------
# Functions to read data
#--------------------------------

get_garden_from_class <- function(classes){
  return(sapply(classes, function(x) strsplit(x, "_")[[1]][1]))
}

get_provenance_from_class <- function(classes){
  return(sapply(classes, function(x) { tt <- strsplit(x, "_")[[1]]; paste0(tt[2:length(tt)], collapse = "_")}))
}

get_names_garden_provenance <- function(all_files, species){
  # Get class names and replicate
  tmp <- as.character(sapply(all_files, function(x) strsplit(x, "_rep")[[1]][1]))
  class <- as.character(sapply(tmp, function(x) strsplit(x, paste0("^", species, "_"))[[1]][2]))
  
  rep <- as.character(sapply(all_files, function(x) strsplit(strsplit(x, "rep_")[[1]][2], "[_\\.]")[[1]][1]))
  
  return(data.frame(f = all_files, rep = rep, class = class))
}

get_final_LL <- function(class, rep, path, species){
  name <- paste0(path, species, "_", class, "_rep_", rep, "_trace.txt")
  if (!file.exists(name) || file.info(name)$size == 0){
    return(NA)
  }
  r <- read.table(name, header = T)
  LL <- r$LL[nrow(r)]
  return(LL)
}

find_better_replicate <- function(c, path, all_files, species){
  LL1 <- get_final_LL(c, 1, path, species)
  LL2 <- get_final_LL(c, 2, path, species)
  if (!is.na(LL1) & !is.na(LL2) & abs(LL1-LL2) > 1){ print(paste(c, LL1, LL2, abs(LL1-LL2))) }
  if (is.na(LL1) & is.na(LL2)){ 
    # Both NA -> just keep one
    print(paste(c, LL1, LL2, abs(LL1-LL2)))
    return(all_files$f[all_files$class == c & all_files$rep == 1 & grepl("trace", all_files$f)])
  }
  
  # keep better
  r <- which.max(c(LL1, LL2))
  return(all_files$f[all_files$class == c & all_files$rep == r & grepl("trace", all_files$f)])
}

read_MLE <- function(all_files, species){
  num_classes <- length(unique(all_files$class))
  estimates <- data.frame(class = unique(all_files$class),
                          garden = get_garden_from_class(unique(all_files$class)),
                          provenance = get_provenance_from_class(unique(all_files$class)),
                          g = rep(NA, num_classes), 
                          alpha = rep(NA, num_classes), 
                          beta = rep(NA, num_classes),
                          delta = rep(NA, num_classes),
                          never_germinated = rep(NA, num_classes))
  
  for (c in unique(all_files$class)){
    name <- all_files$f[all_files$class == c]
    
    # if class germinated: read from trace file
    is_trace <- grepl("trace", name)
    if (any(is_trace)){
      # Get better replicate
      name <- find_better_replicate(c, path, all_files, species)
    } else {
      # Never germinated -> no trace is written. Set g and delta to 0
      name <- character()
      print(paste0(c, ": Never germinated"))
    }
    
    if (length(name) == 0 || file.info(paste0(path, name))$size == 0){
      estimates$g[estimates$class == c] <- 0
      estimates$delta[estimates$class == c] <- 0
      estimates$never_germinated[estimates$class == c] <- TRUE
    } else {
      f <- read.table(paste0(path, name), header = T)
      f <- f[nrow(f),]
      
      estimates$g[estimates$class == c] <- f[,grepl("^g_", names(f))]
      estimates$alpha[estimates$class == c] <- f[,grepl("^alpha_", names(f))]
      estimates$beta[estimates$class == c] <- f[,grepl("^beta_", names(f))]
      estimates$delta[estimates$class == c] <- f[,grepl("^delta_", names(f))]
      estimates$never_germinated[estimates$class == c] <- FALSE
    }
  }
  return(estimates)
}

#--------------------------------
# Functions to calculate
#--------------------------------

calc_prob_death <- function(x, alpha, beta){
  return(alpha * exp(-beta * x))
}

calc_prob_to_survive_until_x <- function(x, alpha, beta, g, delta){
  if (g == 0){ return(rep(0, length(x))) }
  
  num_days <- ceiling(x / delta)
  num_days[num_days > 10000] <- 10000 # cut to avoid huge values
  prob_survive_until_x <- g * sapply(num_days, function(dd) prod(1 - calc_prob_death(delta * 1:dd, alpha, beta)))
  
  return(prob_survive_until_x)
}

calc_prob_to_survive_until_x_classes <- function(MLE, x = 100){
  MLE$p_survive_until_x <- NULL
  for (i in 1:nrow(MLE)){
    MLE$p_survive_until_x[i] <- calc_prob_to_survive_until_x(x, 
                                              alpha = MLE$alpha[i],
                                              beta = MLE$beta[i],
                                              g = MLE$g[i],
                                              delta = MLE$delta[i])
  }
  return(MLE)
}

guess_germination_without_death <- function(MLE, species){
  obs_data <- read.table(paste0("../data/Input_obs_", species, ".txt"), header = T)
  
  if (species == "Abies"){
    stage_0 = "Stage_0"
  } else {
    stage_0 = "Counts_Unobserved"
  }
  
  MLE$g_without_death <- NA
  for (i in 1:nrow(MLE)){
    subset <- obs_data[obs_data$Provenance == MLE$provenance_full[i] & obs_data$Garden_ID == as.numeric(MLE$garden[i]),]
    max_germinated <- subset %>%
      group_by(Garden_block_spot) %>%
      summarise(num_seeds = max(.data[[stage_0]]),
                num_germinated = max(.data[[stage_0]]) - min(.data[[stage_0]])) %>% 
      ungroup()
    
    MLE$g_without_death[i] <- sum(max_germinated$num_germinated) / sum(max_germinated$num_seeds)
  }
  
  
  
  return(MLE)
}

#---------------------------------------
# Functions to plot data as image matrix 
#---------------------------------------

addPlotLetter <- function(num){
  mtext(letters[num],
        at = par("usr")[4], side = 2, padj = 1, adj = 1,
        line = 1, las = 1,
        font = 2)
}

# function for color
getCol <- function(x){
  my_colors <- brewer.pal(n = 9, name = "Blues")
  color_ramp <- colorRamp(my_colors)
  return(rgb(color_ramp(x) / 255))
}

getCol2 <- function(x, col){
  return(alpha(col, x))
}

addPolygon <- function(x, y, col, wx = 1, wy = 1, hatched = FALSE, crossed = FALSE){
  xx <- x + 0.5*wx * c(-1, 1, 1, -1, -1)
  yy <- y + 0.5 * wy * c(-1, -1, 1, 1, -1)
  if (hatched){
    polygon(xx, yy, col = col, xpd = TRUE, border = NA, density = 20, angle = 45)
  } else if (crossed) {
    cx <- mean(xx)
    cy <- mean(yy)
    d <- 0
    
    # draw small "X"
    segments(min(xx), min(yy), max(xx), max(yy), col = col, xpd = TRUE) # \
    segments(min(xx), max(yy), max(xx), min(yy), col = col, xpd = TRUE) # /
  } else {
    polygon(xx, yy, col = col, xpd = TRUE, border = NA)
  }
}

add_scale <- function(xlabs, zlim){
  # coordinates for scale bar
    ncols <- 1000
    scale_x <- seq(0.5, length(xlabs)+0.5, length.out = ncols)
    ##factor <- ifelse(species == "Abies", 8, 10.5)
    scale_y <- (diff(par("usr")[1:2]) )   # place below plot
    
  
    cols <- alpha("black", seq(0, 1, length.out = ncols))
    for (i in 1:(ncols-1)) {
        rect(scale_x[i], scale_y, scale_x[i+1], scale_y - 0.7, 
             col = cols[i], border = cols[i], xpd = T)
    }
    
    at <- seq(0.5, length(xlabs)+0.5, length.out = 5)
    labels <- round(seq(zlim[1], zlim[2], length.out = 5), 1)
    
                                        # add numeric axis
    axis(1, at = at,
         labels = labels,
         line = scale_y, tick = T, xpd = TRUE, outer = TRUE, col = NA, col.ticks = NA)
}

makeImageWithMeans <- function(data, stat = mean, addYLabels = TRUE){
                                        # Takes a data frame with three columns: x-labels, y-labels, values
                                        # get labels in order
    xlabs <- sort(unique(data[,1]))
    ylabs <- sort(unique(data[,2]), decreasing = TRUE)
    mycol <- data[,5]
    
                                        # plot as image
    zlim <- range(data[,3], na.rm = TRUE)
    plot(0, type='n', xlim = c(0.5, length(xlabs) + 0.5 + 1.5), ylim = c(0.5, length(ylabs) + 0.5 + 1.5), axes = FALSE, xlab = "", ylab = "", xaxs = 'i', yaxs = 'i')
    
                                        # add polygons for each entry
    for(i in 1:nrow(data)){
        x <- which(xlabs == data[i,1])
        y <- which(ylabs == data[i,2])
        addPolygon(x, y, col = alpha(mycol[i], (data[i,3] - zlim[1]) / (zlim[2] - zlim[1])))
    }
    
                                        # Add polygons for missing entries
    for(i in 1:length(xlabs)){
        for (j in 1:length(ylabs)){
            if (any(xlabs[i] == data[,1] & ylabs[j] == data[,2])){ next }
            addPolygon(i, j, col = "lightgrey", hatched = T)
        }
    }
    
                                        # Add polygons for entries that never germinated
    for(i in 1:nrow(data)){
        if (data$never_germinated[i]){ 
            x <- which(xlabs == data[i,1])
            y <- which(ylabs == data[i,2])
            addPolygon(x, y, col = "lightgrey", crossed = T)
        }
    }
    
                                        # add axis
    axis(side = 1, at = 1:length(xlabs), labels = xlabs, las = 2)
    if(addYLabels){
        axis(side = 2, at = 1:length(ylabs), labels = ylabs, las = 1)
    }
    
    ## add box
    gardens <- unique(data$probSurvive.garden)
    provenances <- unique(data$probSurvive.provenance)
    rect(0.5, 0.5, length(gardens)+0.5, length(provenances)+0.5, lwd = 1)

    coltoalpha = ifelse(is.na(mean(data[data[,1] == xlabs[i], 3], na.rm=T)), "#00000000", "black")
    
    ## add stat to edges
    for(i in 1:length(xlabs)){
        
        addPolygon(i, length(ylabs) + 1.2, col = alpha(coltoalpha, mean(data[data[,1] == xlabs[i], 3], na.rm=T)))
    }
    for(j in 1:length(ylabs)){
        addPolygon(length(xlabs) + 1.2, j, col = alpha(coltoalpha, mean(data[data[,2] == ylabs[j], 3], na.rm=T)))
    }
    ##add_scale(xlabs, zlim)
}

plot_P_survive_and_delta_matrix <- function(MLE, species){
                                        # open plot
    par(mfrow = c(1,3), oma = c(0, 13, 0, 0), mar = c(7, 1, 1, 1))
    
                                        # plot minimal germination rate
    data <- data.frame(garden = MLE$garden, 
                       provenance = MLE$provenance, 
                       val_to_plot = MLE$g_without_death, 
                       never_germinated = MLE$never_germinated)
    makeImageWithMeans(data, colfunc = getCol, stat = mean, species = species)
    mtext("Gardens", side = 1, line = 3)
    mtext("Provenances", side = 2, line = 12.5)
    addPlotLetter(1)
    
                                        # plot P_survive
    probSurvive <- calc_prob_to_survive_until_x_classes(MLE, 5)
    data <- data.frame(garden = probSurvive$garden, 
                       provenance = probSurvive$provenance, 
                       val_to_plot = probSurvive$p_survive_until_x, 
                       never_germinated = probSurvive$never_germinated)
    makeImageWithMeans(data, colfunc = getCol, stat = mean, species = species, addYLabels = FALSE)
    mtext("Gardens", side = 1, line = 3)
    addPlotLetter(2)
    
                                        # plot delta
    data <- data.frame(garden = MLE$garden, 
                       provenance = MLE$provenance, 
                       val_to_plot = MLE$delta, 
                       never_germinated = MLE$never_germinated)
    makeImageWithMeans(data, colfunc = getCol, stat = mean, species = species, addYLabels = FALSE)
    mtext("Gardens", side = 1, line = 3)
    addPlotLetter(3)
}
