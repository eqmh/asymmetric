####-----------------------------------------------------------------------------------------------
####                      OPTIMIZATION OF ROCKY SHORE SAMPLING FROM P2P                           #
####-----------------------------------------------------------------------------------------------

#####-----------------------------------------------------------------------------------------------

# Load required libraries
library(iNEXT)
library(tidyverse)

# setwd(path = "C:/Users/jslef/OneDrive - Smithsonian Institution/Documents/GitHub/asymmetric/")
# enrique change this for your local siles
setwd("~/asymmetric")

# Load function to compute coverage-based stopping
source("R/covstop.R")

# finch::dwca_read(..., read = TRUE)

## TO DO ##
# generate rarefaction curves
# also run collapsing by locality

#####-----------------------------------------------------------------------------------------------

# Read in data

files <- list.files("./data")

# # remove weird locality
# files <- files[-2]

p2p <- do.call(rbind, lapply(files, function(l) {
  
  # read in data
  dat <- read.csv(paste0("./data/", l))
  
  # take only the present species
  dat <- subset(dat, occurrenceStatus == "present")
  
  # split metadata
  meta <- do.call(rbind.data.frame, strsplit(as.character(dat$eventID), "\\_"))
  
  names(meta) <- c("country", "locality", "site", "date", "strata", "quadrat")
  
  taxa <- do.call(rbind.data.frame, strsplit(as.character(dat$scientificNameID), "\\:"))[, 5]
  
  # bind back in
  data <- cbind.data.frame(eventID = dat$eventID, meta, taxa)
  
  return(data)
  
} ) )

# Remove duplicate rows
p2p <- p2p %>% distinct(country, locality, site, strata, quadrat, taxa, .keep_all = TRUE) %>% ungroup()

# Write csv
write.csv(p2p, "./sites.csv")

# Read in site lat/longs
sites <- read.csv("./sites.csv", header = T)

#####-----------------------------------------------------------------------------------------------

# Generate interpolation/extrapolation curves for each locality (or site)

rare <- do.call(rbind, lapply(unique(p2p$locality), function(i) {
  
  x <- subset(p2p, locality == i)
  
  # do.call(rbind, lapply(unique(p2p$site), function(j) {
  
  # x <- subset(p2p, site == j)
    
    do.call(rbind, lapply(unique(p2p$strata), function(k) {
  
      print(paste(i, k))
      
      x <- subset(x, strata == k)
      
      if(nrow(x) == 0) data.frame() else {
        
        # summarize by locality
        x$presence <- 1
        
        # x <- x %>% group_by(locality, site, strata, taxa) %>% summarize(presence = sum(presence))
          
        # Cast longways
        mat <- x %>% select(site, quadrat, taxa, presence) %>% 
          
          pivot_wider(id_cols = c(site, quadrat), names_from = taxa, values_from = presence) 
        
        mat[is.na(mat)] <- 0
        
        dnames <- list(colnames(mat)[-(1:2)], as.character(mat$site))
        
        mat <- t(as.matrix(mat)[, -(1:2), drop = FALSE])
        
        mat <- apply(mat, 2, as.numeric)
        
        dimnames(mat) <- dnames
        
        # z <- c(sum(mat), apply(mat, 1, sum))
        
        z <- as.incfreq(mat)
        
        out <- iNEXT(z, datatype = "incidence_freq")
        
        ret <- out$iNextEst
        
        idx <- which(ret$method == "observed")
        
        ret <- rbind.data.frame(
          ret[ret$method == "interpolated", ],
          ret[idx, ],
          ret[idx, ],
          ret[idx, ],
          ret[ret$method == "extrapolated", ]
        )
        
        ret[idx, "method"] <- "interpolated"
        
        ret[idx + 2, "method"] <- "extrapolated"
      
        data.frame(
          locality = i,
          strata = k,
          ret[, c(1:2, 4)]
        )
        
      }
      
      } ) )
    
    } ) )
  
  # } ) )

rare$strata <- factor(rare$strata, levels = c("HIGHTIDE", "MIDTIDE", "LOWTIDE"))

# Plot results
(rareplot <- ggplot() +
  geom_line(data = subset(rare, method == "interpolated"), aes(x = t, y = qD, group = paste(locality, strata), col = strata)) +
  geom_line(data = subset(rare, method == "extrapolated"), aes(x = t, y = qD, group = paste(locality, strata), col = strata), lty = 3) +
  geom_point(data = subset(rare, method == "observed"), aes(x = t, y = qD, group = paste(locality, strata), col = strata), size = 2) +
  scale_color_manual(values = c("black", "dodgerblue3", "forestgreen")) +
  labs(x = "Number of samples", y = "Species richness") +
  facet_grid( ~ strata, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
)


ggsave("./output/Rarefaction plot.pdf", rareplot, device = "pdf", width = 10, height = 5, units = "in")

# Plot results with curves colored according to locality  
(rareplot_1 <- ggplot() +
  geom_line(data = subset(rare, method == "interpolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality)) + 
  geom_line(data = subset(rare, method == "extrapolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), lty = 3) + 
  geom_point(data = subset(rare, method == "observed" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), size = 2) + 
  ylim(0, 65) + 
  labs(x = "Number of samples", y = "Species richness") + 
  facet_grid( ~ strata, scales = "free_x") + 
  theme_bw(base_size = 14) +
  theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # legend.position = "none"
  )
)

ggsave("./output/Rarefaction plot 2.pdf", rareplot_1, device = "pdf", width = 10, height = 5, units = "in")  

### This plots individual sites

sel_locality = "NORTHERNMA"

(rareplot_1 <- ggplot() +
    geom_line(data = subset(rare, method == "interpolated" & strata == strata & locality == sel_locality), aes(x = t, y = qD, group = paste(strata))) + 
    geom_line(data = subset(rare, method == "extrapolated" & strata == strata & locality == sel_locality), aes(x = t, y = qD, group = paste(strata)), lty = 3) + 
    geom_point(data = subset(rare, method == "observed" & strata == strata & locality == sel_locality), aes(x = t, y = qD, group = paste(strata)), size = 2) + 
    ylim(0, 65) + 
    xlim(0, 80) + 
    labs(x = "Number of samples", y = "Species richness") + 
    facet_grid( ~ strata, scales = "free_x") + 
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # legend.position = "none"
    )
)

### Delete rows from specific localities in 'rare' data frame
rare2 <- rare[!grepl("MASSACHUSETTS", rare$locality),] ### this allowes to generate 'rare_plot_2'

# Plot results with curves colored according to locality  
(rareplot_1 <- ggplot() +
    geom_line(data = subset(rare2, method == "interpolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality)) + 
    geom_line(data = subset(rare2, method == "extrapolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), lty = 3) + 
    geom_point(data = subset(rare2, method == "observed" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), size = 2) + 
    ylim(0, 65) + 
    labs(x = "Number of samples", y = "Species richness") + 
    facet_grid( ~ strata, scales = "free_x") + 
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # legend.position = "none"
    )
)

### Calculate % coverage based on observed number of species and maximum extrapolated values

sel_site = "ANTARTICA"
tide_stratum = "MIDTIDE"
inter_data <- subset(rare2, method == "interpolated" & strata == tide_stratum & locality == sel_site)
extra_data <- subset(rare2, method == "extrapolated" & strata == tide_stratum & locality == sel_site)
obs_data <- subset(rare2, method == "observed" & strata == tide_stratum & locality == sel_site)

spp_cover <- obs_data$qD/max(extra_data$qD)

#####-----------------------------------------------------------------------------------------------

### Delete rows from specific localities in 'samps' data frame
p2p <- p2p[!grepl("MASSACHUSETTS", p2p$locality),] ### this allowes to generate 'stopping_plot_3'

# Run coverage optimization
samps <- do.call(rbind, lapply(unique(p2p$locality), function(i) {
  
  x <- subset(p2p, locality == as.character(i)) 
  
  do.call(rbind, lapply(unique(x$site), function(j) {
  
    x <- subset(x, site == as.character(j)) 
    
    do.call(rbind, lapply(unique(x$strata), function(k) {
    
      print(paste(i, j, k))
      
      x <- subset(x, strata == as.character(k))
      
      if(nrow(x) == 0) data.frame() else {
    
        x$presence <- 1
      
        # Cast longways
        mat <- x %>% select(quadrat, taxa, presence) %>% 
        
          pivot_wider(id_cols = quadrat, names_from = taxa, values_from = presence) 
        
        mat[is.na(mat)] <- 0
    
        data.frame(
          x[1, 1:6],
          totsamples = nrow(mat),
          minsamples = covstop(mat[, -1])
        )
        
      }
      
      } ) )
    
  } ) )
  
} ) )

samps$strata <- factor(samps$strata, levels = c("HIGHTIDE", "MIDTIDE", "LOWTIDE"))

samps$locality <- factor(samps$locality, levels = c("ANTARTICA", "PUNTAARENAS", "PUERTOMADRYN",  "CONCEPCIÃ“N", "REÃ‘ACA,VIÃ‘ADELMAR",
                                                       "ARRAIALDOCABO", "APACOSTADASALGAS", "SANTACRUZ", "FERNANDODENORONHA", "ISLAGORGONA", 
                                                      "NORTHERNMA", "BIDDEFORD","GIANTSTAIRS", "CHAMBERLAIN", "MAINE", "CENTRALMAINE"))

# Generate summary figures
samps.summary <- samps %>% group_by(locality, strata) %>% 
  
  # mutate(minsamples = minsamples / totsamples ) +
  
  summarize(mean.samples = mean(minsamples), se.samples = plotrix::std.error(minsamples), totsamples = mean(totsamples))

(stopplot <- ggplot(samps.summary, aes(x = locality, y = mean.samples, group = strata, fill = strata)) +
  geom_errorbar(aes(ymax = mean.samples + se.samples, ymin = mean.samples - se.samples), width = 0.3) +
  geom_bar(stat = "identity") +
  geom_point(aes(x = locality, y = totsamples, group = strata), shape = 23, fill = "red", size = 5) +
  facet_grid(~ strata, scale = "free_y") +
  scale_fill_manual(values = c("black", "dodgerblue3", "forestgreen")) +
  labs(x = "", y = "Minimum number of samples") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
)

ggsave("./output/Stopping plot.pdf", stopplot, device = "pdf", width = 11, height = 6, units = "in")

### This extracts maximum values of extrapolated richness from each locality
max_extra <- filter(rare2, method == "extrapolated" & locality == "ISLAGORGONA")
max(max_extra$qD)

### This calculates the difference between minimum number of samples needed for 100% coverage and actual number of samples collected

min_coverage <- samps.summary$mean.samples
collected_samps <- samps.summary$totsamples

diff_samps <- collected_samps - min_coverage
frac_samps <- diff_samps/min_coverage

### covstop vs MultiSE

covstop_data <- samps.summary[, c("locality", "strata", "mean.samples", "se.samples")]

### use this for filtering by range value (need to run multSE_SSP.R first to generate 'multse.samples')
# multse_data <- multse.samples[, c("samples", "mean", "upper", "lower", "Strata", "locality")] %>% 
#   filter(mean < 0.11) %>% filter(mean > 0.09) 

### use this for filtering by single value (need to run multSE_SSP.R first to generate 'multse.samples')
multse_data <- multse.samples[, c("samples", "mean", "upper", "lower", "Strata", "locality")] %>% 
  filter(mean < 0.11)
 
multse_min <- aggregate(multse_data$samples, by = list(multse_data$Strata, multse_data$locality), min)
# multse_mean <- aggregate(multse_data$samples, by = list(multse_data$Strata, multse_data$locality), mean)
# multse_sd <- aggregate(multse_data$samples, by = list(multse_data$Strata, multse_data$locality), sd)



  
  









