####-----------------------------------------------------------------------------------------------
####                      OPTIMIZATION OF ROCKY SHORE SAMPLING FROM P2P                           #
####-----------------------------------------------------------------------------------------------

#####-----------------------------------------------------------------------------------------------

# Load required libraries
library(iNEXT)
library(tidyverse)
library(here)

# setwd(path = "C:/Users/jslef/OneDrive - Smithsonian Institution/Documents/GitHub/asymmetric/")
# enrique change this for your local siles
# setwd("~/asymmetric")

# Load function to compute coverage-based stopping
source("R/covstop.R")

# finch::dwca_read(..., read = TRUE)

## TO DO ##
# generate rarefaction curves
# also run collapsing by locality

#####-----------------------------------------------------------------------------------------------

# Read in data

files <- list.files("./data")

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

# Complete missing replicates
with(p2p, table(site, strata, quadrat))

p2p <- p2p %>% complete(quadrat, nesting(locality, site, strata), fill = list(taxa = 0))

# Find sites where all quadrats only have a taxa code of 0
remove <- p2p %>% group_by(locality, site, strata) %>% summarize(remove = ifelse(all(taxa == "0"), TRUE, FALSE))

subset(remove, remove == TRUE)$site # no sites missing all qaudrats

# Write csv
write.csv(p2p, "./sites.csv")

# Read in site lat/longs
sites <- read.csv("./sites.csv", header = T)

#####-----------------------------------------------------------------------------------------------

# Generate interpolation/extrapolation curves for each locality (or site)

# rare <- do.call(rbind, lapply(unique(p2p$locality), function(i) {
#   
#   x <- subset(p2p, locality == i)
  
rare <- do.call(rbind, lapply(unique(p2p$site), function(j) {

  x <- subset(p2p, site == j)
    
    do.call(rbind, lapply(unique(p2p$strata), function(k) {
  
      # print(paste(i, k))
      print(paste(j, k))
      
      x <- subset(x, strata == k)
      
      if(nrow(x) == 0) data.frame() else {
        
        x$presence <- ifelse(x$taxa == 0, 0, 1)
        
        # x <- x %>% group_by(locality, site, strata, quadrat, taxa) %>% summarize(presence = sum(presence))
          
        # Cast longways
        mat <- x %>% select(site, quadrat, taxa, presence) %>% 
          
          pivot_wider(id_cols = c(site, quadrat), names_from = taxa, values_from = presence, 
                      values_fn = mean) 
        
        mat[is.na(mat)] <- 0
        
        # remove taxa placeholder
        mat <- mat[, colnames(mat) != "0"]
        
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
          locality = unique(x$locality),
          site = j,
          strata = k,
          ret[, c(1:2, 4)]
        )
        
      }
      
      } ) )
    
    } ) )
  
  # } ) )

rare$strata <- factor(rare$strata, levels = c("HIGHTIDE", "MIDTIDE", "LOWTIDE"))

# # Plot results
# (rareplot <- ggplot() +
#   geom_line(data = subset(rare, method == "interpolated"), aes(x = t, y = qD, group = paste(site, strata), col = strata)) +
#   geom_line(data = subset(rare, method == "extrapolated"), aes(x = t, y = qD, group = paste(site, strata), col = strata), lty = 3) +
#   geom_point(data = subset(rare, method == "observed"), aes(x = t, y = qD, group = paste(site, strata), col = strata), size = 2) +
#   scale_color_manual(values = c("black", "dodgerblue3", "forestgreen")) +
#   labs(x = "Number of samples", y = "Species richness") +
#   facet_grid( ~ strata, scales = "free_x") +
#   theme_bw(base_size = 12) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.position = "none"
#   )
# )

# ggsave("./output/Rarefaction plot.pdf", rareplot, device = "pdf", width = 10, height = 5, units = "in")

rare_summ <- rare %>% group_by(method, locality, strata, t) %>% 
  
  mutate(locality = as.factor(as.character(locality))) %>% summarize(qD = mean(qD)) 

# select the order of localities in the legend 
rare_summ$locality <- factor(rare_summ$locality, levels = c("ANTARTICA", "PUNTAARENAS", "PUERTOMADRYN",  "CONCEPCIÃ“N", "REÃ‘ACA,VIÃ‘ADELMAR",
                                                    "ARRAIALDOCABO", "APACOSTADASALGAS", "SANTACRUZ", "FERNANDODENORONHA", "ISLAGORGONA", 
                                                    "NORTHERNMA", "BIDDEFORD","GIANTSTAIRS", "CHAMBERLAIN", "GRINDSTONE", "CENTRALMAINE"))

# here is where you can change the locality names
levels(rare_summ$locality) <- c("Antarctica", "Punta Arenas", "Puerto Madryn", "Concepción", 
                                "Montemar", "Arraial do Cabo", "Costa das Algas", "Santa Cruz", "Fernando de Noronha", 
                                "Isla Gorgona", "Massachusetts", "Biddeford", "Giant Stairs", "Chamberlain", 
                                "Grindstone", "N. Maine") 

# Plot results with curves colored according to locality  
(rareplot_1 <- ggplot() +
  geom_line(data = subset(rare_summ, method == "interpolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), lty = 1, size = 1) + 
  geom_line(data = subset(rare_summ, method == "extrapolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), lty = 3, size = 1) + 
  geom_point(data = subset(rare_summ, method == "observed" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), size = 2) + 
  xlim(0, 25) +
  ylim(0, 40) + 
  labs(x = "Number of samples", y = "Species richness") + 
  facet_grid( ~ strata, scales = "free_x") + 
  theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=20),
      text = element_text(size = 20)
  )
)

# Plot results with curves colored according to site  
rare_summ_2 <- rare %>% group_by(method, site, strata, t) %>% 
  mutate(site = as.factor(as.character(site))) %>% summarize(qD = mean(qD)) 

levels(rare_summ_2$site) <- c() # here is where you can change the site names

# Plot results with curves colored according to site  
(rareplot_1b <- ggplot() +
    geom_line(data = subset(rare_summ_2, method == "interpolated" & strata == strata), aes(x = t, y = qD, group = paste(site, strata), col = site)) + 
    geom_line(data = subset(rare_summ_2, method == "extrapolated" & strata == strata), aes(x = t, y = qD, group = paste(site, strata), col = site), lty = 3) + 
    geom_point(data = subset(rare_summ_2, method == "observed" & strata == strata), aes(x = t, y = qD, group = paste(site, strata), col = site), size = 2) + 
    ylim(0, 40) + 
    labs(x = "Number of samples", y = "Species richness") + 
    facet_grid( ~ strata, scales = "free_x") + 
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # legend.position = "none"
    )
)

# ggsave("./output/Rarefaction plot 2.pdf", rareplot_1, device = "pdf", width = 10, height = 5, units = "in")  

### This plots individual sites

sel_locality = c("Santa Cruz")

(rareplot_2 <- ggplot() +
    geom_line(data = subset(rare_summ, method == "interpolated" & strata == strata & locality == sel_locality), aes(x = t, y = qD, group = paste(strata))) +
    geom_line(data = subset(rare_summ, method == "extrapolated" & strata == strata & locality == sel_locality), aes(x = t, y = qD, group = paste(strata)), lty = 3) +
    geom_point(data = subset(rare_summ, method == "observed" & strata == strata & locality == sel_locality), aes(x = t, y = qD, group = paste(strata)), size = 2) +
    ylim(0, 40) +
    xlim(0, 25) +
    labs(x = "Number of samples", y = "Species richness") +
    facet_grid( ~ strata, scales = "free_x") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # legend.position = "none"
    )
)
# 
# ### Delete rows from specific localities in 'rare' data frame
# rare2 <- rare[!grepl("MASSACHUSETTS", rare$locality),] ### this allows to generate 'rare_plot_2'
# 
# # Plot results with curves colored according to locality  
# (rareplot_2 <- ggplot() +
#     geom_line(data = subset(rare2, method == "interpolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality)) + 
#     geom_line(data = subset(rare2, method == "extrapolated" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), lty = 3) + 
#     geom_point(data = subset(rare2, method == "observed" & strata == strata), aes(x = t, y = qD, group = paste(locality, strata), col = locality), size = 2) + 
#     ylim(0, 65) + 
#     labs(x = "Number of samples", y = "Species richness") + 
#     facet_grid( ~ strata, scales = "free_x") + 
#     theme_bw(base_size = 14) +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       # legend.position = "none"
#     )
# )

### Calculate % coverage based on observed number of species and maximum extrapolated values (for coverage_obs_vs_extrapolated_v2.csv)

sel_loc = "NORTHERNMA"
inter_data_h <- subset(rare_summ, method == "interpolated" & strata == "HIGHTIDE" & locality == sel_loc)
extra_data_h  <- subset(rare_summ, method == "extrapolated" & strata == "HIGHTIDE" & locality == sel_loc)
obs_data_h  <- subset(rare_summ, method == "observed" & strata == "HIGHTIDE" & locality == sel_loc)
spp_cover_h  <- obs_data_h$qD/max(extra_data_h$qD)

inter_data_m <- subset(rare_summ, method == "interpolated" & strata == "MIDTIDE" & locality == sel_loc)
extra_data_m  <- subset(rare_summ, method == "extrapolated" & strata == "MIDTIDE" & locality == sel_loc)
obs_data_m  <- subset(rare_summ, method == "observed" & strata == "MIDTIDE" & locality == sel_loc)
spp_cover_m  <- obs_data_m$qD/max(extra_data_m$qD)

inter_data_l <- subset(rare_summ, method == "interpolated" & strata == "LOWTIDE" & locality == sel_loc)
extra_data_l  <- subset(rare_summ, method == "extrapolated" & strata == "LOWTIDE" & locality == sel_loc)
obs_data_l  <- subset(rare_summ, method == "observed" & strata == "LOWTIDE" & locality == sel_loc)
spp_cover_l  <- obs_data_l$qD/max(extra_data_l$qD)

spp_cover_h
spp_cover_m
spp_cover_l

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
    
        x$presence <- ifelse(x$taxa == 0, 0, 1)
        
        # x <- x %>% group_by(locality, site, strata, quadrat, taxa) %>% summarize(presence = sum(presence))
        
        # Cast longways
        mat <- x %>% select(site, quadrat, taxa, presence) %>% 
          
          pivot_wider(id_cols = c(site, quadrat), names_from = taxa, values_from = presence, 
                      values_fn = mean) 
        
        mat[is.na(mat)] <- 0
        
        # remove taxa placeholder
        mat <- mat[, colnames(mat) != "0"]
    
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
                                                      "NORTHERNMA", "BIDDEFORD","GIANTSTAIRS", "CHAMBERLAIN", "GRINDSTONE", "CENTRALMAINE"))

# Generate summary figures
samps.summary <- samps %>% group_by(locality, strata) %>% 
  
  # mutate(minsamples = minsamples / totsamples ) +
  
  summarize(mean.samples = mean(minsamples), se.samples = plotrix::std.error(minsamples), totsamples = mean(totsamples))

# # this replaces the number of samples collected in specific localities 
samps.summary$totsamples[which(samps.summary$locality == "NORTHERNMA")] <- 5
samps.summary$totsamples[which(samps.summary$locality == "BIDDEFORD")] <- 5
samps.summary$totsamples[which(samps.summary$locality == "GIANTSTAIRS")] <- 5
samps.summary$totsamples[which(samps.summary$locality == "CHAMBERLAIN")] <- 5
samps.summary$totsamples[which(samps.summary$locality == "GRINDSTONE")] <- 5
samps.summary$totsamples[which(samps.summary$locality == "CENTRALMAINE")] <- 5

(stopplot <- ggplot(samps.summary, aes(x = locality, y = mean.samples, group = strata, fill = strata)) +
  geom_hline(yintercept = 10, col = "grey50", linetype = "dashed") +
  geom_errorbar(aes(ymax = mean.samples + se.samples, ymin = mean.samples - se.samples), width = 0.3) +
  geom_bar(stat = "identity") +
  geom_point(aes(x = locality, y = totsamples, group = strata), shape = 23, fill = "red", size = 3) +
  facet_grid(~ strata, scale = "free_y") +
  scale_fill_manual(values = c("black", "dodgerblue3", "forestgreen")) +
  lims(y = c(0, 30)) +
  labs(x = "", y = "Minimum number of samples") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
)

#ggsave("./output/Stopping plot.pdf", stopplot, device = "pdf", width = 11, height = 6, units = "in")

### This extracts maximum values of extrapolated and observed richness from each locality (values are aggregated in lat_vs_ssp_p2p_v3.csv)
sel_locality = "CENTRALMAINE"
max_extra_h <- filter(rare_summ, method == "extrapolated" & locality == sel_locality & strata == "HIGHTIDE")
obs_h <- filter(rare_summ, method == "observed" & locality == sel_locality & strata == "HIGHTIDE")
max_extra_m <- filter(rare_summ, method == "extrapolated" & locality == sel_locality & strata == "MIDTIDE")
obs_m <- filter(rare_summ, method == "observed" & locality == sel_locality & strata == "MIDTIDE")
max_extra_l <- filter(rare_summ, method == "extrapolated" & locality == sel_locality & strata == "LOWTIDE")
obs_l <- filter(rare_summ, method == "observed" & locality == sel_locality & strata == "LOWTIDE")

max(max_extra_h$qD)
obs_h$qD
max(max_extra_m$qD)
obs_m$qD
max(max_extra_l$qD)
obs_l$qD

### Calculate min number of samples needed to cover 90% of maximum extrapolated richness (to generate Table S2)
sel_site = "PUMPHOUSE"
extrapol_h <- subset(samps, site == sel_site & strata == "HIGHTIDE") 
extrapol_m <- subset(samps, site == sel_site & strata == "MIDTIDE")
extrapol_l <- subset(samps, site == sel_site & strata == "LOWTIDE")

min_sample_h <- extrapol_h$minsamples
min_sample_m <- extrapol_m$minsamples
min_sample_l <- extrapol_l$minsamples
buff_h <- min_sample_h*0.2
buff_m <- min_sample_m*0.2
buff_l <- min_sample_l*0.2

min_sample_h
buff_h
min_sample_m
buff_m
min_sample_l
buff_l

### This calculates the difference between minimum number of samples needed for 100% coverage and actual number of samples collected

min_coverage <- samps.summary$mean.samples
collected_samps <- samps.summary$totsamples

diff_samps <- collected_samps - min_coverage
frac_samps <- diff_samps/min_coverage

### covstop vs MultiSE

covstop_data <- samps.summary[, c("locality", "strata", "mean.samples", "se.samples")]

### use this for filtering by range value (need to run multSE_SSP.R first [lines 357-387] to generate 'multse.samples')
# multse_data <- multse.samples[, c("samples", "mean", "upper", "lower", "Strata", "locality")] %>% 
#   filter(mean < 0.11) %>% filter(mean > 0.09) 

### use this for filtering by single value (need to run multSE_SSP.R first to generate 'multse.samples')
multse_data <- multse.samples[, c("samples", "mean", "upper", "lower", "Strata", "locality")] %>% 
  filter(mean <= 0.1)
 
multse_min <- aggregate(multse_data$samples, by = list(multse_data$Strata, multse_data$locality), min)
# multse_mean <- aggregate(multse_data$samples, by = list(multse_data$Strata, multse_data$locality), mean)
# multse_sd <- aggregate(multse_data$samples, by = list(multse_data$Strata, multse_data$locality), sd)

### plot regressions between covstop and multse

csv_table <- 'multse_v_covstop_csv_table_v2.csv'
arrays <-read.csv(csv_table, row.names = 1)

multse_var <- arrays[, 3]
covstop_var <- arrays[, 6]

theme_set(theme_bw())
ggplot(arrays, aes(x = multse_var, y = covstop_var)) + 
  geom_point(size = 3) +
  stat_smooth(method = "lm", col = "red") +
  geom_abline(slope=1, intercept = 0, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 18, by = 2)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5)) +
  coord_cartesian(ylim=c(2,20), xlim=c(2,17)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(text = element_text(size = 20)) +
  xlab("MultSE minimum sample") +
  ylab("Covstop minimum sample")
  

### plot latitude versus observed and extrapolated ssp richness

ssp_lat <-read.csv('lat_vs_spp_p2p_v3.csv') %>% 
  select(latitude, qD_max, obs_max) %>%
  gather(key = "variable", value = "value", -latitude)

theme_set(theme_bw())
ggplot(ssp_lat, aes(x = latitude, y = value)) + 
  geom_line(aes(color = variable)) + 
  scale_color_manual(values = c("darkred", "steelblue"), labels = c("Extrapolated", "Observed")) +
  xlim(-65, 50) +
  ylim(0, 40) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(text = element_text(size = 20)) +
  xlab("Latitude") +
  ylab("Species richness")

### This is the plot used in the manuscript (the observed values are the maximum at low, mid and high tide, which in turn are mean values for each stratum)
ssp_lat <-read.csv('lat_vs_spp_p2p_v3.csv')
ssp_lat2 <- ssp_lat[-c(8, 9), ] # this removes F. de Noronha and Sta Cruz

theme_set(theme_bw())
ggplot() + 
  geom_line(data = ssp_lat, aes(x = latitude, y = max_obs_all, color = 'Observed'), linetype = 1, size = 1.5) + 
  geom_line(data = ssp_lat, aes(x = latitude, y = max_estrapol_all, color = 'Extrapolated'), linetype = 1, size = 1.5) + 
  geom_line(data = ssp_lat2, aes(x = latitude, y = max_obs_all, color = 'Observed'), linetype = 3, size = 1.5) + 
  geom_line(data = ssp_lat2, aes(x = latitude, y = max_estrapol_all, color = 'Extrapolated'), linetype = 3, size = 1.5) + 
  xlim(-65, 50) +
  ylim(0, 40) +
  scale_x_continuous(breaks = seq(-70, 50, by = 10)) +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(text = element_text(size = 20)) +
  xlab("Latitude") +
  ylab("Species richness")

# end

