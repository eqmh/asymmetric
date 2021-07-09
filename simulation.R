####-----------------------------------------------------------------------------------------------
####                      SIMULATION OF ROCKY SHORE SAMPLING EFFORT FROM P2P                           #
####-----------------------------------------------------------------------------------------------

#####-----------------------------------------------------------------------------------------------

# Load required libraries
library(iNEXT)
library(tidyverse)
library(ggplot2)

# setwd(path = "C:/Users/jslef/OneDrive - Smithsonian Institution/Documents/GitHub/asymmetric/")
setwd("~/asymmetric")

# Load function to compute coverage-based stopping
source("R/covstop.R")


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

low.tide <- filter(p2p, country == "COLOMBIA", strata == "LOWTIDE")
mid.tide <- filter(p2p, country == "COLOMBIA", strata == "MIDTIDE")
high.tide <- filter(p2p, country == "COLOMBIA", strata == "HIGHTIDE")

# extracts quadrat number from 'quadrat' strings and turns it into numeric value
# low tide
q_num.low <- as.data.frame(as.numeric(substr(low.tide$quadrat, 3, 4)))
colnames(q_num.low) <- "q_num.low"
low.tide <- bind_cols(low.tide, q_num.low)

# mid tide
q_num.mid <- as.data.frame(as.numeric(substr(mid.tide$quadrat, 3, 4)))
colnames(q_num.mid) <- "q_num.mid"
mid.tide <- bind_cols(mid.tide, q_num.mid)

# high tide
q_num.high<- as.data.frame(as.numeric(substr(high.tide$quadrat, 3, 4)))
colnames(q_num.high) <- "q_num.high"
high.tide <- bind_cols(high.tide, q_num.high)


# simulate taxa richness per number of quadrats collected (with replacement)
iterations <- 100
df <- matrix(ncol = 1, nrow = iterations)  # creates an empty data frame
for ( i in 1:100){

  q_sel.mid <- sample(1:10, i, replace = T) # generate random numbers from 1 to 10 with replacement
  mid.tide.sub <- mid.tide[mid.tide$q_num.mid %in% q_sel.mid, ]  # selects quadrats based on q_sel list
  num.unique.mid <- apply(mid.tide.sub, 2, function(x) length(unique(x)))  # pulls unique values from each column
  sim.vals.mid <- unname(num.unique.mid[8])  # pulls unique taxa values
  
  df[i,] <- cbind(sim.vals.mid)  # concatenates rows

}

df <- as.data.frame(df)





