####-----------------------------------------------------------------------------------------------
####                      SIMULATION OF ROCKY SHORE SAMPLING EFFORT FROM P2P                           #
####-----------------------------------------------------------------------------------------------

#####-----------------------------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(ggplot2)
library(reshape2)

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
iterations <- 1000
q_rep <- 50
df <- matrix(ncol = 3, nrow = q_rep)  # creates an empty data frame
df.low <- matrix(ncol = iterations, nrow = q_rep) 
df.mid <- matrix(ncol = iterations, nrow = q_rep) 
df.high <- matrix(ncol = iterations, nrow = q_rep) 

for (j in 1:iterations){
    for ( i in 1:q_rep){
      
      # low tide
      q_sel.low<- sample(1:10, i, replace = T) # generate random numbers from 1 to 10 with replacement
      low.tide.sub <- low.tide[low.tide$q_num.low %in% q_sel.low, ]  # selects quadrats based on q_sel list
      num.unique.low <- apply(low.tide.sub, 2, function(x) length(unique(x)))  # pulls unique values from each column
      sim.vals.low <- unname(num.unique.low[8])  # pulls unique taxa values
      
      # Mid tide
      q_sel.mid <- sample(1:10, i, replace = T) # generate random numbers from 1 to 10 with replacement
      mid.tide.sub <- mid.tide[mid.tide$q_num.mid %in% q_sel.mid, ]  # selects quadrats based on q_sel list
      num.unique.mid <- apply(mid.tide.sub, 2, function(x) length(unique(x)))  # pulls unique values from each column
      sim.vals.mid <- unname(num.unique.mid[8])  # pulls unique taxa values
      
      # High tide
      q_sel.high<- sample(1:10, i, replace = T) # generate random numbers from 1 to 10 with replacement
      high.tide.sub <- high.tide[high.tide$q_num.high %in% q_sel.high, ]  # selects quadrats based on q_sel list
      num.unique.high <- apply(high.tide.sub, 2, function(x) length(unique(x)))  # pulls unique values from each column
      sim.vals.high <- unname(num.unique.high[8])  # pulls unique taxa values
      
      df[i,] <- cbind(sim.vals.low, sim.vals.mid, sim.vals.high)  # concatenates rows
    }
  df.low[, j] <- cbind(df[, 1])
  df.mid[, j] <- cbind(df[, 2])
  df.high[, j] <- cbind(df[, 3])
  df <- matrix(ncol = 3, nrow = q_rep) 
}

df.low <- as.data.frame(df.low)
df.mid<- as.data.frame(df.mid)
df.high <- as.data.frame(df.high)

mean.low <- rowMeans(df.low)
mean.mid <- rowMeans(df.mid)
mean.high <- rowMeans(df.high)
sd.low <- apply(df.low, 1, sd)
sd.mid <- apply(df.mid, 1, sd)
sd.high <- apply(df.high, 1, sd)
summ.vals <- as.data.frame(cbind(1:q_rep, mean.low, mean.mid, mean.high, 
                                 mean.low - sd.low, 
                                 mean.low + sd.low,
                                 mean.mid - sd.mid,
                                 mean.mid + sd.mid,
                                 mean.high - sd.high,
                                 mean.high + sd.high))
colnames(summ.vals) <- c("samples", "LOW", "MID", "HIGH", "l.low", "u.low", "l.mid", "u.mid", "l.high", "u.high")



p <- ggplot() + 
  geom_pointrange(summ.vals, mapping = aes(1:q_rep, mean.low, ymin = l.low, ymax = u.low), color = "blue") +
  geom_pointrange(summ.vals, mapping = aes(1:q_rep, mean.mid, ymin = l.mid, ymax = u.mid), color = "red") +
  geom_pointrange(summ.vals, mapping = aes(1:q_rep, mean.high, ymin = l.high, ymax = u.high), color = "green") +
  labs(x = "Number of quadrats", y = "Species richness", title = "Gorgona Island") +
  scale_color_manual(labels = c("low", "mid", "high"), values = c("blue", "red", "green")) +
  theme_bw() + 
  theme(legend.position= "right") 
p

mdf.avg <- melt(summ.vals[, 1:4], id.vars = c("samples"))
q <- ggplot() +geom_point(data = mdf.avg, mapping = aes(x = samples, y = value, colour = variable))



