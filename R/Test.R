library(tidyverse)

# finch::dwca_read(..., read = TRUE)

files <- list.files("C:/Users/jslef/OneDrive - Smithsonian Institution/Documents/GitHub/asymmetric/data")

out <- do.call(rbind, lapply(files, function(l) {
  
  # read in data
  dat <- read.csv(paste0("C:/Users/jslef/OneDrive - Smithsonian Institution/Documents/GitHub/asymmetric/data", l))
  
  # split metadata
  meta <- do.call(rbind.data.frame, strsplit(as.character(dat$eventID), "\\_"))
  
  names(meta) <- c("country", "locality", "site", "date", "strata", "quadrat")
  
  taxa <- do.call(rbind.data.frame, strsplit(as.character(dat$scientificNameID), "\\:"))[, 5]
  
  # bind back in
  data <- cbind.data.frame(meta, taxa)
  
  # split by site
  do.call(rbind, lapply(unique(data$site), function(i) {
    
    x <- subset(data, site == as.character(i)) 
    
    # split by strata
    do.call(rbind, lapply(unique(x$strata), function(j) {
      
      comm <- subset(x, strata == as.character(j))
      
      comm$presence <- 1
      
      # remove duplicate rows
      comm <- comm %>% group_by(country, locality, site, date) %>% distinct() %>% ungroup()
      
      # cast longways
      mat <- comm %>% select(quadrat, taxa, presence) %>% 
        
        pivot_wider(id_cols = quadrat, names_from = taxa, values_from = presence) 
        
      mat[is.na(mat)] <- 0
  
      data.frame(
        comm[1, 1:5],
        totsamples = nrow(mat),
        minsamples = covstop(mat[, -1])
      )
      
    } ) )
    
  } ) )
  
} ) )

out.summary <- out %>% group_by(strata) %>% summarize(mean.samples = mean(minsamples), se.samples = plotrix::std.error(minsamples))

ggplot(out.summary, aes(x = strata, y = mean.samples)) +
  geom_errorbar(aes(ymax = mean.samples + se.samples, ymin = mean.samples - se.samples)) +
  geom_bar(stat = "identity")
  