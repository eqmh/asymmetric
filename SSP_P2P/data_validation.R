library(tidyr)
library(stringr)
library(dplyr)

files <- list.files("~/asymmetric/data")

one <- function(x){
    pa <- sum(x)/sum(x)
    return(pa)
  }

for ( i in 1:length(files)){
  
  dat<-read.csv(files[i])
  dat$occurrenceStatus <- 1

  
  
for (i in 1:nrow(dat))  

  df <- pivot_wider(data = dat, id_cols = eventID, 
                    names_from = scientificName, 
                    values_from = occurrenceStatus,
                    values_fill = list(occurrenceStatus = 0),
                    values_fn = list(occurrenceStatus = one))
  
  info <- df$eventID %>% 
          str_split("_", simplify = TRUE)

  info <- as.data.frame(info)
  colnames(info) <- c("Country", "Locality",	"Site", "Date", "Strata", "sample")
  info$Strata <- str_to_lower(info$Strata)
  out <- bind_cols(info[,c(-4, -6)], df[,-1])
  out$Site <- factor(out$Site)
    if (length(levels(out$Site)) == 1) {
      dt <- "C:/Users/Edlin/Documents/asymmetric/SSP_P2P/single"
      loc <- paste0(out$Country[1], "_", out$Locality[1], ".csv")
      path <- paste (dt, loc, sep = "/")
      write.csv(out, file = path)
    } else {
      dt <- "C:/Users/Edlin/Documents/asymmetric/SSP_P2P/multiple"
      loc <- paste0(out$Country[1], "_", out$Locality[1], ".csv")
      path <- paste (dt, loc, sep = "/")
      write.csv(out, file = path)
    }
     
  }
   





