#### This script maps number of records per locality and site 
#### #### from: https://www.r-graph-gallery.com/19-map-leafletr.html
#### Enrique Montes
#### October 22, 2019

# Library
library(leaflet)

setwd("G:/My Drive/MBON/GLOBAL_MBON/P2P/P2P_Project/papers/manuscripts/network_sampling_optimization")

rm(list=ls()) ## clear variables

data <- read.csv("pole2pole_data_summary.csv", header = TRUE)

long <- data$decimalLongitude
lat <- data$decimalLatitude

# Create a color palette with handmade bins.
# for number of records
# mybins <- seq(10, 250, by=60)
# mypalette <- colorBin( palette="YlOrBr", domain=data$total.records, na.color="transparent", bins=mybins)
# 
# # Prepare the text for the tooltip:
# mytext <- paste(
#   "Locality: ", data$Locality, "<br/>", 
#   "Site: ", data$Site, "<br/>", 
#   "Number of records: ", data$total.records, sep="") %>%
#   lapply(htmltools::HTML)
# 
# # Prepare the text for the tooltip:
# mytext <- paste(
#   "Locality: ", data$Locality, "<br/>", 
#   "Site: ", data$Site, "<br/>", 
#   "Number of records: ", data$total.records, sep="") %>%
#   lapply(htmltools::HTML)
# 
# # Final Map
# m <- leaflet(data) %>% 
#   addTiles()  %>% 
#   setView( lat=0, lng=-50 , zoom=3) %>%
#   addProviderTiles("Esri.WorldImagery") %>%
#   addCircleMarkers(~long, ~lat, 
#                    fillColor = ~mypalette(total.records), fillOpacity = 0.7, color="white", radius=12, stroke=FALSE,
#                    label = mytext,
#                    labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), textsize = "13px", direction = "auto")
#   ) %>%
#   addLegend( pal=mypalette, values=~total.records, opacity=0.9, title = "Number of records", position = "bottomright" )
# 
# m 

# for number of taxa (richness)
mybins <- seq(0, 60, by=10)
mypalette <- colorBin( palette="YlOrBr", domain=data$total.spp, na.color="transparent", bins=mybins)

# Prepare the text for the tooltip:
mytext <- paste(
  "Locality: ", data$Locality, "<br/>", 
  "Site: ", data$Site, "<br/>", 
  "Number of species: ", data$total.spp, sep="") %>%
  lapply(htmltools::HTML)

# Final Map
m <- leaflet(data) %>% 
  addTiles()  %>% 
  setView( lat=0, lng=-50 , zoom=3) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~long, ~lat, 
                   fillColor = ~mypalette(total.spp), fillOpacity = 0.7, color="white", radius=12, stroke=FALSE,
                   label = mytext,
                   labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), textsize = "13px", direction = "auto")
  ) %>%
  addLegend( pal=mypalette, values=~total.spp, opacity=0.9, title = "Number of species", position = "bottomright" )

m 