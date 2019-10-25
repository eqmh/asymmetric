#### This script extracts information from DwC-A files: occurrence and event core.
#### Outputs: number of records and species at high, mid and low tide; number of quadrats per stratum
#### Enrique Montes
#### October 21, 2019


setwd("G:/My Drive/MBON/GLOBAL_MBON/P2P/P2P_Project/Biodiv_surveys/00-DataCurated/IPTFiles/rocky")

rm(list=ls()) ## clear variables

data <- read.csv("ARGENTINA-PUERTO MADRYN_ipt_occurrence.csv", header = TRUE)
event <- read.csv("ARGENTINA-PUERTO MADRYN_ipt_event.csv", header = TRUE)

sel_str <- "PUNTALOMA" ## select site

df <- data$occurrenceID[grep(sel_str, data$occurrenceID)] ## extract indices for rows with selected site

df2 <- as.data.frame(df) ## subset data according to df

t_taxa <- nrow(df2)  ## total number of species records for selected site

strat_high <- df2$df[grep("HIGHTIDE", df2$df)] 
strat_mid <- df2$df[grep("MIDTIDE", df2$df)]
strat_low <- df2$df[grep("LOWTIDE", df2$df)]

df_high <- as.data.frame(strat_high)
df_mid <- as.data.frame(strat_mid)
df_low <- as.data.frame(strat_low)

high_rec <- nrow(data[strat_high,])
mid_rec <- nrow(data[strat_mid,])
low_rec <- nrow(data[strat_low,])

tb <- c(t_taxa, high_rec, mid_rec, low_rec)
tb

num_quad <- event$eventID[grep(sel_str, event$eventID)] 
num_quad <- as.data.frame(num_quad)

num_quad_high <- num_quad$num_quad[grep("HIGHTIDE_R", num_quad$num_quad)] 
num_quad_mid <- num_quad$num_quad[grep("MIDTIDE_R", num_quad$num_quad)]
num_quad_low <- num_quad$num_quad[grep("LOWTIDE_R", num_quad$num_quad)]

high_quad <- nrow(event[num_quad_high,])
mid_quad <- nrow(event[num_quad_mid,])
low_quad <- nrow(event[num_quad_low,])

quad_count <- c(high_quad, mid_quad, low_quad)
quad_count

## Count number of unique species 
sub_df_high <- data[strat_high,]
sppList_high <- sub_df_high$scientificName
spp_high <- nrow(as.data.frame(unique(sppList_high)))

sub_df_mid <- data[strat_mid,]
sppList_mid <- sub_df_mid$scientificName
spp_mid <- nrow(as.data.frame(unique(sppList_mid)))

sub_df_low <- data[strat_low,]
sppList_low <- sub_df_low$scientificName
spp_low <- nrow(as.data.frame(unique(sppList_low)))

sub_spp_all <- data[df,]
sppList_all <- sub_spp_all$scientificName
spp_all <- nrow(as.data.frame(unique(sppList_all)))

f_list <- c(spp_high, spp_mid, spp_low, spp_all)
f_list
