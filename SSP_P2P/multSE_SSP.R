####-----------------------------------------------------------------------------------------------
####        OPTIMIZATION OF ROCKY SHORE SAMPLING FROM P2P   / multSE using SSP                   #
####-----------------------------------------------------------------------------------------------

#####-----------------------------------------------------------------------------------------------

# Load required libraries
library(SSP)
library(sampling)
library(vegan)
library(dplyr)
library(ggplot2)

############ Localities with more than one site, all tide levels##########
#Set as working directory the folder "multiple"
setwd("~/asymmetric/SSP_P2P/multiple")

files <- list.files()

n = 100 # 100 
cases = 20 # 100 
sites = 30 # 30 
p.n = 30 # 30 
p.s = 20 # 20 
k = 10 # 100 

for ( i in 1:length(files)){
  
  dat<-read.csv(files[i])
  dat[is.na(dat)]<-0
  
  ##split the data by strata
  dat.h<- dat[dat$Strata == "hightide",]
  dat.m<- dat[dat$Strata == "midtide",]
  dat.l<- dat[dat$Strata == "lowtide",]
  
  ##adjust the data for the format required by SSP
  
  dat.h<-dat.h[,c(3,5:length(dat.h))]
  dat.m<-dat.m[,c(3,5:length(dat.m))]
  dat.l<-dat.l[,c(3,5:length(dat.l))]
  
  ##SSP for hightide
  
  #parameters for simulation
  
  h.par<-assempar(data = dat.h, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  h.sim<-simdata(Par = h.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set units  and 30 sites
  h.samp<-sampsd(dat.sim = h.sim, Par = h.par, transformation = "fourth root",
                 method = "bray", multi.site = TRUE, n = n, p.n = p.n,
                 sites = sites, p.s = p.s, k = k)
  
  # average of multse for each potential sampling design
  h.sum<-summary_ssp(h.samp)
  h.sum$Strata<-rep("hightide", nrow(h.sum))
  
  
  ##SSP for midtide
  #parameters for simulation
  
  m.par<-assempar(data = dat.m, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  m.sim<-simdata(Par = m.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set 
  m.samp<-sampsd(dat.sim = m.sim, Par = m.par, transformation = "fourth root",
                 method = "bray", multi.site = TRUE, n = n, p.n = p.n,
                 sites = sites, p.s = p.s, k = k)
  
  
  # average of multse for each potential sampling design
  m.sum<-summary_ssp(m.samp)
  m.sum$Strata<-rep("midtide", nrow(m.sum))
  
  ##SSP for lowtide
  
  #parameters for simulation
  
  l.par<-assempar(data = dat.l, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  l.sim<-simdata(Par = l.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set 
  l.samp<-sampsd(dat.sim = l.sim, Par = l.par, transformation = "fourth root",
                 method = "bray", multi.site = TRUE, n = n, p.n = p.n,
                 sites = sites, p.s = p.s, k = k)
  
  # average of multse for each potential sampling design
  l.sum<-summary_ssp(l.samp)
  l.sum$Strata<-rep("lowtide", nrow(l.sum))
  
  #Combine multse from simulations in a single data frame
  
  
  result<-rbind(h.sum, m.sum, l.sum)
  result$locality<-rep(as.character(dat$Locality[1]), nrow(result))
  write.csv(result, file = paste(result$locality[1], ".csv"))
  print(paste("loc", i))
  
}
#cut and paste all the result files in the folder "SSP.output"

############# Localities with one site ###################################
#Set as working directory the folder "single"
setwd("~/asymmetric/SSP_P2P/single")

files <- list.files()

n = 100 # 100 
cases = 20 # 100 
sites = 30 # 30 
p.n = 30 # 30 
p.s = 20 # 20 
k = 10 # 100

for ( i in 1:length(files)){
  
  dat<-read.csv(files[i])
  dat[is.na(dat)]<-0
  
  ##split the data by strata
  dat.h<- dat[dat$Strata == "hightide",]
  dat.m<- dat[dat$Strata == "midtide",]
  dat.l<- dat[dat$Strata == "lowtide",]
  
  ##adjust the data for the format required by SSP
  
  dat.h<-dat.h[,c(3,5:length(dat.h))]
  dat.m<-dat.m[,c(3,5:length(dat.m))]
  dat.l<-dat.l[,c(3,5:length(dat.l))]
  
  ##SSP for hightide
  
  #parameters for simulation
  
  h.par<-assempar(data = dat.h, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  h.sim<-simdata(Par = h.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set 
  h.samp<-sampsd(dat.sim = h.sim, Par = h.par, transformation = "fourth root",
                 method = "bray", multi.site = FALSE, n = n, p.n = p.n, k = k)
  
  # average of multse for each potential sampling design
  h.sum<-summary_ssp(h.samp, multi.site = FALSE)
  h.sum$Strata<-rep("hightide", nrow(h.sum))
  
  
  ##SSP for midtide
  #parameters for simulation
  
  m.par<-assempar(data = dat.m, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  m.sim<-simdata(Par = m.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set 
  m.samp<-sampsd(dat.sim = m.sim, Par = m.par, transformation = "fourth root",
                 method = "bray", multi.site = FALSE, n = n, p.n = p.n, k = k)
  
  
  # average of multse for each potential sampling design
  m.sum<-summary_ssp(m.samp, multi.site = FALSE)
  m.sum$Strata<-rep("midtide", nrow(m.sum))
  
  ##SSP for lowtide
  
  #parameters for simulation
  
  l.par<-assempar(data = dat.l, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  l.sim<-simdata(Par = l.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set 
  l.samp<-sampsd(dat.sim = l.sim, Par = l.par, transformation = "fourth root",
                 method = "bray", multi.site = FALSE, n = n, p.n = p.n, k = k)
  
  # average of multse for each potential sampling design
  l.sum<-summary_ssp(l.samp, multi.site = FALSE)
  l.sum$Strata<-rep("lowtide", nrow(l.sum))
  
  #Combine multse from simulations in a single data frame
  
  
  result<-rbind(h.sum, m.sum, l.sum)
  result$locality<-rep(as.character(dat$Locality[1]), nrow(result))
  result$sv<-rep("samples", nrow(result))
  write.csv(result, file = paste(result$locality[1], ".csv"))
  print(paste("loc", i))
  
}
#cut and paste all the result files in the folder "SSP.output"

############# US sites####################################################
#Set as working directory the folder "USA"
setwd("~/asymmetric/SSP_P2P/USA")

files <- list.files()
 
n = 100 # 100 
cases = 20 # 100 
sites = 30 # 30 
p.n = 30 # 30 
p.s = 20 # 20 
k = 10 # 100


for ( i in 1:length(files)){
  
  dat<-read.csv(files[i])
  dat[is.na(dat)]<-0
  
  ##split the data by strata
  dat.h<- dat[dat$Strata == "hightide",]
  dat.m<- dat[dat$Strata == "midtide",]
  #dat.l<- dat[dat$Strata == "lowtide",]
  
  ##adjust the data for the format required by SSP
  
  dat.h<-dat.h[,c(3,5:length(dat.h))]
  dat.m<-dat.m[,c(3,5:length(dat.m))]
  #dat.l<-dat.l[,c(3,5:length(dat.l))]
  
  ##SSP for hightide
  
  #parameters for simulation
  
  h.par<-assempar(data = dat.h, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  h.sim<-simdata(Par = h.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set 
  h.samp<-sampsd(dat.sim = h.sim, Par = h.par, transformation = "fourth root",
                 method = "bray", multi.site = TRUE, n = n, p.n = p.n,
                 sites = sites, p.s = p.s, k = k)
  
  # average of multse for each potential sampling design
  h.sum<-summary_ssp(h.samp)
  h.sum$Strata<-rep("hightide", nrow(h.sum))
  
  
  ##SSP for midtide
  #parameters for simulation
  
  m.par<-assempar(data = dat.m, type = "cover", Sest.method = "chao")
  
  # Simulation of 100 data sets, each one with 100 potential sampling units  and 30 sites
  
  m.sim<-simdata(Par = m.par, cases = cases, n = n, sites = sites)
  
  # Sampling and estimation of multse for each data set 
  m.samp<-sampsd(dat.sim = m.sim, Par = m.par, transformation = "fourth root",
                 method = "bray", multi.site = TRUE, n = n, p.n = p.n,
                 sites = sites, p.s = p.s, k = k)
  
  
  # average of multse for each potential sampling design
  m.sum<-summary_ssp(m.samp)
  m.sum$Strata<-rep("midtide", nrow(m.sum))
  
  #Combine multse from simulations in a single data frame
  
  
  result<-rbind(h.sum, m.sum)
  result$locality<-rep(as.character(dat$Locality[1]), nrow(result))
  write.csv(result, file = paste(result$locality[1], ".csv"))
  print(paste("loc", i))
  
}

#cut and paste all the result files in the folder "SSP.output"

########## Plot MultSE ###################################################
#Set as working directory the folder "SSP.output"
setwd("~/asymmetric/SSP_P2P/SSP.output")

output<-list.files()

multse<-vector("list", length = length(output))

for (i in 1:length(output)){
  multse[[i]]<-read.csv(output[i]) 
}

multse.all<-do.call(rbind, multse)
multse.all$Strata<-factor(multse.all$Strata, levels = c("hightide","midtide", "lowtide"))

#Plot for sampling effort
multse.samples<-filter(multse.all, sv == "samples")

plot.1<-ggplot(multse.samples, aes(x=samples, y=mean, colour = locality))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
  facet_grid(.~Strata)+
  theme_bw(base_size=16) +
  ylab ("Multivariate pseudo SE")+ 
  xlab("Sampling effort")+
  scale_y_continuous(breaks=seq(0.0, 0.4, 0.025))+
  scale_x_continuous(breaks=seq(2, 50, 2))+
  theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
        axis.text.y = element_text(colour="black", size=rel(0.7)),
        axis.title.x = element_text(colour="black", size=rel(0.9)),
        axis.title.y = element_text(colour="black", size=rel(0.9)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4),
        axis.ticks= element_line(size=0.2))

plot.1  + geom_abline(slope = 0, intercept = 0.1, colour = "blue")

#Plot for number of sites
multse.sites<-filter(multse.all, sv == "sites")

plot.2<-ggplot(multse.sites, aes(x=samples, y=mean, colour = locality))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
  facet_grid(.~Strata)+
  theme_bw(base_size=16) +
  ylab ("Multivariate pseudo SE")+ 
  xlab("Number of sites")+
  scale_y_continuous(breaks=seq(0.0, 0.4, 0.025))+
  scale_x_continuous(breaks=seq(2, 20, 1))+
  theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
        axis.text.y = element_text(colour="black", size=rel(0.7)),
        axis.title.x = element_text(colour="black", size=rel(0.9)),
        axis.title.y = element_text(colour="black", size=rel(0.9)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4),
        axis.ticks= element_line(size=0.2))

plot.2+ geom_abline(slope = 0, intercept = 0.1, colour = "blue")

