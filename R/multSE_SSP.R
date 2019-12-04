####-----------------------------------------------------------------------------------------------
####        OPTIMIZATION OF ROCKY SHORE SAMPLING FROM P2P   / multSE using SSP                   #
####-----------------------------------------------------------------------------------------------

#####-----------------------------------------------------------------------------------------------

# Load required libraries
library(SSP)
library(sampling)
library(vegan)
library(readxl)

#### Data from Argentina - Puerto Madryn
ARG.PMADRYN <- read.csv("C:/Users/Edlin/OneDrive/Documents/MBON - P2P/mbon_pole2pole/data/ARG-PMADRYN.csv")

ARG.PMADRYN[is.na(ARG.PMADRYN)]<-0

##split the data by strata
dat.h<-filter(ARG.PMADRYN, Strata == "hightide")
dat.m<-filter(ARG.PMADRYN, Strata == "midtide")
dat.l<-filter(ARG.PMADRYN, Strata == "lowtide")

##adjust the data for the format required by SSP

dat.h<-dat.h[,c(3,5:length(dat.h))]
dat.m<-dat.m[,c(3,5:length(dat.m))]
dat.l<-dat.l[,c(3,5:length(dat.l))]

##SSP for hightide

#parameters for simulation

h.par<-assempar(data = dat.h, type = "cover", Sest.method = "chao")

# Simulation of 1o data sets, each one with 100 potential sampling units  and 10 sites

h.sim<-simdata(Par = h.par, cases = 10, n = 100, sites = 10)

# Sampling and estimation of multse for each data set (this step will take 6-8 minutes)
h.samp<-sampsd(dat.sim = h.sim, Par = h.par, transformation = "fourth root",
               method = "bray", multi.site = TRUE, n = 100, p.n = 30,
               sites = 10, p.s = 6, k = 10)

# average of multse for each potential sampling design
h.sum<-summary_ssp(h.samp)
h.sum$Strata<-rep("hightide", nrow(h.sum))


##SSP for midtide
#parameters for simulation

m.par<-assempar(data = dat.m, type = "cover", Sest.method = "chao")

# Simulation of 1o data sets, each one with 100 potential sampling units  and 10 sites

m.sim<-simdata(Par = m.par, cases = 10, n = 100, sites = 10)

# Sampling and estimation of multse for each data set (this step will take 6-8 minutes)
m.samp<-sampsd(dat.sim = m.sim, Par = m.par, transformation = "fourth root",
               method = "bray", multi.site = TRUE, n = 100, p.n = 30,
               sites = 10, p.s = 6, k = 10)


# average of multse for each potential sampling design
m.sum<-summary_ssp(m.samp)
m.sum$Strata<-rep("midtide", nrow(m.sum))

##SSP for lowtide

#parameters for simulation

l.par<-assempar(data = dat.l, type = "cover", Sest.method = "chao")

# Simulation of 1o data sets, each one with 100 potential sampling units  and 10 sites

l.sim<-simdata(Par = l.par, cases = 10, n = 100, sites = 10)

# Sampling and estimation of multse for each data set (this step will take 6-8 minutes)
l.samp<-sampsd(dat.sim = l.sim, Par = l.par, transformation = "fourth root",
               method = "bray", multi.site = TRUE, n = 100, p.n = 30,
               sites = 10, p.s = 6, k = 10)

# average of multse for each potential sampling design
l.sum<-summary_ssp(l.samp)
l.sum$Strata<-rep("lowtide", nrow(l.sum))

#Combine multse from simulations in a single data frame
multse.pmad<-rbind(h.sum, m.sum, l.sum)

multse.pmad$locality<-rep("Puerto Madryn", nrow(multse.pmad))

######################################################################################3

#### Data from Brazil - Arraial do Cabo
BRA.ARRAIALDOCABO <- read.csv("C:/Users/Edlin/OneDrive/Documents/MBON - P2P/mbon_pole2pole/data/BRA-ARRAIALDOCABO.csv")

BRA.ARRAIALDOCABO[is.na(BRA.ARRAIALDOCABO)]<-0

##split the data by strata
dat.h<-filter(BRA.ARRAIALDOCABO, Strata == "hightide")
dat.m<-filter(BRA.ARRAIALDOCABO, Strata == "midtide")
dat.l<-filter(BRA.ARRAIALDOCABO, Strata == "lowtide")

##adjust the data for the format required by SSP

dat.h<-dat.h[,c(3,5:length(dat.h))]
dat.m<-dat.m[,c(3,5:length(dat.m))]
dat.l<-dat.l[,c(3,5:length(dat.l))]

##SSP for hightide

#parameters for simulation

h.par<-assempar(data = dat.h, type = "cover", Sest.method = "chao")

# Simulation of 1o data sets, each one with 100 potential sampling units  and 10 sites

h.sim<-simdata(Par = h.par, cases = 10, n = 100, sites = 10)

# Sampling and estimation of multse for each data set (this step will take 6-8 minutes)
h.samp<-sampsd(dat.sim = h.sim, Par = h.par, transformation = "fourth root",
               method = "bray", multi.site = TRUE, n = 100, p.n = 30,
               sites = 10, p.s = 6, k = 10)

# average of multse for each potential sampling design
h.sum<-summary_ssp(h.samp)
h.sum$Strata<-rep("hightide", nrow(h.sum))


##SSP for midtide
#parameters for simulation

m.par<-assempar(data = dat.m, type = "cover", Sest.method = "chao")

# Simulation of 1o data sets, each one with 100 potential sampling units  and 10 sites

m.sim<-simdata(Par = m.par, cases = 10, n = 100, sites = 10)

# Sampling and estimation of multse for each data set (this step will take 6-8 minutes)
m.samp<-sampsd(dat.sim = m.sim, Par = m.par, transformation = "fourth root",
               method = "bray", multi.site = TRUE, n = 100, p.n = 30,
               sites = 10, p.s = 6, k = 10)


# average of multse for each potential sampling design
m.sum<-summary_ssp(m.samp)
m.sum$Strata<-rep("midtide", nrow(m.sum))

##SSP for lowtide

#parameters for simulation

l.par<-assempar(data = dat.l, type = "cover", Sest.method = "chao")

# Simulation of 1o data sets, each one with 100 potential sampling units  and 10 sites

l.sim<-simdata(Par = l.par, cases = 10, n = 100, sites = 10)

# Sampling and estimation of multse for each data set (this step will take 6-8 minutes)
l.samp<-sampsd(dat.sim = l.sim, Par = l.par, transformation = "fourth root",
               method = "bray", multi.site = TRUE, n = 100, p.n = 30,
               sites = 10, p.s = 6, k = 10)

# average of multse for each potential sampling design
l.sum<-summary_ssp(l.samp)
l.sum$Strata<-rep("lowtide", nrow(l.sum))

#Combine multse from simulations in a single data frame
multse.adc<-rbind(h.sum, m.sum, l.sum)

multse.adc$locality<-rep("Arraial do Cabo", nrow(multse.adc))

##########################################################################




multse.all<-rbind(multse.pmad, multse.adc)

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
  
plot.1  + geom_abline(slope = 0, intercept = 0.1)
