
#Gorgonia: muy baja riqueza, una especie ubicua. La reportada en el trabajo no coincide con los datos
COL_GORGONA <- read_csv("SSP_P2P/multiple/COL-GORGONA.csv")
dat <- COL_GORGONA
dat[is.na(dat)]<-0

dat <- dat[,-c(1,2,3,4,5,11)]

library(vegan)
library(dplyr)
library(iNEXT)
library(betapart)

specpool(dat, pool = COL_GORGONA$Strata)
dat.pa <- decostand(dat, "pa")
apply(dat.pa,2,sum)

gorgona <- dat.pa %>% 
  group_by(COL_GORGONA$Strata) %>% 
  summarise_all(sum)


#ECO_GALAPAGOS: muy pocas especies
dat2 <- ECO_GALAPAGOS[,-c(1,2,3,4,14,15)]
dat2[is.na(dat2)]<-0
specpool(dat2, pool = ECO_GALAPAGOS$Strata)
dat2.pa <- decostand(dat2, "pa")

galapagos <- dat2.pa %>% 
  group_by(ECO_GALAPAGOS$Strata) %>% 
    summarise_all(sum)

#Viña del Mar

dat3 <- CHI.VINADELMAR_single[,-(1:4)]
dat3[is.na(dat3)]<-0
specpool(dat3, pool =  CHI.VINADELMAR_single$Strata)
dat3.pa <- decostand(dat2, "pa")

##############################################

#localidades con varios sitios
setwd("~/asymmetric/SSP_P2P/multiple")
files <- list.files()
datos<-NULL

for ( i in 1:length(files)){
  datos[[i]]<-read.csv(files[i])}

dat1 <- bind_rows(datos)

#localidades con un solo sitio
setwd("~/asymmetric/SSP_P2P/single")
files2 <- list.files()
datos2<-NULL

for ( i in 1:length(files2)){
  datos2[[i]]<-read.csv(files2[i])}

dat2 <- bind_rows(datos2)

#localidades USA
setwd("~/asymmetric/SSP_P2P/USA")
files3 <- list.files()
datos3<-NULL

for ( i in 1:length(files3)){
  datos3[[i]]<-read.csv(files3[i])}

dat3 <- bind_rows(datos3)

datos <- bind_rows(dat1, dat2, dat3)

datos$Locality <- as.factor(datos$Locality)
datos$Site <- as.factor(datos$Site)
datos$Strata <- as.factor(datos$Strata)

strata.l <- levels(datos$Strata)

#####################################################

ht <- datos[datos$Strata==strata.l[1],]
lt <- datos[datos$Strata==strata.l[2],]
mt <- datos[datos$Strata==strata.l[3],]








