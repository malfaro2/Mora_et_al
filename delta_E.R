library(gdata)
library(tidyverse)

fn1 <- function(t){
  ifelse(t > (6/29)^3, t^(1/3),(t/(3*(6/29)^2))+ (4/29))
}

Lab <- function(x,y,z,xn,yn,zn){
  L <- 116 * fn1(y/yn)-16
  a <- 500 * (fn1(x/xn)-fn1(y/yn))
  b <- 200 * (fn1(y/yn)-fn1(z/zn))
  return(c(L,a,b))
}

delta_E <- function(espectro1,espectro2,iluminante,xout,cie){

  # step 1: interpolate ilum, cie y spectrum on the same 
  # range, and then plot to check.

  all_inter<-cbind(spline(cie$V1,cie$V2,xout=xout)[[1]],
                   spline(cie$V1,cie$V2,xout=xout)[[2]],
                   spline(cie$V1,cie$V3,xout=xout)[[2]],
                   spline(cie$V1,cie$V4,xout=xout)[[2]],
                   spline(iluminante$V1,iluminante$V2,xout=xout)[[2]],
                   spline(espectro1$X,espectro1$Y, xout=xout)[[2]],
                   spline(espectro2$X,espectro2$Y, xout=xout)[[2]])
  colnames(all_inter) <- c("lambda","cieX","cieY","cieZ","illum","spectra1","spectra2")
  all_inter <- data.frame(all_inter)
  
  # par(mfrow=c(1,2))
  # matplot(all_inter[,1],all_inter[,c(2:4)], type='l', ylab="CIE: x,y,z")
  # legend(700,1.2 , c("x", "y","z"),pch = "ooo", col = c(1,2,3))
  # 
  # matplot(all_inter[,1],all_inter[,c(6:7)], type='l', ylab="Spectra")
  # legend(500,30, c("A", "B"),pch = "ooo", col = c(1,2))
  
  # step 2: Calculate k
  k = 100 / sum(all_inter$illum*all_inter$cieY)
  
  # step 3: calculate x,y,z for each spectra and the base
  Xn <- k * sum(all_inter$illum*all_inter$cieX) 
  Yn <- k * sum(all_inter$illum*all_inter$cieY)
  Zn <- k * sum(all_inter$illum*all_inter$cieZ)
  totn <- Xn + Yn + Zn
  
  xn <- Xn/totn
  yn <- Yn/totn
  zn <- Zn/totn
  
  X1 <- k * sum(all_inter$spectra1*all_inter$illum*all_inter$cieX/100)
  Y1 <- k * sum(all_inter$spectra1*all_inter$illum*all_inter$cieY/100)
  Z1 <- k * sum(all_inter$spectra1*all_inter$illum*all_inter$cieZ/100)
  tot1 <- X1 + Y1 + Z1

  x1 <- X1/tot1
  y1 <- Y1/tot1
  z1 <- Z1/tot1
  
  X2 <- k * sum(all_inter$spectra2*all_inter$illum*all_inter$cieX/100)
  Y2 <- k * sum(all_inter$spectra2*all_inter$illum*all_inter$cieY/100)
  Z2 <- k * sum(all_inter$spectra2*all_inter$illum*all_inter$cieZ/100)
  tot2 <- X2 + Y2 + Z2
  
  x2 <- X2/tot2
  y2 <- Y2/tot2
  z2 <- Z2/tot2
  
  # step 4: calculate L,a,b for each spectra
  sp1<-Lab(x1,y1,z1,xn,yn,zn)
  sp2<-Lab(x2,y2,z2,xn,yn,zn)
  
  # step 5: calculate the squared difference
  
  deltaE <- sqrt(sum((sp1-sp2)^2))
  
  return(deltaE)
}

# Data needed:

cie<-read.xls("ciexyz31_1.xls", header=FALSE)
iluminante<-read.xls("illuminant.xls", header=FALSE)
xout<- seq(420,800,0.1)
load(file="dat.Rdata")

# Test:
caso1 <- which(dat$especie=="a"&dat$BOB=="black"&dat$spot=="s1"&dat$muestra==1)
caso2 <- which(dat$especie=="b"&dat$BOB=="black"&dat$spot=="s1"&dat$muestra==1)
espectro_ej1<-dat[caso1,c("X","Y")]
espectro_ej2<-dat[caso2,c("X","Y")]
delta_E(espectro_ej1,espectro_ej2,iluminante,xout,cie)

## Create the data set with all the pairs:

datg <- dat %>% mutate(especie = factor(especie, 
        labels = c("AC","BA","CR","EV","LA",
                   "MA","SC","SM",'TR')),
        BOB = factor(BOB, 
                         labels = c("BL","OR")))

datg <- datg %>% filter(especie!="EV") %>% 
  mutate(treat = paste0(especie,BOB,spot,muestra))

### All the combinations:

all <- data.frame(t(combn(unique(datg$treat),2)))
names(all) <- c('spectra1','spectra2')

### All delta E: 18336
### Check if they have the same structure:

# fac <-list(unique(datg$especie),unique(datg$spot),
#            unique(datg$muestra),unique(datg$BOB))

datos <- tibble(deltaE = sapply(1:dim(all)[1],function(i){
  delta_E(
 (datg %>% filter(treat==all[i,1]))[,c("X","Y")],
 (datg %>% filter(treat==all[i,2]))[,c("X","Y")],
 iluminante,xout,cie)}),
      spectra1 = as.character(all[1:dim(all)[1],1]),
      spectra2 = as.character(all[1:dim(all)[1],2]))

#plot(sort(datos$deltaE))
#abline(h=2.3, col="red",lty=2)

length(which(datos$deltaE < 2.3)) / (dim(datos)[1])

# spe1: specie 1
# spe2: specie 2
# speq: are specie 1 and 2 the same? 0 no 1 yes
# col1: BOB 1
# col2: BOB 2
# colq: are BOB 1 and 2 the same? 0 no 1 yes
# spo1: spot 1
# spo2: spot 2
# spoq: are spot 1 and 2 the same? 0 no 1 yes
# sam1: sample 1
# sam2: sample 2
# samq: are sample 1 and 2 the same? 0 no 1 yes

all_data <- datos %>% mutate(spe1 = substr(datos$spectra1,1,2),
                 spe2 = substr(datos$spectra2,1,2),
                 col1 = substr(datos$spectra1,3,4),
                 col2 = substr(datos$spectra2,3,4),
                 spo1 = substr(datos$spectra1,5,6),
                 spo2 = substr(datos$spectra2,5,6),
                 sam1 = substr(datos$spectra1,7,7),
                 sam2 = substr(datos$spectra2,7,7)) %>% 
  mutate(speq = ifelse(spe1==spe2,0,1),
         colq = ifelse(col1==col2,0,1),
         spoq = ifelse(spo1==spo2,0,1),
         samq = ifelse(sam1==sam2,0,1))
         
all_datac <- all_data
cambio <- which(all_data$spe1=="SM"&all_data$spe2=="SC")
all_datac$spe1[cambio] <- all_data$spe2[cambio]
all_datac$spe2[cambio] <- all_data$spe1[cambio] 
all_datac$col1[cambio] <- all_data$col2[cambio]
all_datac$col2[cambio] <- all_data$col1[cambio] 
all_datac$spo1[cambio] <- all_data$spo2[cambio]
all_datac$spo2[cambio] <- all_data$spo1[cambio] 
all_datac$sam1[cambio] <- all_data$sam2[cambio]
all_datac$sam2[cambio] <- all_data$sam1[cambio] 

all_data <- all_datac
save(all_data, file="deltaEdata.Rdata")


