## Code to read data

library(dplyr)
library(ggplot2)

temp = list.files(path="Datos_April/DIA1",pattern="*.csv")
## 51 files -> eliminate refscan y darkscan
temp = temp[c(1:24,26:37,39:56)]
for (i in 1:length(temp)) 
  assign(paste0(gsub(" ", "", temp[i], fixed = TRUE)), 
      read.csv(paste0("Datos_April/DIA1/",temp[i]),
                header=FALSE, skip=9))

temp = list.files(path="Datos_April/DIA2",pattern="*.csv")
## 51 files -> eliminate refscan y darkscan
temp = temp[c(1:24,26:37,39:56)]
for (i in 1:length(temp)) 
  assign(paste0(gsub(" ", "", temp[i], fixed = TRUE)), 
         read.csv(paste0("Datos_April/DIA2/",temp[i]),
                  header=FALSE, skip=9))

temp = list.files(path="Datos_April/DIA3",pattern="*.csv")
## 51 files -> eliminate refscan y darkscan
temp = temp[c(1:24,26:37,39:56)]
for (i in 1:length(temp)) 
  assign(paste0(gsub(" ", "", temp[i], fixed = TRUE)), 
         read.csv(paste0("Datos_April/DIA3/",temp[i]),
                  header=FALSE, skip=9))

temp = list.files(path="Datos_April/DIA4",pattern="*.csv")
## 51 files -> eliminate refscan y darkscan
temp = temp[c(1:24,26:37,39:56)]
for (i in 1:length(temp)) 
  assign(paste0(gsub(" ", "", temp[i], fixed = TRUE)), 
         read.csv(paste0("Datos_April/DIA4/",temp[i]),
                  header=FALSE, skip=9))

## Create a table with x, y, dia, especimen, BOB (B1, O, B2),especie:
rm(i, temp)
file_names <- ls()  
length(file_names)

dia     <- substring(file_names, 
                     regexpr("d", file_names), 
                     regexpr("d",file_names)+1)
especie <- substring(file_names, 1,
                     regexpr(".csv", file_names)-7)
spot    <- substring(file_names,
                     regexpr(".csv", file_names)-5,
                     regexpr(".csv", file_names)-4)
muestra <- substring(file_names,
                     regexpr(".csv", file_names)-3,
                     regexpr(".csv", file_names)-3)
BOB     <- substring(file_names,
                     regexpr(".csv", file_names)-6,
                     regexpr(".csv", file_names)-6)

BOB[BOB=="b"] <- "black"
BOB[BOB=="o"] <- "orange"

## Extract from 420 to 620 lambda (in total 1185 obs)
data <- list()
for(i in 1:length(file_names)){
  a<-get(file_names[i])
  d<-a[a$V1>420&a$V1<920,]
  n<-dim(d)[1]
  data[[i]]<-cbind(d,dia=rep(dia[i],n),
                   especie=rep(especie[i],n),
                   spot=rep(spot[i],n),
                   muestra=rep(muestra[i],n), 
                   BOB =rep(BOB[i],n))
  names(data[[i]]) <- c("X", "Y", "dia", "especie", 
                        "spot", "muestra","BOB")   
}
factores <- list(dia, especie, muestra,BOB)

length(data)
unique(unlist(lapply(data, dim)))

# There are $1185$ observations of X, Y, for each dia, 
# especie, spot, muestra and BOB (in total $216$ tables 
# or curves: $9*3*4*2$ 
# (especies+control)x(muestra)x(dia)x(BOB)).

# Complete table: 216 x 1185 x 7

dat <- do.call(rbind,data)
summary(dat)
save(dat, file="dat.Rdata")
save(factores, file="factores.Rdata")




