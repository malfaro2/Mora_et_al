# Analysis of distance per component - RGB - :

components <- function(espectro1,espectro2,xout,cie){
  
  # step 1: interpolate cie y spectra in the same range
  # plot to check
  
  all_inter<-cbind(spline(cie$V1,cie$V2,xout=xout)[[1]],
                   spline(cie$V1,cie$V2,xout=xout)[[2]],
                   spline(cie$V1,cie$V3,xout=xout)[[2]],
                   spline(cie$V1,cie$V4,xout=xout)[[2]],
                   spline(espectro1$X,espectro1$Y, xout=xout)[[2]],
                   spline(espectro2$X,espectro2$Y, xout=xout)[[2]])
  colnames(all_inter) <- c("lambda","cieX","cieY","cieZ","spectra1","spectra2")
  all_inter <- data.frame(all_inter)
  
  # step 2: calculate components x,y,z for each spectra (R G B)
  
  X1rojo <- all_inter$spectra1*all_inter$cieX / 
    sum(all_inter$spectra1*all_inter$cieX)
  Y1verde <- all_inter$spectra1*all_inter$cieY / 
    sum(all_inter$spectra1*all_inter$cieY)
  Z1azul <- all_inter$spectra1*all_inter$cieZ / 
    sum(all_inter$spectra1*all_inter$cieZ)
  
  X2rojo <- all_inter$spectra2*all_inter$cieX / 
    sum(all_inter$spectra2*all_inter$cieX)
  Y2verde <- all_inter$spectra2*all_inter$cieY / 
    sum(all_inter$spectra2*all_inter$cieY)
  Z2azul <- all_inter$spectra2*all_inter$cieZ / 
    sum(all_inter$spectra2*all_inter$cieZ)
  
  rojo  <- sum(abs(X1rojo-X2rojo))
  verde <- sum(abs(Y1verde-Y2verde))
  azul  <- sum(abs(Z1azul-Z2azul))
  
  out <- c(rojo, verde, azul)
  #names(out) <- c("difROJO","difVERDE","difAZUL")
  return(out)
}

# Data:

cie<-read.xls("ciexyz31_1.xls", header=FALSE)
xout<- seq(420,800,0.1)
load(file="dat.Rdata")

## Create the data set with all the pairs:

datg <- dat %>% mutate(especie = factor(especie, 
        labels = c("AC","BA","CR","EV","LA",
                    "MA","SC","SM",'TR')),
        BOB = factor(BOB,labels = c("BL","OR")))
datg <- datg %>%  filter(especie != "EV")
datg <- datg %>% mutate(treat = 
          paste0(especie,BOB,spot,muestra))

### All the combinations:
all <- data.frame(t(combn(unique(datg$treat),2)))
names(all) <- c('spectra1','spectra2')

## Distance between curves:

col_comp <- as.data.frame(t(sapply(1:dim(all)[1],function(i){
  components(
    (datg %>% filter(treat==all[i,1]))[,c("X","Y")],
    (datg %>% filter(treat==all[i,2]))[,c("X","Y")],
    xout,cie)})))

data_dist <- data.frame(cbind(col_comp,
  spectra1 = as.character(all[1:dim(all)[1],1]),
  spectra2 = as.character(all[1:dim(all)[1],2])))

names(data_dist) <- c("dist_rojo", "dist_verde", "dist_azul", 
                      "spectra1", "spectra2")
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

all_data <- data_dist %>% mutate(spe1 = substr(data_dist$spectra1,1,2),
                             spe2 = substr(data_dist$spectra2,1,2),
                             col1 = substr(data_dist$spectra1,3,4),
                             col2 = substr(data_dist$spectra2,3,4),
                             spo1 = substr(data_dist$spectra1,5,6),
                             spo2 = substr(data_dist$spectra2,5,6),
                             sam1 = substr(data_dist$spectra1,7,7),
                             sam2 = substr(data_dist$spectra2,7,7)) %>% 
  mutate(speq = ifelse(spe1==spe2,0,1),
         colq = ifelse(col1==col2,0,1),
         spoq = ifelse(spo1==spo2,0,1),
         samq = ifelse(sam1==sam2,0,1))

save(all_data, file="dist_data.Rdata")




