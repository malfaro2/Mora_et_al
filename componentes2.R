
components <- function(espectro1,xout,cie,treat){
  
  # step 1: interpolate cie y spectra on the same range
  # plot to check
  
  all_inter<-cbind(spline(cie$V1,cie$V2,xout=xout)[[1]],
                   spline(cie$V1,cie$V2,xout=xout)[[2]],
                   spline(cie$V1,cie$V3,xout=xout)[[2]],
                   spline(cie$V1,cie$V4,xout=xout)[[2]],
                   spline(espectro1$X,espectro1$Y, xout=xout)[[2]])
  colnames(all_inter) <- c("lambda","cieX","cieY","cieZ","spectra1")
  all_inter <- data.frame(all_inter)
  
  # step 2: calculate componentes x,y,z for each spectra
  
  X1rojo <- all_inter$spectra1*all_inter$cieX 
  Y1verde <- all_inter$spectra1*all_inter$cieY 
  Z1azul <- all_inter$spectra1*all_inter$cieZ 
  
  out <- tibble(Spectra=as.double(c(X1rojo, Y1verde, Z1azul)),
                Component=rep(c("R","G","B"),each=length(xout)),
                Treat=rep(treat,(length(xout)*3)))
  return(out)
}

cie<-read.xls("ciexyz31_1.xls", header=FALSE)
xout<- seq(420,800,0.1)
load(file="dat.Rdata")

## Create the data set with all the treatments and without EV:

datg <- dat %>% mutate(especie = factor(especie, 
         labels = c("AC","BA","CR","EV","LA",
                    "MA","SC","SM",'TR')),
         BOB = factor(BOB, labels = c("BL","OR")))
#datg <- datg %>%  filter(especie != "EV")
#&                          especie != "CR")
datg <- datg %>% mutate(treat = paste0(especie,BOB,
                                       spot,muestra))

### For all cases:

all <- unique(datg$treat)
col_comp <- (lapply(1:length(all),function(i){
  components((datg %>% filter(treat==all[i]))
             [,c("X","Y")],xout,cie, all[i])}))

datall<-do.call(rbind,col_comp)
datall<-data.frame(cbind(datall,X= rep(xout,216)))
summary(datall)

#"Spectra"   
#"Component" 
#"Treat" 
# "X"

datall <- datall %>% mutate(SP = substr(Treat, start = 1, stop = 2),
                  CO = substr(Treat, start = 3, stop = 4))

datall<-as_tibble(datall)
head(datall)
