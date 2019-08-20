library(tidyverse)
library(reshape2)
library(fda.usc)
library(fdANOVA)
library(ERP)
library(fields)
library(ggforce)
library(xtable)
library(gdata)
library(ggpubr)

# Results

## Read the data:
load(file="dat.Rdata")
load(file="factores.Rdata")
aa<-dat %>% group_by(especie,BOB) %>% summarise(mean=mean(Y)) %>% arrange(especie,BOB) %>% dcast(especie ~ BOB, value.var = "mean") 
dat0 <- dat %>% dcast(X+especie+muestra+spot+dia~BOB, value.var = "Y") %>% mutate(dif=black-orange) %>% arrange(especie,muestra,spot)

# Change Laphita (LA) for Opisthacantha (OP)
# And, when possible use most of the name:
# \textit{Acanthoscelio} (AC), \textit{Baryconus} (BA), \textit{Chromoteleia} (CR), \textit{Macroteleia} (MA), \textit{Opisthacantha} (OP), \textit{Scelio} (SC), \textit{Sceliomorpha} (SM) and \textit{Triteleia}

orord2<-c("AC","BA","CR","EV",
          "OP","MA","SC","SM",'TR')
alord2<-c("AC","BA","CR","MA","OP",
          "SC","SM",'TR',"EV")

orord4<-c("Acanthoscelio","Baryconus","Chromoteleia",
          "Evaniella","Opisthacantha","Macroteleia",
          "Scelio","Sceliomorpha",'Triteleia')
alord4<-c("Acanthoscelio","Baryconus","Chromoteleia",
          "Macroteleia","Opisthacantha","Scelio",
          "Sceliomorpha",'Triteleia',"Evaniella")

## Figure 2:
datg <- dat %>% mutate(especie = 
        factor(especie, labels = orord4))%>% 
  mutate(especie = factor(especie, levels = alord4))
pp <- ggplot(datg, aes(X, Y, colour = BOB))+ 
  geom_line(size=0.9) + stat_summary(aes(group=BOB),
    fun.y = "mean", color="black",size = 0.5, 
    geom = "line") + 
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  theme_bw() + xlab(expression(lambda)) +
  ylab(expression(Z~"("~lambda~")")) +
  theme(axis.text.x = element_text(angle = 90))

tiff("Figures/Fig2.tiff", units="in", width=9, 
     height=3.2, res=300)
pp + facet_grid(BOB ~ especie) +
  theme(legend.position = "none") 
dev.off()

## Figure 3:
dat0g <- dat0 %>% mutate(genera = 
         factor(especie, labels = orord4))%>% 
  mutate(especie = factor(especie, levels = alord4))
pp1 <- ggplot(dat0g, aes(X, dif, colour = genera))+ 
  geom_line(size=0.01) +
  stat_summary(aes(group=genera, colour=genera),
    fun.y = "mean",size = 1, geom = "line") + 
  ylab("Difference Curves") +
  stat_summary(aes(group=genera, colour=genera),
    fun.y = "max",size = 0.15, geom = "line", linetype = "dashed") + 
  stat_summary(aes(group=genera, colour=genera),
    fun.y = "min",size = 0.15, geom = "line", linetype = "dashed") + 
  ylab("Difference Curves")
tiff("Figures/Fig3.tiff", units="in", width=7, 
     height=3, res=300)
pp1 + geom_hline(yintercept=0, linetype="dashed", 
                 color="black", size=1)+ 
  theme_bw() + xlab(expression(lambda)) +
  ylab(expression(Y~"("~lambda~")"))
dev.off()

## Calculate $\Delta E$ and calculate Table 2
#source(delta_E.R)
load(file="deltaEdata.Rdata")
aa<-lm(deltaE ~ as.factor(speq) * as.factor(colq) + 
         as.factor(samq) + as.factor(spoq), data=all_data)
summary(aa)
xtable((aa))

all_data$spe1<-factor(all_data$spe1,
                      labels=orord2[-4])
all_data$spe1<-factor(all_data$spe1,
                      levels=alord2[-9])

all_data$spe2<-factor(all_data$spe2,
                      labels=orord2[-4])
all_data$spe2<-factor(all_data$spe2,
                      levels=alord2[-9])

## Figure 4:
p <- ggplot(all_data, aes(x=as.factor(colq), 
                          y=deltaE, 
                          #color=as.factor(colq),
                          fill= as.factor(colq))) + 
  geom_violin() + facet_grid(spe1~spe2)+ theme_bw() +
  geom_boxplot(width=0.1) +
  labs(y= expression(Delta~"E"),
       fill= "",
       x= "From different color? No = 0, Yes = 1" ) + theme_bw() +
  geom_hline(yintercept=2.3, linetype="dashed", color = "red")
tiff("Figures/Fig4.tiff", units="in", width=8, 
     height=8, res=200)
p + theme(legend.position = "none")
dev.off()

## Figure 5:
cambio2 <- which(all_data$col1=="OR"&all_data$col2=="BL")
all_data2 <- all_data
all_data2$col1[cambio2]<-all_data$col2[cambio2]
all_data2$col2[cambio2]<-all_data$col1[cambio2]



all_data2 <-all_data2 %>% mutate(colone = 
factor(col1, labels = c("Black","Orange"))) %>%    
mutate(coltwo = 
factor(col2, labels = c("Black","Orange")))

q <- ggplot(all_data2, aes(x=as.factor(speq), 
                           y=deltaE, 
                           fill= as.factor(speq))) + 
  geom_violin() + facet_grid(colone~coltwo)+ theme_bw() +
  geom_boxplot(width=0.1) +
  labs(y=expression(Delta~"E"),
       fill= "",
       x= "From different genera? No = 0, Yes = 1" ) + theme_bw() +
  geom_hline(yintercept=2.3, linetype="dashed", color = "red")

tiff("Figures/Fig5.tiff", units="in", width=5, 
     height=5, res=300)
q + theme(legend.position = "none")
dev.off()

## Figure 8
dat0c <- dat0 %>% 
  mutate(genera = factor(especie, labels = orord4))%>% 
  mutate(especie = factor(especie, levels = alord4))
cc0 <- dat0c %>% group_by(genera, muestra, spot) %>% 
  summarise(Y1 = mean(dif))
bp <- ggplot(cc0, aes(x=genera, y=Y1, colour=genera)) + 
  geom_boxplot()+
  labs(title="",x="Genera", y = "Difference in Mean 
       Reflectance between Orange and Black") + 
  theme_bw() 
tiff("Figures/Fig8.tiff", units="in", width=7.5, 
     height=4, res=300)
bp + geom_hline(yintercept=0, linetype="dashed", 
                color = "black", size=1)+ theme(legend.position = "none")
dev.off()

# Table 5

cc <- dat0g %>% group_by(genera, muestra,spot) %>% 
  summarize(O = mean(orange), B = mean(black)) 
amm <- cc %>% group_by(genera) %>% 
  summarize(O = mean(O), B = mean(B)) 

cbind(amm,Difference=abs(amm$O-amm$B))[order(abs(amm$O-amm$B)),]


# FUNCTIONAL ANOVA:
## Reshape the data:
new<-dat0g
arvals <- unique(new$X)
mdata <- t(matrix(new$dif,nrow=length(arvals)))
fdatad <- fdata(mdata, arvals)
factores_new <- as.data.frame(  
  cbind(spot=dat0$spot[seq(1,127980,1185)],
        especie=dat0$especie[seq(1,127980,1185)],
        muestra=rep(1:36,each=3)))
image.plot(cov(mdata)/dim(mdata)[1],axes=F)
llab <- round(unique((dat0$X)[seq(1,1185,131)]),1)
axis(2, at=seq(0,1,1/9),labels=llab)
axis(1, at=seq(0,1,1/9),labels=llab)

## Functional ANOVA:

design  = model.matrix( ~muestra+spot+especie,data=factores_new)
design0 = model.matrix( ~muestra+especie,data=factores_new)
group <- as.factor(factores_new[,2])
ffdata <- (fdatad$data) 
test1 = erpFtest(dta=ffdata,design=design,design0=design0, svd.method=c("fast.svd"))
test1$pval
## no difference between spots / no variability to explain.

design  = model.matrix( ~muestra+especie,data=factores_new)
design0 = model.matrix( ~especie,data=factores_new)
test2 = erpFtest(dta=ffdata,design=design,design0=design0)
test2$pval
## difference intra-genera controlling for everything else

group <- as.factor(factores_new[,2])
design  = model.matrix( ~muestra+especie,data=factores_new)
design0 = model.matrix( ~muestra,data=factores_new)
test3 = erpFtest(dta=ffdata,design=design,design0=design0)
test3$pval
## difference inter-genera controlling for everything else

#Another (very slow) alternative:
#aa<-(fanova <- fanova.tests(x = t(mdata), group.label=group))
#save(aa,file="resFANOVA.Rdata")

load(file="resFANOVA.Rdata")
summary(aa)
#p.value = 0.00000

## Constrasts:
dat0 <- dat0 %>% mutate(genera = factor(especie, 
  labels = orord4))
group <- rep(unique(dat0$genera),each=12)
a<-plotFANOVA(x = t(mdata), group.label=group, means=TRUE)
a 
m0=data.frame(group)
cr5=contr.treatment(9, contrast=TRUE,base=9)
resul03c1=anova.RPm(fdatad,~group,m0,contrast=list(group=cr5))
summary.anova(resul03c1)


## Ranking Tables (3,4,6,7)

table1 <-all_data %>% 
  group_by(spectra1, spectra2) %>%
  summarise(mean1=mean(deltaE)) %>%arrange(mean1) %>% 
  spread(spectra1, mean1, fill = 0)


aam1 <- all_data %>% 
  group_by(spe1, spe2, col1, col2) %>% summarize(mdistance=mean(deltaE)) %>% 
  arrange(desc(mdistance),spe1,col1)

mismos_col <- which(aam1$col1 == aam1$col2)
mismos_gen <- which(aam1$spe1 == aam1$spe2)

summary(aam1$mdistance)
xtable(head(aam1[mismos_col,], n=10))
xtable(head(aam1[-mismos_col,], n=10))
xtable(head(aam1[mismos_gen,], n=10))
xtable(head(aam1[-mismos_gen,], n=10))

## Component analysis

# Calculate the distance between curves per component
# source(file=components.R)
load(file="dist_data.Rdata")
orord4
selec <- all_data$spe1%in%c("AC", "BA", "TR")&
  all_data$spe2%in%c("AC", "BA", "TR")

# FIGURE 6

p <- ggplot(all_data[selec,], aes(x=as.factor(colq), 
                                  y=dist_verde, 
                                  fill= as.factor(colq))) + 
  geom_violin() + facet_grid(spe1~spe2)+ theme_bw() +
  geom_boxplot(width=0.1) + ylim(0, 0.4) + 
  labs(y="GREEN distance",
       fill= "",
       x= "From different color? NO = 0, YES = 1" ) + 
  theme_bw() + theme(legend.position = "none")

q <- ggplot(all_data[selec,], aes(x=as.factor(colq), 
                                  y=dist_rojo, 
                                  fill= as.factor(colq))) + 
  geom_violin() + facet_grid(spe1~spe2)+ theme_bw() +
  geom_boxplot(width=0.1) + ylim(0, 0.4) + 
  labs(y="RED distance",
       fill= "",
       x= "From different color? NO = 0, YES = 1" ) + 
  theme_bw() + theme(legend.position = "none")

r <- ggplot(all_data[selec,], aes(x=as.factor(colq), 
                                  y=dist_azul, 
                                  fill= as.factor(colq))) + 
  geom_violin() + facet_grid(spe1~spe2)+ theme_bw() +
  geom_boxplot(width=0.1) + ylim(0, 0.4) + 
  labs(y="BLUE distance",
       fill= "",
       x= "From different color? NO = 0, YES = 1" ) + 
  theme_bw() + theme(legend.position = "none")

## Gráfico de RGB: 

source(file="componentes2.R")
dat<-aa %>% gather() %>% mutate(x=rep(xout,12)) %>% 
  mutate(RGB=factor(key, 
      levels=c("ROJO1","ROJO2","ROJO3","ROJO4",
               "VERDE1", "VERDE2","VERDE3","VERDE4",
               "AZUL1", "AZUL2","AZUL3", "AZUL4"),
labels=c("Red AC-B","Red AC-O","Red TR-B","Red TR-O", 
"Green AC-B", "Green AC-O","Green TR-B","Green TR-O",
"Blue AC-B", "Blue AC-O","Blue TR-B", "Blue TR-O")))

datO<-dat %>% mutate(linea = 
      str_extract(RGB, "[^ ]+$"))

s <- ggplot(data=datO, aes(x=x,y=value, color=RGB)) +
  geom_line( ) + 
  theme_bw() +
  labs(x=expression(lambda), 
       y=expression(RGB(lambda))) + 
  theme(legend.justification=c("left", "top"),
        legend.position = c(.76, .95))

tiff("Figures/Fig6.tiff", units="in", width=12, 
     height=8, res=160)
ggarrange(q, p, r, s,
          labels = c("(A)", "(B)", "(C)", "(D)"),
          ncol = 2, nrow = 2)
dev.off()

# Area under the curve and maximum per component

# lambda in which we find maximum spectra:

datamean<- datall %>% 
  group_by(Component,SP,CO, Treat) %>% 
  filter(Spectra == max(Spectra))

# area under the curve for all:

datamean2<- datall %>% 
  group_by(Component,SP,CO, Treat) %>% 
  summarise(mSpec=sum(Spectra))

# Scatter plots

# maximum per mean genera:
datamean3 <- datall %>% 
  group_by(Component,SP,CO, Treat) %>% 
  mutate(areaSpec=sum(Spectra)) %>% 
  filter(Spectra == max(Spectra)) 

datamean4 <- datamean3%>%
  group_by(SP, Component, CO) %>%
  mutate(X=median(X),
  areaSpec=median(areaSpec)) %>%
  arrange(Component)

datamean4$genera<-factor(datamean4$SP,labels=orord2)
datamean4$genera<-factor(datamean4$genera,labels=orord4)

datamean4 <- datamean4 %>% 
mutate(genera = factor(genera, levels = alord4))

datamean3$genera<-factor(datamean3$SP,labels=orord2)
datamean3$genera<-factor(datamean3$genera,labels=orord4)

datamean3 <- datamean3 %>% 
  mutate(genera = factor(genera, levels = alord4))

## Figure 7

p1<-ggplot(datamean3[datamean3$CO=="BL",],
       aes(x=X, y=areaSpec, col=genera)) + 
  geom_point(alpha=1/3) + 
  facet_grid(.~Component, scales="free_x")+
  stat_summary(data=datamean4[datamean4$CO=="BL",],
               fun.y = "mean",size = 3, 
               geom = "point", shape=16)+
  scale_y_log10() +
  theme_bw()+
  ylab("log (Area under the curve)")+
  xlab(expression(lambda~"["~max~Spectra~"]"))
#  ggtitle("Scatter Plot, for all BLACK observations")

p2<-ggplot(datamean3[datamean3$CO=="OR",],
       aes(x=X, y=areaSpec, col=genera)) + 
  geom_point(alpha=1/3) + 
  facet_grid(.~Component, scales="free_x") +
  stat_summary(data=datamean4[datamean4$CO=="OR",],
               fun.y = "mean",size = 3, 
               geom = "point", shape=16)+
  scale_y_log10() +
  theme_bw() +
  ylab("log (Area under the curve)") +
  xlab(expression(lambda~"["~max~Spectra~"]"))
#  ggtitle("Scatter Plot, for all ORANGE observations")

p12<-ggplot(datamean3,
       aes(x=X, y=areaSpec, col=genera)) + 
  geom_point(alpha=1/5) + 
  facet_grid(CO~Component, scales="free")+
  stat_summary(data=datamean4,
               fun.y = "mean",size = 3, 
               geom = "point", shape=16)+
  scale_y_log10() +
  theme_bw()+
  ylab("log (Area under the curve)")+
  xlab(expression(lambda~"["~max~Spectra~"]"))


tiff("Figures/Fig7.tiff", units="in", width=5.8, 
     height=4.3, res=300)
p12
dev.off()
