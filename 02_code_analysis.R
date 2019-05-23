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

## Figure 2:
tiff("Figures/Fig2.tiff", units="in", width=15, 
     height=5, res=300)
datg <- dat %>% mutate(especie = 
        factor(especie, labels = 
  c("AC","BA","CR","EV","LA","MA","SC","SM",'TR')))%>% 
  mutate(especie = factor(especie, levels = 
  c("AC","BA","CR","LA","MA","SC","SM",'TR',"EV")))
pp <- ggplot(datg, aes(X, Y, colour = BOB))+ 
  geom_line(size=0.9) + stat_summary(aes(group=BOB),
    fun.y = "mean", color="black",size = 0.5, 
    geom = "line") + 
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  theme_bw() + xlab(expression(lambda)) +
  ylab(expression(Z~"("~lambda~")")) +
  theme(axis.text.x = element_text(angle = 90))
pp + facet_grid(BOB ~ especie)
#pp +  facet_grid_paginate(BOB ~ especie, ncol = 3, 
#                          nrow = 2)#, page = 1)
#pp +  facet_grid_paginate(BOB ~ especie, ncol = 3, 
#                          nrow = 2)#, page = 2)
#pp +  facet_grid_paginate(BOB ~ especie, ncol = 3, 
#                          nrow = 2)#, page = 3)
dev.off()

## Figure 3:
dat0g <- dat0 %>% mutate(genera = factor(especie, labels = 
      c("AC","BA","CR","EV","LA","MA","SC","SM",'TR')))%>% 
  mutate(genera = factor(genera, levels = 
      c("AC","BA","CR","LA","MA","SC","SM",'TR',"EV")))
tiff("Figures/Fig3.tiff", units="in", width=10, 
     height=5, res=300)
pp1 <- ggplot(dat0g, aes(X, dif, colour = genera))+ 
  geom_line(size=0.09) +
  stat_summary(aes(group=genera, colour=genera),
    fun.y = "mean",size = 1, geom = "line") + 
  ylab("Difference Curves") 
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
tiff("Figures/Fig4.tiff", units="in", width=15, 
     height=15, res=300)
p
dev.off()

## Figure 5:
cambio2 <- which(all_data$col1=="OR"&all_data$col2=="BL")
all_data2 <- all_data
all_data2$col1[cambio2]<-all_data$col2[cambio2]
all_data2$col2[cambio2]<-all_data$col1[cambio2]

q <- ggplot(all_data2, aes(x=as.factor(speq), 
                           y=deltaE, 
                           #color=as.factor(colq),
                           fill= as.factor(speq))) + 
  geom_violin() + facet_grid(col1~col2)+ theme_bw() +
  geom_boxplot(width=0.1) +
  labs(y=expression(Delta~"E"),
       fill= "",
       x= "From different genera? No = 0, Yes = 1" ) + theme_bw() +
  geom_hline(yintercept=2.3, linetype="dashed", color = "red")

tiff("Figures/Fig5.tiff", units="in", width=5, 
     height=5, res=300)
q
dev.off()

## Figure 8
dat0c <- dat0 %>% 
  mutate(genera = factor(especie, labels = 
  c("AC","BA","CR","EV","LA","MA","SC","SM",'TR')))%>% 
  mutate(genera = factor(genera, levels = 
  c("AC","BA","CR","LA","MA","SC","SM",'TR',"EV")))
cc0 <- dat0c %>% group_by(genera, muestra, spot) %>% 
  summarise(Y1 = mean(dif))
bp <- ggplot(cc0, aes(x=genera, y=Y1, colour=genera)) + 
  geom_boxplot()+
  labs(title="",x="Genera", y = "Difference in Mean 
       Reflectance between Orange and Black") + 
  theme_bw() 
tiff("Figures/Fig8.tiff", units="in", width=10, 
     height=5, res=300)
bp + geom_hline(yintercept=0, linetype="dashed", 
                color = "black", size=1)
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
  labels = c("AC","BA","CR","EV","LA","MA","SC","SM",'TR')))
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
  theme_bw() 

q <- ggplot(all_data[selec,], aes(x=as.factor(colq), 
                                  y=dist_rojo, 
                                  fill= as.factor(colq))) + 
  geom_violin() + facet_grid(spe1~spe2)+ theme_bw() +
  geom_boxplot(width=0.1) + ylim(0, 0.4) + 
  labs(y="RED distance",
       fill= "",
       x= "From different color? NO = 0, YES = 1" ) + 
  theme_bw() 

r <- ggplot(all_data[selec,], aes(x=as.factor(colq), 
                                  y=dist_azul, 
                                  fill= as.factor(colq))) + 
  geom_violin() + facet_grid(spe1~spe2)+ theme_bw() +
  geom_boxplot(width=0.1) + ylim(0, 0.4) + 
  labs(y="BLUE distance",
       fill= "",
       x= "From different color? NO = 0, YES = 1" ) + 
  theme_bw() 

## Gráfico de RGB: 

dat<-aa %>% gather() %>% mutate(x=rep(xout,6)) %>% 
  mutate(RGB=factor(key, 
      levels=c("ROJO1","ROJO2", "VERDE1", 
               "VERDE2","AZUL1", "AZUL2"),
      labels=c("Red AC-B","Red SM-O", "Green AC-B", 
               "Green SM-O","Blue AC-B", "Blue SM-O")))

datO<-dat %>% mutate(linea = 
      str_extract(RGB, "[^ ]+$")) %>% 
  mutate(color = 
      substr(RGB, start = 1, stop = 4)) %>% 
  mutate(color = factor(color, levels=
      c("Blue","Red","Green")))

s <- ggplot(data=datO, aes(x=x,y=value, color=RGB)) +
  geom_line( ) + 
  theme_bw() +
  labs(x=expression(lambda), 
       y=expression(RGB(lambda))) + 
  theme(legend.justification=c("left", "top"),
        legend.position = c(.70, .85))

tiff("Figures/Fig6.tiff", units="in", width=12, 
     height=8, res=300)
ggarrange(p, q, r, s,
          labels = c("(A)", "(B)", "(C)", "(D)"),
          ncol = 2, nrow = 2)
dev.off()

# Area under the curve and maximum per component

source(file="componentes2.R")
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
  summarise(X=median(X),
            areaSpec=median(areaSpec)) %>%
  arrange(Component)

datamean4 <- datamean4 %>% 
mutate(genera = factor(SP, levels = 
       levels(factor(datamean4$SP)))) %>% 
mutate(genera = factor(genera, levels = 
       c("AC","BA","CR","LA","MA",
         "SC","SM",'TR',"EV")))

datamean3 <- datamean3 %>% 
  mutate(genera = factor(SP, levels = 
      levels(factor(datamean3$SP)))) %>% 
  mutate(genera = factor(genera, levels = 
       c("AC","BA","CR","LA","MA",
         "SC","SM",'TR',"EV")))

## Figure 7

p1<-ggplot(datamean3[datamean3$CO=="BL",],
       aes(x=X, y=areaSpec, col=genera)) + 
  geom_point(alpha=1/3) + 
  facet_grid(.~Component, scales="free_x")+
  stat_summary(data=datamean4[datamean4$CO=="BL",],
               fun.y = "mean",size = 4, 
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
               fun.y = "mean",size = 4, 
               geom = "point", shape=16)+
  scale_y_log10() +
  theme_bw() +
  ylab("log (Area under the curve)") +
  xlab(expression(lambda~"["~max~Spectra~"]"))
#  ggtitle("Scatter Plot, for all ORANGE observations")

tiff("Figures/Fig7.tiff", units="in", width=8, 
     height=8, res=300)
ggarrange(p1, p2,
          labels = c("(A)", "(B)"),
          ncol = 1, nrow = 2)
dev.off()
