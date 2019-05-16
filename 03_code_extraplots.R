## Read the data:
load(file="dat.Rdata")
load(file="factores.Rdata")
aa<-dat %>% group_by(especie,BOB) %>% 
  summarise(mean=mean(Y)) %>% 
  arrange(especie,BOB) %>% 
  dcast(especie ~ BOB, value.var = "mean") 
dat0 <- dat %>% 
  dcast(X+especie+muestra+spot+dia~BOB, value.var = "Y") %>% 
  mutate(dif=black-orange) %>% 
  arrange(especie,muestra,spot)

## Figure X:

dat0g <- dat0 %>% mutate(genera = factor(especie, labels = 
        c("AC","BA","CR","EV","LA","MA","SC","SM",'TR')))%>% 
  mutate(genera = factor(genera, levels = 
        c("AC","BA","CR","LA","MA","SC","SM",'TR',"EV")))
pp1 <- ggplot(dat0g, aes(X, dif, colour = genera))+ 
  stat_summary(aes(group=genera, colour=genera),
               fun.y = "mean",size = 1, geom = "line") +
  ylab("Difference Curves") + theme_bw()
pp1 + geom_hline(yintercept=0, linetype="dashed", 
                 color="black", size=1)

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


## Functional box plots:
## Functional box plot:
library(fda)
fbplot(t(mdata),method='MBD')
abline(h=0,col="red", type="dashed")
fbplot(t(mdata),method='MBD', xlim=c(400,800))
abline(h=0,col="red", type="dashed")
fbplot(t(mdata),method='MBD', xlim=c(800,1200))
abline(h=0,col="red", type="dashed")

# Component analysis:

source(file=componentes2.R)

### Visualize components

ggplot(datall, 
       aes(X, Spectra, col=SP))+ geom_line(size=0.9)+
  stat_summary(aes(group=SP),fun.y = "mean", 
               color="black",size = 0.5, geom = "line") +
  theme_bw() + xlab(expression(lambda)) +  
  facet_grid(Component ~ CO)

### Visualize components
# Per component:
ggplot(datall[datall$Component=="R",], 
       aes(X, Spectra, col=SP))+ geom_line(size=0.9)+
  stat_summary(aes(group=SP),fun.y = "mean", 
               color="black",size = 0.5, geom = "line") + 
  theme_bw() + xlab(expression(lambda)) +  facet_grid(SP ~ CO)

ggplot(datall[datall$Component=="G",], 
       aes(X, Spectra, col=SP))+ geom_line(size=0.9)+
  stat_summary(aes(group=SP),fun.y = "mean", 
               color="black",size = 0.5, geom = "line") + 
  theme_bw() + xlab(expression(lambda)) +  facet_grid(SP ~ CO)

ggplot(datall[datall$Component=="B",], 
       aes(X, Spectra, col=SP))+ geom_line(size=0.9)+
  stat_summary(aes(group=SP),fun.y = "mean", 
               color="black",size = 0.5, geom = "line") + 
  theme_bw() + xlab(expression(lambda)) +  facet_grid(SP ~ CO)

# lambda in which we find maximum spectra:

datamean<- datall %>% 
  group_by(Component,SP,CO, Treat) %>% 
  filter(Spectra == max(Spectra))

ggplot(datamean[datamean$Component=="R"&
                  datamean$CO=="BL",], 
       aes(reorder(SP, X), X)) + 
  geom_point() +
  stat_summary(fun.y = "mean",size = 2, 
               geom = "point", col="red")+
  ggtitle("Red Component, Black Spot")

ggplot(datamean[datamean$Component=="R"&
                  datamean$CO=="OR",], 
       aes(reorder(SP, X), X)) + 
  geom_point() +
  stat_summary(fun.y = "mean",size = 2, 
               geom = "point", col="red")+
  ggtitle("Red Component, Orange Spot")

# area under the curve for all:

datamean2<- datall %>% 
  group_by(Component,SP,CO, Treat) %>% 
  summarise(mSpec=sum(Spectra))

ggplot(datamean2[datamean2$Component=="R"&
                   datamean2$CO=="OR",], 
       aes(reorder(SP, mSpec), mSpec)) + 
  geom_point() +
  stat_summary(fun.y = "mean",size = 2, 
               geom = "point", col="red")+
  ggtitle("Red Component, Orange Spot")

ggplot(datamean2[datamean2$Component=="R"&
                   datamean2$CO=="BL",], 
       aes(reorder(SP, mSpec), mSpec)) + 
  geom_point() +
  stat_summary(fun.y = "mean",size = 2, 
               geom = "point", col="red")+
  ggtitle("Red Component, Black Spot")


# Scatter plots:
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

ggplot(datamean3[datamean3$CO=="BL",],
       aes(x=1/X, y=areaSpec, col=SP)) + 
  geom_point(alpha=1/3) + facet_grid(.~Component, scales="free_x")+
  stat_summary(data=datamean4[datamean4$CO=="BL",],
               fun.y = "mean",size = 4, 
               geom = "point", shape=16)+
  scale_y_log10() +
  ylab("log (Area under the curve)")+
  xlab("1/Lambda[max Spectra]")+
  ggtitle("Scatter Plot, for all BLACK observations")

ggplot(datamean3[datamean3$CO=="OR",],
       aes(x=1/X, y=areaSpec, col=SP)) + 
  geom_point(alpha=1/3) + facet_grid(.~Component, scales="free_x")+
  stat_summary(data=datamean4[datamean4$CO=="OR",],
               fun.y = "mean",size = 4, 
               geom = "point", shape=16)+
  scale_y_log10() +
  ylab("log (Area under the curve)")+
  xlab("1/Lambda[max Spectra]")+
  ggtitle("Scatter Plot, for all ORANGE observations")


ggplot(datamean3,
       aes(x=1/X, y=areaSpec, col=SP)) + 
  geom_point(alpha=1/3) + facet_grid(CO~Component, scales="free_x")+
  stat_summary(data=datamean4,
               fun.y = "mean",size = 4, 
               geom = "point", shape=16)+
  scale_y_log10() +
  ylab("log (Area under the curve)")+
  xlab("Lambda[max Spectra]")+
  ggtitle("Scatter Plot, for all observations")


plotdat<-datamean3 %>%  
  gather(key=variable,value=value, -CO, -SP, -Component, -Treat) %>% 
  filter(variable != "areaSpec") %>%
  spread(Component,value=value,drop=TRUE)

plotdat<-datamean3 %>%  
  gather(key=variable,value=value, -CO, -SP, -Component, -Treat) %>% 
  filter(variable != "LmaxSpec") %>%
  spread(Component,value=value,drop=TRUE)

#### Delta color vs detal lambda[max(Spectra)]

datamean<- datall %>% 
  group_by(Component,SP,CO,Treat) %>% 
  mutate(areaSpec=sum(Spectra)*(0.1)) %>% 
  filter(Spectra == max(Spectra)) %>% 
  arrange(Component,SP,CO)

newdata<-tibble(
  diffSpec = filter(datamean, CO == "BL")$areaSpec  - 
    filter(datamean, CO == "OR")$areaSpec,
  difflamb = filter(datamean, CO == "BL")$X  - 
    filter(datamean, CO == "OR")$X,
  Component = filter(datamean, CO == "BL")$Component,
  SP = filter(datamean, CO == "BL")$SP)
newdata <- as_tibble(newdata)

datamean4 <- newdata%>%
  group_by(SP, Component) %>%
  summarise(difflamb=median(difflamb),
            diffSpec=median(diffSpec)) %>%
  arrange(Component)

ggplot(newdata,
       aes(x=difflamb, y=diffSpec, col=SP)) + 
  geom_point(alpha=1/3) + facet_grid(.~Component, scales="free_x")+
  stat_summary(data=datamean4,
               fun.y = "median",size = 4, 
               geom = "point", shape=16)+
  theme_bw()+
  ylab(expression(Delta*"Area under the Spectra BL vs OR"))+
  xlab(expression(Delta*lambda*"[max Spectra] BL vs OR"))



