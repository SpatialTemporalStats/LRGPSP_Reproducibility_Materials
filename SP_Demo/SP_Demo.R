################################################################################
# This file includes all steps to reproduce Figures 1 and S2                   #
################################################################################
# Necessary packages and functions
# library(fields)
# library(mvtnorm)
# library(support)
# library(SPlit)
# library(magrittr)
# library(ggplot2)
# source(here("Rep_Points/k_DPP","helper-functions.R"))
# source(here("Rep_Points/k_DPP","kdpp-sample.R"))
# source(here("Rep_Points","SP-sample.R"))
# source(here("Rep_Points","Rand-sample.R"))
# source(here("Rep_Points","Grid-sample.R"))


### Generate full data for the first location set
# t1=proc.time()[[3]]
set.seed(200)
X=cbind(runif(2e4),runif(2e4))
thetat=c(1.5,0.168639,0.27,1.5)
C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
set.seed(1)
Y=rmvnorm(1,mean=rep(0,2e4),sigma=C,method="chol")
dataF=data.frame(lon=X[,1],lat=X[,2],z=c(Y))
# t2=proc.time()[[3]]
# t2-t1=78.288


### Select k=100 SPs and plot Figure 1(a)
# t1=proc.time()[[3]]
spid=sp.sample(100,X,1)
# t2=proc.time()[[3]]
# t2-t1=1.815
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=dataF[spid,])+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Select k=100 Rands and plot Figure 1(b)
# t1=proc.time()[[3]]
rid=rand.sample(100,X,1)
# t2=proc.time()[[3]]
# t2-t1=0.004
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=dataF[rid,])+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Select k=100 k-DPPs and plot Figure S2(a)
# t1=proc.time()[[3]]
L=decompose.matrix(C)
set.seed(1)
dppid=k.dpp.sample(L,100)
# t2=proc.time()[[3]]
# t2-t1=1073.491
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=dataF[dppid,])+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Select k=100 Grids and plot Figure S2(b)
# t1=proc.time()[[3]]
gr=seq(0,1,length=12)[-c(1,12)]
Gr=cbind(rep(gr,each=10),rep(gr,times=10))
# t2=proc.time()[[3]]
# t2-t1=0.008
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=data.frame(lon=Gr[,1],lat=Gr[,2]))+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Generate full data for the second location set
set.seed(2)
X=rbind(cbind(runif(round((2e4)*2/3),0,0.5),runif(round((2e4)*2/3),0,0.5)),
            cbind(runif((2e4)*1/3+1),runif((2e4)*1/3+1)))
thetat=c(1.5,0.168639,0.27,1.5)
C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
set.seed(1)
Y=rmvnorm(1,mean=rep(0,2e4),sigma=C,method="chol")
dataF=data.frame(lon=X[,1],lat=X[,2],z=c(Y))


### Select k=100 SPs and plot Figure 1(c)
spid=sp.sample(100,X,1)
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=dataF[spid,])+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Select k=100 Grids and plot Figure 1(d)
gr=seq(0,1,length=12)[-c(1,12)]
Gr=cbind(rep(gr,each=10),rep(gr,times=10))
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=data.frame(lon=Gr[,1],lat=Gr[,2]))+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Select k=100 k-DPPs and plot Figure S2(c)
L=decompose.matrix(C)
set.seed(1)
dppid=k.dpp.sample(L,100)
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=dataF[dppid,])+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Select k=100 Rands and plot Figure S2(d)
rid=rand.sample(100,X,1)
PT=ggplot()+xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=17,data=dataF[rid,])+
  scale_colour_gradientn(colors = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "none",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


