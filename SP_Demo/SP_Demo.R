########### Demonstration of support points ############
library(fields)
library(mvtnorm)
library(SPlit)
library(ggplot2)
source("k_DPP/helper-functions.R")
source("k_DPP/kdpp-sample.R")
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")

### Generate full data for location set 1
set.seed(100)
X=cbind(runif(2e4),runif(2e4))
thetat=c(1.5,0.168639,0.27,1.5)
C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
set.seed(1)
Y=rmvnorm(1,mean=rep(0,2e4),sigma=C,method="chol")
dataF=data.frame(lon=X[,1],lat=X[,2],z=c(Y))

# SP
set.seed(1)
spres=sp(100,2,dist.samp=X)
spid=subsample(X,spres$sp)
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

# Rand
set.seed(1)
rid=sample(1:2e4,100,replace=FALSE)
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

# Grid
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

# k-DPP
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
# 4.30*4.30


### Generate full data for location set 2
set.seed(100)
X=rbind(cbind(runif(round((2e4)*2/3),0,0.5),runif(round((2e4)*2/3),0,0.5)),
            cbind(runif((2e4)*1/3+1),runif((2e4)*1/3+1)))
thetat=c(1.5,0.168639,0.27,1.5)
C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
set.seed(1)
Y=rmvnorm(1,mean=rep(0,2e4),sigma=C,method="chol")
dataF=data.frame(lon=X[,1],lat=X[,2],z=c(Y))

# SP
set.seed(1)
spres=sp(100,2,dist.samp=X)
spid=subsample(X,spres$sp)
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


# Rand
set.seed(1)
rid=sample(1:2e4,100,replace=FALSE)
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

# Grid
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

# k-DPP
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
# 4.30*4.30
