########################################################################################
# This file reproduces Table 3 and Figure S3                                           #
########################################################################################    
# Necessary packages and functions
# library(GpGp)
# library(sp)
# library(FRK)
# library(autoFRK)
# library(LatticeKrig)
# library(support)
# library(SPlit)
# library(fields)
# library(assertthat)
# library(magrittr)
# library(purrr)
# library(rlist)
# library(geoR)
# library(ggplot2)
# library(foreach)
# library(doParallel)
# source(here("Rep_Points","SP-sample.R"))
# source(here("S2_VersusK","helper-functions.R"))



########    Reproduce Table 3      ##########################
# Real data 1: annual total precipitation anomalies data    #
# Compare the predictive performance (MSPE), computational  #
# time, number of bases of various low-rank approximation   #
# methods, including: SP(210), SP(500), SP(750), SP(1000),  #
# FRK(2), FRK(3), autoFRK(210), autoFRK(500), autoFRK(750), #
# autoFRK(1000), and LatticeKrig. Results based on full data#
# are also provided for reference.                          #
#############################################################
##############################################################################
### Load and illustrate the data (Figure S3)
load(here("Real_Data","anom1962.RData"))
Data=data.frame(lon=loc[,1],lat=loc[,2],z=z)
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
PN=ggplot(Data,aes(x=lon,y=lat,colour=z))+xlab("Longitude")+ylab("Latitude")+
  geom_point()+scale_color_gradientn(values = seq(0,1,0.125),colours = Col,name="Precipitation anomalies")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position = "right",
                   #legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PN


### Build the function for reproducing results of one replicate in Table 3
# Input: seed 
# Output: a matrix store the outputs of mspe, computational time, and bases number
Main=function(seed){
  ### Split the data into training and testing datasets
  # t1=proc.time()[[3]]
  set.seed(seed)
  it=sample(1:7352,7e3)
  X=loc[it,]
  Y=z[it]
  Xt=loc[-it,]
  Yt=z[-it]
  # t2=proc.time()[[3]]
  # t2-t1=0.016

  
  ### Evaluate the Gaussian process prediction \hat{f}_{\hat{c}} for reference, where
  ### \hat{c} is obtained by a software called ExaGeoStat. 
  t1=proc.time()[[3]]
  theta.exact=c(0.80560249,4.75562447,0.01,0.23764869)
  fullres.ExaGeoStat=hatf.hatc.rmse(X,Y,t(as.matrix(theta.exact)),Xt,Yt)
  t2=proc.time()[[3]]
  fullres.ExaGeoStat=c(fullres.ExaGeoStat^2,t2-t1,7e3)
  
  
  ### Evaluate the Gaussian process prediction \hat{f}_{\hat{c}} for reference, where
  ### \hat{c} is obtained by the Vecchia approximation in R package GpGp.
  t1=proc.time()[[3]]
  gpres=fit_model(Y,X,as.matrix(rep(1,nrow(X))),"matern_isotropic")
  theta.hat=c(gpres$covparms[1],theta=gpres$covparms[2],gpres$covparms[1]*gpres$covparms[4],gpres$covparms[3])
  t2=proc.time()[[3]]
  fullres.GpGp=hatf.hatc.rmse(X,Y,t(as.matrix(theta.hat)),Xt,Yt)
  t3=proc.time()[[3]]
  fullres.GpGp=c(fullres.GpGp^2,t3-t1,7e3)
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with SPs and \hat{c}
  fsp=function(k){
    # Select k SPs
    tsp1=proc.time()[[3]]
    spid=sp.sample(k,X,seed)
    tsp2=proc.time()[[3]]
    
    # Calculate using theta.hat
    tsp3=proc.time()[[3]]
    mspe.hat=tildef.hatc.rmse(X,Y,X[spid,],t(as.matrix(theta.hat)),Xt,Yt)
    tsp4=proc.time()[[3]]
    
    # Calculate using theta.exact
    tsp5=proc.time()[[3]]
    mspe.exact=tildef.hatc.rmse(X,Y,X[spid,],t(as.matrix(theta.exact)),Xt,Yt)
    tsp6=proc.time()[[3]]
    
    return(cbind(c(mspe.hat^2,tsp4-tsp1+t2-t1,k),c(mspe.exact^2,tsp6-tsp5+tsp2-tsp1+t2-t1,k)))
  }
  SPres=sapply(c(210,500,750,1000),fsp)
  
  
  ### Evaluate the FRK method with nres=2
  t4=proc.time()[[3]]
  Dfrk=data.frame(x1=X[,1],x2=X[,2],Y=Y)
  coordinates(Dfrk)=~x1+x2
  GridBAUs1=auto_BAUs(manifold = plane(),cellsize = c(1,0.5),type="grid",data=Dfrk,convex=-0.05,nonconvex_hull = FALSE)
  GridBAUs1$fs=1
  G=auto_basis(manifold = plane(),data=Dfrk,nres=2,type="exp",regular=1)
  f=Y~1
  S=SRE(f=f,data=list(Dfrk),BAUs=GridBAUs1,basis=G,est_error=TRUE,average_in_BAU = FALSE)
  S=SRE.fit(SRE_model=S,n_EM=3,tol=0.01)
  Dfrkt=data.frame(x1=Xt[,1],x2=Xt[,2],Y=Yt)
  coordinates(Dfrkt)=~x1+x2
  ythat=predict(S,newdata=Dfrkt)$mu
  t5=proc.time()[[3]]
  Frkres1=c(mean((Yt-ythat)^2),t5-t4,length(G@fn))
  
  
  ### Evaluate the FRK method with nres=3
  t6=proc.time()[[3]]
  Dfrk=data.frame(x1=X[,1],x2=X[,2],Y=Y)
  coordinates(Dfrk)=~x1+x2
  GridBAUs1=auto_BAUs(manifold = plane(),cellsize = c(1,0.5),type="grid",data=Dfrk,convex=-0.05,nonconvex_hull = FALSE)
  GridBAUs1$fs=1
  G=auto_basis(manifold = plane(),data=Dfrk,nres=3,type="exp",regular=1)
  f=Y~1
  S=SRE(f=f,data=list(Dfrk),BAUs=GridBAUs1,basis=G,est_error=TRUE,average_in_BAU = FALSE)
  S=SRE.fit(SRE_model=S,n_EM=3,tol=0.01)
  Dfrkt=data.frame(x1=Xt[,1],x2=Xt[,2],Y=Yt)
  coordinates(Dfrkt)=~x1+x2
  ythat=predict(S,newdata=Dfrkt)$mu
  t7=proc.time()[[3]]
  Frkres2=c(mean((Yt-ythat)^2),t7-t6,length(G@fn))
  
  
  ### Evaluate the autoFRK method
  fauto=function(k){
    tauto1=proc.time()[[3]]
    autores=autoFRK(Data=Y,loc=X,maxK=k)
    ythat=predict.FRK(autores,newloc=Xt)
    tauto2=proc.time()[[3]]
    mspe=mean((Yt-ythat$pred.value)^2)
    return(c(mspe,tauto2-tauto1,length(autores$w)))
  }
  Autores=sapply(c(210,500,750),fauto)
  
  
  ### Evaluate the LatticeKrig method
  t8=proc.time()[[3]]
  LKinfol=LKrigSetup(X,NC=15,nlevel=2,a.wght=6,nu=1.0)
  LKres=LatticeKrig(X,Y,LKinfo=LKinfol)
  ythat=predict(LKres,xnew=Xt)
  t9=proc.time()[[3]]
  Lkres=c(mean((Yt-ythat)^2),t9-t8,1258)
  
  return(cbind(fullres.GpGp,fullres.ExaGeoStat,SPres[1:3,],SPres[4:6,],Frkres1,Frkres2,Autores,Lkres))
}
# Do one replicate to save time
# t1=proc.time()[[3]]
result=Main(1)
# t2=proc.time()[[3]]
# t2-t1=808.371

# Repeat 100 times and calculate the average
result=sapply(1:100,Main(666*i+888))
# write.csv(result,here("Real_Data","Res_DP.csv"))
result=as.matrix(read.csv(here("Real_Data","Res_DP.csv"))[,-1])
res=matrix(apply(result,1,mean),3,16)

# > res
#      fullres.GpGp fullres.ExaGeoStat SP(210)   SP(500)     SP(750)     SP(1000)     SP.exa(210)
# [1,]    0.2209736    0.2209601   0.3169398   0.2789732   0.2613839    0.2515936   0.3171139
# [2,]  159.0795500  161.1304900   8.6683800   9.9063200  11.2450600   12.6980300   8.6683800
# [3,] 7000.0000000 7000.0000000 210.0000000 500.0000000 750.0000000 1000.0000000 210.0000000
# [,8]  SP.exa(500) SP.exa(750)  SP.exa(1000)  FRK(2)      FRK(3)   autoFRK(210)  autoFRK(500)
# [1,]   0.2787938   0.2610819    0.251244   0.2867042    0.2780921   0.3046087   0.268832
# [2,]   9.8850600  11.2134900   12.577280  15.1097500  382.9467800   5.7422000  13.003680
# [3,] 500.0000000 750.0000000 1000.000000 210.0000000 1936.6500000 205.0000000 499.490000
#      autoFRK(750)   LatticeKrig
# [1,]   0.2597549    0.2995753
# [2,]  21.1422800   51.0432000
# [3,] 747.6600000 1258.0000000
