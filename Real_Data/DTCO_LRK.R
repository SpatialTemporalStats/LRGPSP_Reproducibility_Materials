########################################################################################
# This file reproduces Table 4                                                         #
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



########    Reproduce Table 4      ##########################
# Real data 2: level-2 total column ozone data              #
# Compare the predictive performance (MSPE), computational  #
# time, number of bases of various low-rank approximation   #
# methods, including: SP(189), SP(500), SP(1361), SP(1755), #
# FRK(2), FRK(3), LatticeKrig(2), and LatticeKrig(3).       #
#############################################################
##############################################################################
### Load the data 
Ozone_dat=read.csv(here("Real_Data","Ozone_dat.csv"))
Data=data.frame(loc.x1=Ozone_dat$lon,loc.x2=Ozone_dat$lat,y=Ozone_dat$ozone)


### Estimate the covariance parameters with GpGp
t1=proc.time()[[3]]
gpres=fit_model(Data$y,cbind(Data$loc.x1,Data$loc.x2),as.matrix(rep(1,nrow(Data))),"matern_anisotropic2D")
t2=proc.time()[[3]]
MLEt=t2-t1
AL=matrix(gpres$covparms[c(2,3,3,4)],2,2)


### Build the function for reproducing results of one replicate in Table 4
# Input: seed 
# Output: a matrix store the outputs of mspe, computational time, and bases number
Main=function(seed){
  ### Split the data into training and testing datasets
  set.seed(seed)
  it=sample(1:nrow(Data),15e4)
  X=as.matrix(Data[it,-3])
  Y=Data[it,3]-mean(Data$y)
  Xt=as.matrix(Data[-it,-3])
  Yt=Data[-it,3]-mean(Data$y)
  
  ### Evaluate the LatticeKrig with nlevel=2
  t1=proc.time()[[3]]
  LKinfol=LKrigSetup(X,NC=15,nlevel=2,a.wght=4.1,nu=0.5)
  LKres=LatticeKrig(X,Y,LKinfo=LKinfol)
  ythat=predict(LKres,xnew=Xt)
  t2=proc.time()[[3]]
  Lkres1=c(mean((Yt-ythat)^2),t2-t1,LKinfol$latticeInfo$m)
  
  ### Evaluate the LatticeKrig with nlevel=3
  t1=proc.time()[[3]]
  LKinfol=LKrigSetup(X,NC=15,nlevel=3,a.wght=4.1,nu=0.5)
  LKres=LatticeKrig(X,Y,LKinfo=LKinfol)
  ythat=predict(LKres,xnew=Xt)
  t2=proc.time()[[3]]
  Lkres2=c(mean((Yt-ythat)^2),t2-t1,LKinfol$latticeInfo$m)
  
  ### Evaluate the FRK with nres=2
  t1=proc.time()[[3]]
  Dfrk=data.frame(x1=X[,1],x2=X[,2],Y=Y)
  coordinates(Dfrk)=~x1+x2
  GridBAUs1=auto_BAUs(manifold = plane(),cellsize = c(1,0.5),type="grid",data=Dfrk,convex=-0.05,nonconvex_hull = FALSE)
  GridBAUs1$fs=1
  G=auto_basis(manifold = plane(),data=Dfrk,nres=2,type="Gaussian",regular=1)
  f=Y~1
  S=SRE(f=f,data=list(Dfrk),BAUs=GridBAUs1,basis=G,est_error=TRUE,average_in_BAU = FALSE)
  S=SRE.fit(SRE_model=S,n_EM=3,tol=0.01)
  Dfrkt=data.frame(x1=Xt[,1],x2=Xt[,2],Y=Yt)
  coordinates(Dfrkt)=~x1+x2
  ythat=predict(S,newdata=Dfrkt)$mu
  t2=proc.time()[[3]]
  Frkres1=c(mean((Yt-ythat)^2),t2-t1,189)
  
  ### Evaluate the FRK with nres=3
  t1=proc.time()[[3]]
  Dfrk=data.frame(x1=X[,1],x2=X[,2],Y=Y)
  coordinates(Dfrk)=~x1+x2
  GridBAUs1=auto_BAUs(manifold = plane(),cellsize = c(1,0.5),type="grid",data=Dfrk,convex=-0.05,nonconvex_hull = FALSE)
  GridBAUs1$fs=1
  G=auto_basis(manifold = plane(),data=Dfrk,nres=3,type="Gaussian",regular=1)
  f=Y~1
  S=SRE(f=f,data=list(Dfrk),BAUs=GridBAUs1,basis=G,est_error=TRUE,average_in_BAU = FALSE)
  S=SRE.fit(SRE_model=S,n_EM=3,tol=0.01)
  Dfrkt=data.frame(x1=Xt[,1],x2=Xt[,2],Y=Yt)
  coordinates(Dfrkt)=~x1+x2
  ythat=predict(S,newdata=Dfrkt)$mu
  t2=proc.time()[[3]]
  Frkres2=c(mean((Yt-ythat)^2),t2-t1,1755)
  
  # prediction based on SP
  fsp=function(k){
    tsp1=proc.time()[[3]]
    set.seed(seed)
    spres=sp(k,2,dist.samp=X)
    spX=spres$sp
    Ch=(gpres$covparms[1])*stationary.cov(X%*%t(AL),spX%*%t(AL),Covariance = "Matern",smoothness=gpres$covparms[5],theta=1)
    Cth=(gpres$covparms[1])*stationary.cov(Xt%*%t(AL),spX%*%t(AL),Covariance = "Matern",smoothness=gpres$covparms[5],theta=1)
    Ckk=(gpres$covparms[1])*stationary.cov(spX%*%t(AL),Covariance = "Matern",smoothness=gpres$covparms[5],theta=1)
    cres=svd(t(Ch)%*%Ch+(gpres$covparms[1]*gpres$covparms[6])*Ckk)
    ythat=gpres$betahat+Cth%*%cres$u%*%diag(1/cres$d)%*%t(cres$v)%*%t(Ch)%*%(Y-gpres$betahat)
    tsp2=proc.time()[[3]]
    mspe=mean((Yt-ythat)^2)
    return(c(mspe,tsp2-tsp1+MLEt,k))
  }
  SPres=sapply(c(189,500,1361,1755),fsp)
  
  return(cbind(Lkres1,Lkres2,SPres,Frkres1,Frkres2))
}
# Do one replicate to save time
# t1=proc.time()[[3]]
result=Main(1)
# t2=proc.time()[[3]]
# t2-t1+MLEt=8098.03/3600=2.25 hours

# cl<- makeCluster(4) 
# registerDoParallel(cl) 
# result= foreach(i=1:100,
#                 .combine=cbind,
#                 .packages=c("geoR","SPlit","support","LatticeKrig","fields","FRK","sp")) %dopar% Main(666*i+888)
# stopCluster(cl)


