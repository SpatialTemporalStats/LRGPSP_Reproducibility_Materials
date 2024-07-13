########################################################################################
# This file reproduces Figure 6                                                        #
######################################################################################## 
# Necessary packages and functions
# library(support)
# library(SPlit)
# library(fields)
# library(mvtnorm)
# library(assertthat)
# library(magrittr)
# library(purrr)
# library(rlist)
# library(ggplot2)
# library(foreach)
# library(doParallel)
# source(here("Rep_Points/k_DPP","helper-functions.R"))
# source(here("Rep_Points/k_DPP","kdpp-sample.R"))
# source(here("Rep_Points","SP-sample.R"))
# source(here("Rep_Points","Rand-sample.R"))
# source(here("Rep_Points","Grid-sample.R"))
# source(here("S2_VersusK","helper-functions.R"))



########    Reproduce Figure 6      #############################
# True covariance thetat=c(1.5,0.168639,0.27,1.5)               #
# Estimated covariance nu=1 (over) nu=3 (under)                 #
# Rep-points Support points and Random                          #
# n=5000                                                        #
# k= 1.5*n^(2/2.9), e.g., n=5000 -> k=534                       #
# tau^2=0.27*10^{-3,2,1,0,0.5,1,1.5,2}                          #
#################################################################
#################################################################################
### Build the function for reproducing results of one replicate in Figure 6
# Input: seed 
# Output: a matrix store the outputs  
Main=function(seed){
  ### Set parameters
  nt=n=5e3
  tauseq=0.27*10^c(-3,-2,-1,0,0.5,1,1.5,2)
  thetat=c(1.5,0.168639,0.27,1.5)
  k=1.5*round(n^(2/2.9))
  
  ### generate the training and test data
  set.seed(seed)
  loc=cbind(runif(n+nt),runif(n+nt))
  C=thetat[1]*stationary.cov(loc,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  y.t=rmvnorm(1,mean=rep(0,n+nt),sigma=C,method="chol")
  Xt=loc[1:nt,]
  Yt=y.t[1:nt]
  X=loc[(nt+1):(nt+n),]
  Y=y.t[(nt+1):(nt+n)]+rnorm(n,0,sqrt(thetat[3]))
  
  
  ### Select k SPs
  spid=sp.sample(k,X,seed)
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}}
  Thetaset=cbind(rep(thetat[1],length(tauseq)),rep(thetat[2],length(tauseq)),
                 tauseq,rep(thetat[4]*2/3,length(tauseq)))
  resnu1=tildef.hatc.rmse(X,Y,X[spid,],Thetaset,Xt,Yt)
  Thetaset=cbind(rep(thetat[1],length(tauseq)),rep(thetat[2],length(tauseq)),
                 tauseq,rep(thetat[4],length(tauseq)))
  resnu2=tildef.hatc.rmse(X,Y,X[spid,],Thetaset,Xt,Yt)
  Thetaset=cbind(rep(thetat[1],length(tauseq)),rep(thetat[2],length(tauseq)),
                 tauseq,rep(thetat[4]*2,length(tauseq)))
  resnu3=tildef.hatc.rmse(X,Y,X[spid,],Thetaset,Xt,Yt)
  out=cbind(resnu1^2,resnu2^2,resnu3^2)
  return(out)
}
# t1=proc.time()[[3]]
# a=Main(1)
# t2=proc.time()[[3]]
# t2-t1=67.107
cl<- makeCluster(4) 
registerDoParallel(cl) 
result= foreach(i=1:100,
                 .combine=cbind,
                 .packages=c("geoR","fields","support","SPlit","mvtnorm")) %dopar% Main(666*i+108)
stopCluster(cl)
# result=as.matrix(read.csv(here("S4_Tau","Tau.csv"))[,-1])
res=matrix(0,8,3)
for(i in 1:100){
  res=res+result[,(3*i-2):(3*i)]
}
res=res/100
dataT=data.frame(MSPE=c(res),
                 nu=as.factor(rep(c("1.0","1.5","3.0"),each=8)),
                 taur=rep(c(1e-3,1e-2,1e-1,1,sqrt(10),10,10^(1.5),1e2),times=3))
p=ggplot(data=dataT,aes(x=log10(taur),y=log(MSPE),colour=nu))+
  geom_point()+geom_line()+
  #geom_hline(yintercept = fmse,linetype=3)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(0,1),
                   legend.position =c(0,1),
                   legend.title = element_text(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_linetype_manual(values=c(1,2,3))+
  scale_colour_manual(values=c("#4DAF4A","#E41A1C","#984EA3"),
                      name=expression(nu))+
  xlab(" ")+ylab("log(MSPE)")
p
