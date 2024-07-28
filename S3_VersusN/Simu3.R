########################################################################################
# This file reproduces Figures 4, 5, and S1                                            #
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
# source(here("Rep_Points","SP-sample.R"))
# source(here("Rep_Points","Rand-sample.R"))
# source(here("S2_VersusK","helper-functions.R"))



########    Reproduce Figures 4 and 5      ######################
# True covariance thetat=c(1.5,0.168639,0.27,1.5)               #
# Estimated covariance nu=0.5 and 1 (over) nu=2,3 (under)       #
# Rep-points Support points and Random                          #
# n=1000, 1500, 2000, ..., 6000, 6500, 7000                     #
# k= 1.5*n^(2/2.9), e.g., n=5000 -> k=534                       #
#################################################################
#################################################################################
### For a fixed n, get the MSPE of \tilde{f}_{\hat{c}} with SPs and Rands of sufficient k.
# Input: n       size of the training data
#        seed    seed
# Output: a matrix store the outputs of MSPE
Main=function(n,seed){
  ### Set parameters
  rseq=c(1/3,2/3,1,4/3,2)
  thetat=c(1.5,0.168639,0.27,1.5)
  Thetaset=cbind(rep(thetat[1],5),rep(thetat[2],5),rep(thetat[3],5),thetat[4]*rseq) # Candidate of \hat{c}
  k=1.5*round(n^(2/2.9))
  
  
  ### Generate the training and test data
  set.seed(seed)
  nt=5e3
  loc=cbind(runif(n+nt),runif(n+nt))
  C=thetat[1]*stationary.cov(loc,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  y.t=rmvnorm(1,mean=rep(0,n+nt),sigma=C,method="chol")
  Xt=loc[1:nt,]
  Yt=y.t[1:nt]
  X=loc[(nt+1):(nt+n),]
  Y=y.t[(nt+1):(nt+n)]+rnorm(n,0,sqrt(thetat[3]))
  
  
  ### Select k SPs and Rands
  spid=sp.sample(k,X,seed)
  rid=rand.sample(k,X,seed)
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}}
  spres=tildef.hatc.rmse(X,Y,X[spid,],Thetaset,Xt,Yt)
  randres=tildef.hatc.rmse(X,Y,X[rid,],Thetaset,Xt,Yt)
  return(rbind(spres^2,randres^2))
}


### Function for reproducing results of one replicate
# Input: seed 
# Output: a matrix store the MSPEs for all ns and the negative slopes (NSs), see Section 3.3
nseq=seq(1e3,7e3,5e2)    # Candidates of n
Main2=function(seed){
  ### Calculate the MSPEs for all ns
  afuc=function(n){
    return(Main(n,seed))
  }
  out=sapply(nseq,afuc)    # (2*length(rseq)) by length(nseq) matrix, storing theoutputs
  
  ### Calculate the NSs
  bfuc=function(x){
    return(lm(log(x)~log(nseq))$coefficients[2])
  }
  res=apply(out,1,bfuc)
  return(cbind(out,-res))
}
# t1=proc.time()[[3]]
# a=Main2(1)
# t2=proc.time()[[3]]
# t2-t1=387.594/60=6.4599 minutes
cl<- makeCluster(4)
registerDoParallel(cl)
resultn= foreach(i=1:100,
                 .combine=cbind,
                 .packages=c("geoR","fields","support","SPlit","mvtnorm")) %dopar% Main2(666*i+108)
stopCluster(cl)
# resultn=as.matrix(read.csv(here("S3_VersusN","Setting1.csv"))[,-1]) 
B=matrix(0,10,100)
A=matrix(0,10,1300)
for(i in 1:100){
  B[,i]=resultn[,14*i]
  A[,(13*i-12):(13*i)]=resultn[,(14*i-13):(14*i-1)]
}


### Plot Figure 4
# Calculate the average of log(MSPE) among 100 replicates
A=log(A)
spu1=apply(matrix(A[1,],13,100),1,mean)
spu2=apply(matrix(A[3,],13,100),1,mean)
spt=apply(matrix(A[5,],13,100),1,mean)
spo1=apply(matrix(A[7,],13,100),1,mean)
spo2=apply(matrix(A[9,],13,100),1,mean)
ru1=apply(matrix(A[2,],13,100),1,mean)
ru2=apply(matrix(A[4,],13,100),1,mean)
rt=apply(matrix(A[6,],13,100),1,mean)
ro1=apply(matrix(A[8,],13,100),1,mean)
ro2=apply(matrix(A[10,],13,100),1,mean)
data15=data.frame(MSPE=c(spt,spu1,spu2,spo1,spo2,rt,ru1,ru2,ro1,ro2),
                  MSPE1=c(lm(spt~log(nseq))$fitted.values,lm(spu1~log(nseq))$fitted.values,lm(spu2~log(nseq))$fitted.values,
                          lm(spo1~log(nseq))$fitted.values,lm(spo2~log(nseq))$fitted.values,lm(rt~log(nseq))$fitted.values,
                          lm(ru1~log(nseq))$fitted.values,lm(ru2~log(nseq))$fitted.values,lm(ro1~log(nseq))$fitted.values,lm(ro2~log(nseq))$fitted.values),
                  n=rep(nseq,times=10),
                  Method=as.factor(rep(c("SP(1.5)","SP(0.5)","SP(1.0)","SP(2.0)","SP(3.0)","Rand(1.5)","Rand(0.5)","Rand(1.0)","Rand(2.0)","Rand(3.0)"),each=13)))
# write.csv(data15,here("S3_VersusN","Set1_MSPE.csv"))
p=ggplot(data=data15,aes(x=log(n),y=MSPE,colour=Method,shape=Method))+
  geom_point()+geom_line(mapping = aes(x=log(n),y=MSPE1,colour=Method,linetype=Method))+
  #geom_hline(yintercept = fmse,linetype=3)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position ="right",
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_linetype_manual(values=c(2,2,2,2,2,1,1,1,1,1))+
  scale_shape_manual(values=c(2,2,2,2,2,1,1,1,1,1))+
  scale_colour_manual(values=c("orange","#4DAF4A","#E41A1C","blue","#984EA3","orange","#4DAF4A","#E41A1C","blue","#984EA3"))+
  xlab("log(n)")+ylab("log(MSPE)")
p



### Plot Figure 5
dataG=data.frame(Ga=c(c(t(1/(1-B[1:6,]))),c(t(rep(1,2)%*%t(rep(2.9-1,100))/B[c(7,9),])),
                      c(t(rep(1,2)%*%t(rep(2.9-1,100))/B[c(8,10),]))),
                 Method=as.factor(rep(c("SP(0.5)","Rand(0.5)","SP(1.0)","Rand(1.0)","SP(1.5)","Rand(1.5)",
                                        "SP(2.0)","SP(3.0)","Rand(2.0)","Rand(3.0)"),each=100)))
# write.csv(dataG,here("S3_VersusN","Set1_GS.csv"))
p=ggplot(data = dataG,aes(x=Method,y=Ga,color=Method))+geom_boxplot()+
  stat_summary(fun.y=mean,size=0.2,shape=21)+
  geom_linerange(y=1.43,xmin=0.65,xmax=1.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.22,xmin=1.65,xmax=2.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.90,xmin=2.65,xmax=3.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=3.55,xmin=3.65,xmax=4.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=4.75,xmin=4.65,xmax=5.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=1.43,xmin=5.65,xmax=6.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.22,xmin=6.65,xmax=7.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.90,xmin=7.65,xmax=8.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=3.55,xmin=8.65,xmax=9.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=4.75,xmin=9.65,xmax=10.43,linetype=3,color="black",size=0.5)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position ="right",
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_color_manual(values=c("orange","#4DAF4A","#E41A1C","blue","#984EA3","orange","#4DAF4A","#E41A1C","blue","#984EA3"))+
  xlab(" ")+ylab(expression(gamma))
p




########    Reproduce Figure S1            ######################
# True covariance thetat=c(0.5,0.168639,0.27,1.5)               #
# Estimated covariance nu=0.5 and 1 (over) nu=2,3 (under)       #
# Rep-points Support points and Random                          #
# n=4000, 1500, 2000, ..., 6000, 6500, 10000                    #
# k= 1.5*n^(2/2.9), e.g., n=5000 -> k=534                       #
#################################################################
#################################################################################
### For a fixed n, get the MSPE of \tilde{f}_{\hat{c}} with SPs and Rands of sufficient k.
# Input: n       size of the training data
#        seed    seed
# Output: a matrix store the outputs of MSPE
Main=function(n,seed){
  ### Set parameters
  rseq=c(1/3,2/3,1,4/3,2)
  thetat=c(0.5,0.168639,0.27,1.5)
  Thetaset=cbind(rep(thetat[1],5),rep(thetat[2],5),rep(thetat[3],5),thetat[4]*rseq) # Candidate of \hat{c}
  k=1.5*round(n^(2/2.9))
  
  
  ### Generate the training and test data
  set.seed(seed)
  nt=5e3
  loc=cbind(runif(n+nt),runif(n+nt))
  C=thetat[1]*stationary.cov(loc,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  y.t=rmvnorm(1,mean=rep(0,n+nt),sigma=C,method="chol")
  Xt=loc[1:nt,]
  Yt=y.t[1:nt]
  X=loc[(nt+1):(nt+n),]
  Y=y.t[(nt+1):(nt+n)]+rnorm(n,0,sqrt(thetat[3]))
  
  
  ### Select k SPs and Rands
  spid=sp.sample(k,X,seed)
  rid=rand.sample(k,X,seed)
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}}
  spres=tildef.hatc.rmse(X,Y,X[spid,],Thetaset,Xt,Yt)
  randres=tildef.hatc.rmse(X,Y,X[rid,],Thetaset,Xt,Yt)
  return(rbind(spres^2,randres^2))
}



### Function for reproducing results of one replicate
# Input: seed 
# Output: a matrix store the MSPEs for all ns and the negative slopes (NSs), see Section 3.3
nseq=seq(4e3,1e4,5e2)    # Candidates of n
Main2=function(seed){
  ### Calculate the MSPEs for all ns
  afuc=function(n){
    return(Main(n,seed))
  }
  out=sapply(nseq,afuc)    # (2*length(rseq)) by length(nseq) matrix, storing theoutputs
  
  ### Calculate the NSs
  bfuc=function(x){
    return(lm(log(x)~log(nseq))$coefficients[2])
  }
  res=apply(out,1,bfuc)
  return(cbind(out,-res))
}
# t1=proc.time()[[3]]
# a=Main2(1)
# t2=proc.time()[[3]]
# t2-t1=714.145/60=11.9 minutes
cl<- makeCluster(4)
registerDoParallel(cl)
resultn= foreach(i=1:100,
                 .combine=cbind,
                 .packages=c("geoR","fields","support","SPlit","mvtnorm")) %dopar% Main2(666*i+108)
stopCluster(cl)
# write.csv(resultn,here("S3_VersusN","SettingS.csv"))
# resultn=as.matrix(read.csv(here("S3_VersusN","SettingS.csv"))[,-1])
B=matrix(0,10,100)
A=matrix(0,10,1300)
for(i in 1:100){
  B[,i]=resultn[,14*i]
  A[,(13*i-12):(13*i)]=resultn[,(14*i-13):(14*i-1)]
}


### Plot Figure S1
# Calculate the average of log(MSPE) among 100 replicates
A=log(A)
spu1=apply(matrix(A[1,],13,100),1,mean)
spu2=apply(matrix(A[3,],13,100),1,mean)
spt=apply(matrix(A[5,],13,100),1,mean)
spo1=apply(matrix(A[7,],13,100),1,mean)
spo2=apply(matrix(A[9,],13,100),1,mean)
ru1=apply(matrix(A[2,],13,100),1,mean)
ru2=apply(matrix(A[4,],13,100),1,mean)
rt=apply(matrix(A[6,],13,100),1,mean)
ro1=apply(matrix(A[8,],13,100),1,mean)
ro2=apply(matrix(A[10,],13,100),1,mean)
datas=data.frame(MSPE=c(spt,spu1,spu2,spo1,spo2,rt,ru1,ru2,ro1,ro2),
                 MSPE1=c(lm(spt~log(nseq))$fitted.values,lm(spu1~log(nseq))$fitted.values,lm(spu2~log(nseq))$fitted.values,
                         lm(spo1~log(nseq))$fitted.values,lm(spo2~log(nseq))$fitted.values,lm(rt~log(nseq))$fitted.values,
                         lm(ru1~log(nseq))$fitted.values,lm(ru2~log(nseq))$fitted.values,lm(ro1~log(nseq))$fitted.values,lm(ro2~log(nseq))$fitted.values),
                 n=rep(nseq,times=10),
                 Method=as.factor(rep(c("SP(1.5)","SP(0.5)","SP(1.0)","SP(2.0)","SP(3.0)","Rand(1.5)","Rand(0.5)","Rand(1.0)","Rand(2.0)","Rand(3.0)"),each=13)))
# write.csv(datas,here("S3_VersusN","SetS_MSPE.csv"))
p=ggplot(data=datas,aes(x=log(n),y=MSPE,colour=Method,shape=Method))+
  geom_point()+geom_line(mapping = aes(x=log(n),y=MSPE1,colour=Method,linetype=Method))+
  #geom_hline(yintercept = fmse,linetype=3)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position ="right",
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_linetype_manual(values=c(2,2,2,2,2,1,1,1,1,1))+
  scale_shape_manual(values=c(2,2,2,2,2,1,1,1,1,1))+
  scale_colour_manual(values=c("orange","#4DAF4A","#E41A1C","blue","#984EA3","orange","#4DAF4A","#E41A1C","blue","#984EA3"))+
  xlab("log(n)")+ylab("log(MSPE)")
p
dataG=data.frame(Ga=c(c(t(1/(1-B[1:6,]))),c(t(rep(1,2)%*%t(rep(2.9-1,100))/B[c(7,9),])),
                      c(t(rep(1,2)%*%t(rep(2.9-1,100))/B[c(8,10),]))),
                 Method=as.factor(rep(c("SP(0.5)","Rand(0.5)","SP(1.0)","Rand(1.0)","SP(1.5)","Rand(1.5)",
                                        "SP(2.0)","SP(3.0)","Rand(2.0)","Rand(3.0)"),each=100)))
# write.csv(dataG,here("S3_VersusN","SetS_GS.csv"))
p=ggplot(data = dataG,aes(x=Method,y=Ga,color=Method))+geom_boxplot()+
  stat_summary(fun.y=mean,size=0.2,shape=21)+
  geom_linerange(y=1.43,xmin=0.65,xmax=1.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.22,xmin=1.65,xmax=2.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.90,xmin=2.65,xmax=3.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=3.55,xmin=3.65,xmax=4.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=4.75,xmin=4.65,xmax=5.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=1.43,xmin=5.65,xmax=6.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.22,xmin=6.65,xmax=7.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=2.90,xmin=7.65,xmax=8.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=3.55,xmin=8.65,xmax=9.43,linetype=3,color="black",size=0.5)+
  geom_linerange(y=4.75,xmin=9.65,xmax=10.43,linetype=3,color="black",size=0.5)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position ="right",
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_color_manual(values=c("orange","#4DAF4A","#E41A1C","blue","#984EA3","orange","#4DAF4A","#E41A1C","blue","#984EA3"))+
  xlab(" ")+ylab(expression(gamma))
p







