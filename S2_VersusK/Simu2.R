########################################################################################
# This file reproduces Figure 3 and Table 2                                            #
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



########    Reproduce Figure 3(a)      ##########################
# True covariance thetat=c(1.5,0.168639,0.27,1.5)               #
# Estimated covariance nu=0.5 and 1 (over) nu=2,3 (under)       #
# Rep-points: SPs, Grids, Rands, and k-DPPs                     #
# k= (0.10, 0.18, 0.28, 0.40, 0.55, 0.81, 1.12, 1.36)*n^(2/2.9) #
#      36 ,  64 , 100 , 144,  196 , 289 , 400 , 484             #
#################################################################
#################################################################################
### Build the function for reproducing results of one replicate in Figure 3(a)
# Input: seed 
# Output: a matrix store the outputs of energy distance and rmse 
Main=function(seed){
  ### Generate the training and test data
  # t1=proc.time()[[3]]
  set.seed(seed)
  n=nt=5e3
  thetat=c(1.5,0.168639,0.27,1.5)
  loc=cbind(runif(n+nt),runif(n+nt))
  C=thetat[1]*stationary.cov(loc,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  y.t=rmvnorm(1,mean=rep(0,n+nt),sigma=C,method="chol")
  Xt=loc[1:nt,]
  Yt=y.t[1:nt]
  X=loc[(nt+1):(nt+n),]
  Y=y.t[(nt+1):(nt+n)]+rnorm(n,0,sqrt(thetat[3]))
  # t2=proc.time()[[3]]
  # t2-t1=16.884
  
  
  ### Evaluate the Gaussian process prediction \hat{f}_{\hat{c}} for reference
  # t1=proc.time()[[3]]
  rseq=c(1/3,2/3,1,4/3,2)
  Thetaset=cbind(rep(thetat[1],5),rep(thetat[2],5),rep(thetat[3],5),thetat[4]*rseq) # Candidate of \hat{c}
  fullres=hatf.hatc.rmse(X,Y,Thetaset,Xt,Yt)
  # t2=proc.time()[[3]]
  # t2-t1=321.339
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b SPs and \hat{c}
  # Input: size of SPs
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.sp=function(b){
    # Select b SPs
    idsp=sp.sample(b,X,seed)
    # Caculate their energy distance
    Edist=energy.dist(X,idsp)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    spres=tildef.hatc.rmse(X,Y,X[idsp,],Thetaset,Xt,Yt)
    return(c(Edist,spres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with SPs for different b values
  bseq=c(36,64,100,144,196,289,400,484)    # Candidate of sizes of rep-points
  # t1=proc.time()[[3]]
  res.sp=sapply(bseq,tildef.sp)            # (1+nrow(Thetaset)) by length(bseq) matrix, 
  # t2=proc.time()[[3]]                      # storing energy distances and rmses for all bs.
  # t2-t1=43.199
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b Grids and \hat{c}
  # Input: size of Grids
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.grid=function(b){
    # Select b Grids
    idgrid=grid.sample(b,X)
    # Caculate their energy distance
    Edist=energy.dist(X,idgrid)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    gridres=tildef.hatc.rmse(X,Y,X[idgrid,],Thetaset,Xt,Yt)
    return(c(Edist,gridres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with Grids for different b values
  # t1=proc.time()[[3]]
  res.grid=sapply(bseq,tildef.grid)
  # t2=proc.time()[[3]]
  # t2-t1=31.229
  

  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b Rands and \hat{c}
  # Input: size of Rands
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.rand=function(b){
    # Select b Rands
    idrand=rand.sample(b,X,seed)
    # Caculate their energy distance
    Edist=energy.dist(X,idrand)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    randres=tildef.hatc.rmse(X,Y,X[idrand,],Thetaset,Xt,Yt)
    return(c(Edist,randres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with Rands for different b values
  # t1=proc.time()[[3]]
  res.rand=sapply(bseq,tildef.rand)
  # t2=proc.time()[[3]]
  # t2-t1=29.445
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with a given \hat{c} and various bs
  # Input: ratio of nu
  # Output: (2*length(bseq)) by length(rseq) matrix storing storing energy distances and rmses for all bs and rs.
  tildef.dpp=function(r){
    # The method of k-DPP can not be available for all b values. 
    # Therefore, the first step is to test which b values are okay for k-DPP method and which ones are not.
    C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=r*thetat[4],theta=thetat[2])
    L.decomposed=decompose.matrix(C)
    lambda=L.decomposed$lambda
    n.lambda=length(lambda)
    # There will be numerical values for available b and NaNs for not available ones.
    b.can=rep(0,length(bseq))   
    for(i in 1:length(bseq)){
      E=compute.elementary.symmetric.polynomial(lambda,bseq[i])
      b.can[i]=lambda[n.lambda]*E[bseq[i],n.lambda]/E[bseq[i]+1,n.lambda+1]
    }
    
    # For b values which are okay for k-DPP, evaluate \tilde{f}_{\hat{c}} with b k-DPPs and the given \hat{c}.
    # Input: id of b value
    # Output: a 2 by 1 vector of energy distance and the rmses of \tilde{f}_{\hat{c}} for the given \hat{c}.
    sres=function(ib){
      if(!is.na(b.can[ib])){
        # Select b k-DPPs 
        set.seed(seed)
        iddpp=k.dpp.sample(L.decomposed,bseq[ib])
        # Caculate their energy distance
        Edist=energy.dist(X,iddpp)
        # Evaluate the \tilde{f}_{\hat{c}} based on them 
        dpp.res=tildef.hatc.rmse(X,Y,X[iddpp,],t(as.matrix(thetat*c(1,1,1,r))),Xt,Yt)
        return(c(Edist,dpp.res))
      }
      else{
        return(c(0,0))
      }
    }
    # t1=proc.time()[[3]]
    dppRes=sapply(1:length(bseq),sres)     # 2 by length(bseq) matrix, 
                                           # storing energy distances and rmses for all bs. 
    # t2=proc.time()[[3]]
    # t2-t1 is about 256.214 seconds for r=1/3 and 2, 782.278 second for r=4/3, 4640.816 seconds for r=2/3 and 1
    return(dppRes)
  }
  res.dpp=sapply(rseq,tildef.dpp)          # (2*length(bseq)) by length(rseq) matrix. Put c(dppRes) in each colume.
                                           # The total computational time for k-DPP is about 10576.34/3600=2.938 hours 
  main.res=cbind(res.sp,res.grid,res.rand,c(0,fullres),rbind(res.dpp[c(1:8)*2-1,3],t(res.dpp[c(1:8)*2,])),
                 rbind(rep(0,8),t(res.dpp[c(1:8)*2-1,])))    # resort the outputs
  return(main.res)                         # Total computational is about 11018.44/3600=3.0607 hours.
}


### Three options for reproducibility
# Option 1. Do one replicate to save time
result=Main(1)
res=result

# Option 2. Do 100 replicates as what we did 
# cl<- makeCluster(4) 
# registerDoParallel(cl) 
# result=foreach(i=1:100,
#                .combine=cbind,
#                .packages=c("geoR","fields","support","SPlit","mvtnorm",
#                             "assertthat","magrittr","purrr","rlist")) %dopar% Main(666*i+108)
# stopCluster(cl)
# res=matrix(0,6,41)
# for(i in 1:100){
#   res=res+result[,(41*i-40):(41*i)]
# }
# res=res/100

# Option 3. Load the outputs obtained by us
# result=as.matrix(read.csv(here("S2_VersusK","Setting1.csv"))[,-1])
# res=matrix(0,6,41)
# for(i in 1:100){
#   res=res+result[,(41*i-40):(41*i)]
# }
# res=res/100


### Plot Figure 3(a)
# Check available b values for the method k-DPP in each replicate
# dpk=matrix(0,5,100)
# for(i in 1:100){
#   aa=result[,(41*i-40):(41*i)]
#   for(j in 1:5){
#     dpk[j,i]=which.max(which(aa[j+1,34:41]>0))
#   }
# }
# min(dpk)
Set1Sum=data.frame(Edist=c(res[1,1:24],res[1,26:33]),
                   RMSPE=c(t(res[-1,c(1:24,26:33)])),
                   Method=as.factor(rep(rep(c("SP","Grid","Rand","DPP"),each=8),times=5)),
                   nu=as.factor(rep(c("0.5","1.0","1.5","2.0","3.0"),each=32)),
                   k=rep(bseq,times=20))
# write.csv(Set1Sum,here("S2_VersusK","Set1Sum.csv"))
id1=which(Set1Sum$nu=="0.5")
id2=which(Set1Sum$nu=="2.0")
id3=which(Set1Sum$RMSPE==0)
Set1Sum$RMSPE[id3]=rep(NA,length(id3))
p=ggplot(data=Set1Sum[-c(id1,id2),],aes(x=k,y=RMSPE,colour=nu,linetype=Method))+
  geom_point()+geom_line()+
  geom_hline(yintercept=res[3,25],linetype=4,colour="#4DAF4A")+
  geom_hline(yintercept=res[4,25],linetype=4,colour="#E41A1C")+
  geom_hline(yintercept=res[6,25],linetype=4,colour="#984EA3")+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position =c(1,1),
                   #legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_colour_manual(values=c("#4DAF4A","#E41A1C","#984EA3"),name=expression(nu))+
  scale_linetype_manual(values=c(5,2,3,1),name="Rep-points")+
  #scale_y_continuous(limits = c(0.1,0.3),breaks = seq(0.1,0.3,0.05))+
  xlab("k")+ylab("RMSPE")
p




#######   Reproduce Figure 3(b)       ###########################
# True covariance thetat=c(1.5,0.063240,0.27,1.5)               #
# Estimated covariance nu=0.5 and 1 (over) nu=2,3 (under)       #
# Rep-points: SPs, Grids, Rands, and k-DPPs                     #
# k= (0.66, 1.21, 1.66, 2.05, 2.63, 3.12, 3.64, 4.02)*n^(2/2.8) #
#     289 ,  529, 729 , 900 , 1156, 1369, 1600, 1764            #
#################################################################
#################################################################################
### Build the function for reproducing results of one replicate in Figure 3(b)
# Input: seed 
# Output: a matrix store the outputs of energy distance and rmse 
Main=function(seed){
  # Generate the training and test data
  # t1=proc.time()[[3]]
  set.seed(seed)
  n=nt=5e3
  thetat=c(1.5,0.063240,0.27,1.5)
  loc=cbind(runif(n+nt),runif(n+nt))
  C=thetat[1]*stationary.cov(loc,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  y.t=rmvnorm(1,mean=rep(0,n+nt),sigma=C,method="chol")
  Xt=loc[1:nt,]
  Yt=y.t[1:nt]
  X=loc[(nt+1):(nt+n),]
  Y=y.t[(nt+1):(nt+n)]+rnorm(n,0,sqrt(thetat[3]))
  # t2=proc.time()[[3]]
  # t2-t1=18.993
  
  
  ### Evaluate the Gaussian process prediction \hat{f}_{\hat{c}} for reference
  # t1=proc.time()[[3]]
  rseq=c(1/3,2/3,1,4/3,2)
  Thetaset=cbind(rep(thetat[1],5),rep(thetat[2],5),rep(thetat[3],5),thetat[4]*rseq) # Candidate of \hat{c}
  fullres=hatf.hatc.rmse(X,Y,Thetaset,Xt,Yt)
  # t2=proc.time()[[3]]
  # t2-t1=360.715
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b SPs and \hat{c}
  # Input: size of SPs
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.sp=function(b){
    # Select b SPs
    idsp=sp.sample(b,X,seed)
    # Caculate their energy distance
    Edist=energy.dist(X,idsp)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    spres=tildef.hatc.rmse(X,Y,X[idsp,],Thetaset,Xt,Yt)
    return(c(Edist,spres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with SPs for different b values
  bseq=c(289,529,729,900,1156,1369,1600,1764)    # Candidate of sizes of rep-points
  # t1=proc.time()[[3]]
  res.sp=sapply(bseq,tildef.sp)              # (1+nrow(Thetaset)) by length(bseq) matrix, 
  # t2=proc.time()[[3]]                      # storing energy distances and rmses for all bs.
  # t2-t1=179.104
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b Grids and \hat{c}
  # Input: size of Grids
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.grid=function(b){
    # Select b Grids
    idgrid=grid.sample(b,X)
    # Caculate their energy distance
    Edist=energy.dist(X,idgrid)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    gridres=tildef.hatc.rmse(X,Y,X[idgrid,],Thetaset,Xt,Yt)
    return(c(Edist,gridres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with Grids for different b values
  # t1=proc.time()[[3]]
  res.grid=sapply(bseq,tildef.grid)
  # t2=proc.time()[[3]]
  # t2-t1=145.738
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b Rands and \hat{c}
  # Input: size of Rands
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.rand=function(b){
    # Select b Rands
    idrand=rand.sample(b,X,seed)
    # Caculate their energy distance
    Edist=energy.dist(X,idrand)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    randres=tildef.hatc.rmse(X,Y,X[idrand,],Thetaset,Xt,Yt)
    return(c(Edist,randres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with Rands for different b values
  # t1=proc.time()[[3]]
  res.rand=sapply(bseq,tildef.rand)
  # t2=proc.time()[[3]]
  # t2-t1=133.285
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with a given \hat{c} and various bs
  # Input: ratio of nu
  # Output: (2*length(bseq)) by length(rseq) matrix storing storing energy distances and rmses for all bs and rs.
  tildef.dpp=function(r){
    # The method of k-DPP can not be available for all b values. 
    # Therefore, the first step is to test which b values are okay for k-DPP method and which ones are not.
    # t1=proc.time()[[3]]
    C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=r*thetat[4],theta=thetat[2])
    L.decomposed=decompose.matrix(C)
    lambda=L.decomposed$lambda
    n.lambda=length(lambda)
    # There will be numerical values for available b and NaNs for not available ones.
    b.can=rep(0,length(bseq))   
    for(i in 1:length(bseq)){
      E=compute.elementary.symmetric.polynomial(lambda,bseq[i])
      b.can[i]=lambda[n.lambda]*E[bseq[i],n.lambda]/E[bseq[i]+1,n.lambda+1]
    }
    # t2=proc.time()[[3]]
    # t2-t1=34.452
    
    # For b values which are okay for k-DPP, evaluate \tilde{f}_{\hat{c}} with b k-DPPs and the given \hat{c}.
    # Input: id of b value
    # Output: a 2 by 1 vector of energy distance and the rmses of \tilde{f}_{\hat{c}} for the given \hat{c}.
    sres=function(ib){
      if(!is.na(b.can[ib])){
        # Select b k-DPPs 
        set.seed(seed)
        iddpp=k.dpp.sample(L.decomposed,bseq[ib])
        # Caculate their energy distance
        Edist=energy.dist(X,iddpp)
        # Evaluate the \tilde{f}_{\hat{c}} based on them 
        dpp.res=tildef.hatc.rmse(X,Y,X[iddpp,],t(as.matrix(thetat*c(1,1,1,r))),Xt,Yt)
        return(c(Edist,dpp.res))
      }
      else{
        return(c(0,0))
      }
    }
    dppRes=sapply(1:length(bseq),sres)     # 2 by length(bseq) matrix, 
                                           # storing energy distances and rmses for all bs. 
    return(dppRes)
  }
  res.dpp=sapply(rseq,tildef.dpp)          # (2*length(bseq)) by length(rseq) matrix. Put c(dppRes) in each colume.
                                           # The computtaional time is about 31462.95/3600=8.7397 hours
                              
  main.res=cbind(res.sp,res.grid,res.rand,c(0,fullres),rbind(res.dpp[c(1:8)*2-1,3],t(res.dpp[c(1:8)*2,])),
                 rbind(rep(0,8),t(res.dpp[c(1:8)*2-1,])))    # resort the outputs
  return(main.res)                         # Total computational is about 32300.78/3600=8.97 hours.
}


### Three options for reproducibility
# Option 1. Do one replicate to save time
result=Main(1)
res=result

# Option 2. Do 100 replicates as what we did 
# cl<- makeCluster(4) 
# registerDoParallel(cl) 
# result=foreach(i=1:100,
#                .combine=cbind,
#                .packages=c("geoR","fields","support","SPlit","mvtnorm",
#                             "assertthat","magrittr","purrr","rlist")) %dopar% Main(666*i+108)
# stopCluster(cl)
# res=matrix(0,6,41)
# for(i in 1:100){
#   res=res+result[,(41*i-40):(41*i)]
# }
# res=res/100

# Option 3. Load the outputs obtained by us
# result=as.matrix(read.csv(here("S2_VersusK","Setting2.csv"))[,-1])
# res=matrix(0,6,41)
# for(i in 1:100){
#   res=res+result[,(41*i-40):(41*i)]
# }
# res=res/100


### Plot Figure 3(b)
# Check available b values for the method k-DPP in each replicate
# dpk=matrix(0,5,100)
# for(i in 1:100){
#   aa=result[,(41*i-40):(41*i)]
#   for(j in 1:5){
#     if(sum(aa[j+1,34:41])>0){
#       dpk[j,i]=which.max(which(aa[j+1,34:41]>0))
#     }
#     else{
#       dpk[j,i]=NA
#     }
#   }
# }
Set2Sum=data.frame(Edist=c(res[1,1:24],res[1,26:33]),
                   RMSPE=c(t(res[-1,c(1:24,26:33)])),
                   Method=as.factor(rep(rep(c("SP","Grid","Rand","DPP"),each=8),times=5)),
                   nu=as.factor(rep(c("0.5","1.0","1.5","2.0","3.0"),each=32)),
                   k=rep(bseq,times=20))
# write.csv(Set2Sum,here("S2_VersusK","Set2Sum.csv"))
id1=which(Set2Sum$nu=="0.5")
id2=which(Set2Sum$nu=="2.0")
id3=which(Set2Sum$RMSPE==0)
Set2Sum$RMSPE[id3]=rep(NA,length(id3))
p=ggplot(data=Set2Sum[-c(id1,id2),],aes(x=k,y=RMSPE,colour=nu,linetype=Method))+
  geom_point()+geom_line()+
  geom_hline(yintercept=res[3,25],linetype=4,colour="#4DAF4A")+
  geom_hline(yintercept=res[4,25],linetype=4,colour="#E41A1C")+
  geom_hline(yintercept=res[6,25],linetype=4,colour="#984EA3")+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position ="none",
                   #legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_colour_manual(values=c("#4DAF4A","#E41A1C","#984EA3"),name=expression(nu))+
  scale_linetype_manual(values=c(5,2,3,1),name="Rep-points")+
  #scale_y_continuous(limits = c(0.1,0.3),breaks = seq(0.1,0.3,0.05))+
  xlab("k")+ylab("RMSPE")
p





#######   Reproduce Figure 3(c)         ##################################
# True covariance thetat=c(1.5,0.168639,0.27,1.5)                        #
# Consistent parameters c(1,0.147320,0.27,1.5) and c(2,0.185611,0.27,1.5)#
# Rep-points: SPs, Grids, Rands, and k-DPPs                              #
# k= (0.10, 0.18, 0.28, 0.40, 0.55, 0.81, 1.12, 1.36)*n^(2/2.9)          #
#      36 ,  64 , 100 , 144,  196 , 289 , 400 , 484                      #
##########################################################################
###############################################################################
thetat=c(1.5,0.168639,0.27,1.5)
theta1=c(1,0.147320,0.27,1.5) 
theta2=c(2,0.185611,0.27,1.5) 
Thetaset=rbind(thetat,theta1,theta2)


### Build the function for reproducing results of one replicate in Figure 3(c)
# Input: seed 
# Output: a matrix store the outputs of energy distance and rmse 
Main=function(seed){
  ### Generate the training and test data
  # t1=proc.time()[[3]]
  set.seed(seed)
  n=nt=5e3
  loc=cbind(runif(n+nt),runif(n+nt))
  C=thetat[1]*stationary.cov(loc,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  y.t=rmvnorm(1,mean=rep(0,n+nt),sigma=C,method="chol")
  Xt=loc[1:nt,]
  Yt=y.t[1:nt]
  X=loc[(nt+1):(nt+n),]
  Y=y.t[(nt+1):(nt+n)]+rnorm(n,0,sqrt(thetat[3]))
  # t2=proc.time()[[3]]
  # t2-t1=18.155
  
  
  ### Evaluate the Gaussian process prediction \hat{f}_{\hat{c}} for reference
  # t1=proc.time()[[3]]
  fullres=hatf.hatc.rmse(X,Y,Thetaset,Xt,Yt)
  # t2=proc.time()[[3]]
  # t2-t1=168.291
  
  
  ### Calculate some values for k-DPPs in advance. It could save time.
  # t1=proc.time()[[3]]
  C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  L.decomposed=decompose.matrix(C)
  C1=theta1[1]*stationary.cov(X,Covariance = "Matern",smoothness=theta1[4],theta=theta1[2])
  L1.decomposed=decompose.matrix(C1)
  C2=theta2[1]*stationary.cov(X,Covariance = "Matern",smoothness=theta2[4],theta=theta2[2])
  L2.decomposed=decompose.matrix(C2)
  # t2=proc.time()[[3]]
  # t2-t1=58.538
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b SPs and \hat{c}
  # Input: size of SPs
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.sp=function(b){
    # Select b SPs
    idsp=sp.sample(b,X,seed)
    # Caculate their energy distance
    Edist=energy.dist(X,idsp)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    spres=tildef.hatc.rmse(X,Y,X[idsp,],Thetaset,Xt,Yt)
    return(c(Edist,spres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with SPs for different b values
  bseq=c(36,64,100,144,196,289,400,484)    # Candidate of sizes of rep-points
  # t1=proc.time()[[3]]
  res.sp=sapply(bseq,tildef.sp)            # (1+nrow(Thetaset)) by length(bseq) matrix, 
  # t2=proc.time()[[3]]                   # storing energy distances and rmses for all bs.
  # t2-t1=19.501
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b Rands and \hat{c}
  # Input: size of Rands
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.rand=function(b){
    # Select b Rands
    idrand=rand.sample(b,X,seed)
    # Caculate their energy distance
    Edist=energy.dist(X,idrand)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    randres=tildef.hatc.rmse(X,Y,X[idrand,],Thetaset,Xt,Yt)
    return(c(Edist,randres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with Rands for different b values
  # t1=proc.time()[[3]]
  res.rand=sapply(bseq,tildef.rand)
  # t2=proc.time()[[3]]
  # t2-t1=9.544
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b Grids and \hat{c}
  # Input: size of Grids
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.grid=function(b){
    # Select b Grids
    idgrid=grid.sample(b,X)
    # Caculate their energy distance
    Edist=energy.dist(X,idgrid)
    # Evaluate the \tilde{f}_{\hat{c}} based on them 
    gridres=tildef.hatc.rmse(X,Y,X[idgrid,],Thetaset,Xt,Yt)
    return(c(Edist,gridres))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with Grids for different b values
  # t1=proc.time()[[3]]
  res.grid=sapply(bseq,tildef.grid)
  # t2=proc.time()[[3]]
  # t2-t1=8.285
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b k-DPPs and \hat{c}
  # Input: size of k-DPPs
  # Output: a vector of energy distance and the rmses of \tilde{f}_{\hat{c}}
  tildef.dpp=function(b){
    # For thetat
    set.seed(seed)
    iddpp=k.dpp.sample(L.decomposed,b)
    Edist=energy.dist(X,iddpp)
    dppres=tildef.hatc.rmse(X,Y,X[iddpp,],t(as.matrix(thetat)),Xt,Yt)
    # For theta1
    set.seed(seed)
    iddpp1=k.dpp.sample(L1.decomposed,b)
    dppres1=tildef.hatc.rmse(X,Y,X[iddpp1,],t(as.matrix(theta1)),Xt,Yt)
    # For thetat2
    set.seed(seed)
    iddpp2=k.dpp.sample(L2.decomposed,b)
    dppres2=tildef.hatc.rmse(X,Y,X[iddpp2,],t(as.matrix(theta2)),Xt,Yt)
    return(c(Edist,dppres,dppres1,dppres2))
  }
  ### Evaluate \tilde{f}_{\hat{c}} with k-DPPs for different b values
  # t1=proc.time()[[3]]
  res.dpp=sapply(bseq,tildef.dpp)
  # t2=proc.time()[[3]]
  # t2-t1=10904.11/3600=3.03 hours
  
  
  ### Resort the outputs
  main.res=cbind(res.sp,res.rand,res.grid,res.dpp,c(0,fullres))
  # Total computational is about 11186.42/3600=3.11 hours.
  return(main.res)
}




### Three options for reproducibility
# Option 1. Do one replicate to save time
result=Main(1)
res=result

# Option 2. Do 100 replicates as what we did 
# cl<- makeCluster(4) 
# registerDoParallel(cl) 
# result=foreach(i=1:100,
#                .combine=cbind,
#                .packages=c("geoR","fields","support","SPlit","mvtnorm",
#                             "assertthat","magrittr","purrr","rlist")) %dopar% Main(666*i+108)
# stopCluster(cl)
# res=matrix(0,4,33)
# for(i in 1:100){
#   res=res+result[,(33*i-32):(33*i)]
# }
# res=res/100


# Option 3. Load the outputs obtained by us
# result=as.matrix(read.csv(here("S2_VersusK","Setting3.csv"))[,-1])
# res=matrix(0,4,33)
# for(i in 1:100){
#   res=res+result[,(33*i-32):(33*i)]
# }
# res=res/100


### Plot Figure 3(c)
Set3Sum=data.frame(Edist=c(res[1,-33],res[1,-33],res[1,-33]),
                   RMSPE=c(t(res[-1,-33])),
                   Method=as.factor(rep(rep(c("SP","Grid","Rand","DPP"),each=8),times=3)),
                   theta=as.factor(rep(c("TP","CP1","CP2"),each=32)),
                   k=rep(bseq,times=12))
# write.csv(Set3Sum,here("S2_VersusK","Set3Sum.csv"))
p=ggplot(data=Set3Sum,aes(x=k,y=RMSPE,colour=theta,linetype=Method))+
  geom_point()+geom_line()+
  geom_hline(yintercept=res[2,33],linetype=4,colour="#4DAF4A")+
  geom_hline(yintercept=res[3,33],linetype=4,colour="#984EA3")+
  geom_hline(yintercept=res[4,33],linetype=4,colour="#E41A1C")+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position =c(1,1),
                   #legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_y_continuous(limits = c(0.1,0.3),breaks = seq(0.1,0.3,0.05))+
  scale_linetype_manual(values=c(5,2,3,1),name="Rep-points")+
  scale_colour_manual(values=c("#4DAF4A","#984EA3","#E41A1C"),name="Settings")+
  xlab("k")+ylab("RMSPE")
p



########  Reproduce Figure 3(d) and Table 2    #####################
# 75% data are in [0,0.5]^2, and 25% data are in [0,1]^2\[0,0.5]^2 #
# True covariance thetat=c(1.5,0.168639,0.27,1.5)                  #
# Rep-points: SPs, SP_Us, Grid, Rands, k-DPPs                      #
# k= (0.10, 0.18, 0.28, 0.40, 0.55, 0.81, 1.12, 1.36)*n^(2/2.9)    #
#      36 ,  64 , 100 , 144,  196 , 289 , 400 , 484                #
####################################################################
#######################################################################
### Function for selecting various kinds of rep-points and then calculate their energy distance and rmse
# Name: fsub
# Purpose: Select various kinds of rep-points of size b, calculate their energy distance and rmse
# Inputs: X             n by p training location matrix, n and p are the size and dimension of the training locations, respectively
#         Y             n by 1 training measurements
#         thetat        parameters in covariance
#         L.decomposed  used in k-DPP
#         Xt            nt by p testing location matrix, nt is the size of the testing locations
#         Yt            nt by 1 testing measurements
#         b             size of rep-points
#         seed          seed
# Output: a 2 by 5 matrix with the first and second rows storing energy distances and rmse for rep-points, respectively  
fsub=function(X,Y,thetat,L.decomposed,Xt,Yt,b,seed){
  n=nrow(X)
  Sid=matrix(0,5,b)
  set.seed(seed+100)
  # Select b SPs
  Sid[1,]=sp.sample(b,X,seed+100)
  # Select b Rands  
  Sid[2,]=rand.sample(b,X,seed+100)
  # Select b SP_Us
  spures=sp(b,2,dist.str = rep("uniform",2))
  Sid[3,]=subsample(X,spures$sp)
  # Select b Grids
  Sid[4,]=grid.sample(b,X)
  # Select b k-DPP
  Sid[5,]=k.dpp.sample(L.decomposed,b)
  
  ### For each kind of rep-points, calculate its energy distance and rmse
  fpred=function(id){
    # Caculate the energy distance
    id=unique(id)
    Edist=energy.dist(X,id)
    # Calculate the rmse
    rmse=tildef.hatc.rmse(X,Y,X[id,],t(as.matrix(thetat)),Xt,Yt)
    return(c(Edist,rmse))
  }
  out=apply(Sid,1,fpred)
  return(out)
}

### Build the function for reproducing results of one replicate in Figure 3(d) and Table 2
# Input: seed 
# Output: a matrix store the outputs of energy distance and rmse 
Main=function(seed){
  ### generate the training and test data
  # t1=proc.time()[[3]]
  set.seed(seed)
  n=nt=5e3
  thetat=c(1.5,0.168639,0.27,1.5)
  loc=rbind(cbind(runif(round((n+nt)*2/3),0,0.5),runif(round((n+nt)*2/3),0,0.5)),
            cbind(runif((n+nt)*1/3),runif((n+nt)*1/3)))
  C=thetat[1]*stationary.cov(loc,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  y.t=rmvnorm(1,mean=rep(0,n+nt),sigma=C,method="chol")
  Xt=loc[1:nt,]
  Yt=y.t[1:nt]
  X=loc[(nt+1):(nt+n),]
  Y=y.t[(nt+1):(nt+n)]+rnorm(n,0,sqrt(thetat[3]))
  # t2=proc.time()[[3]]
  # t2-t1=17.134
  
  
  ### Evaluate the Gaussian process prediction \hat{f}_{\hat{c}} for reference
  # t1=proc.time()[[3]]
  fullres=hatf.hatc.rmse(X,Y,t(as.matrix(thetat)),Xt,Yt)
  # t2=proc.time()[[3]]
  # t2-t1=55.734
  
  
  ### Calculate some values for k-DPPs in advance. It could save time.
  # t1=proc.time()[[3]]
  C=thetat[1]*stationary.cov(X,Covariance = "Matern",smoothness=thetat[4],theta=thetat[2])
  L.decomposed=decompose.matrix(C)
  # t2=proc.time()[[3]]
  # t2-t1=19.039
  
  
  ### Evaluate the predictive process prediction \tilde{f}_{\hat{c}} with b rep-points
  main1=function(b){
    return(fsub(X,Y,thetat,L.decomposed,Xt,Yt,b,seed))
  }
  # t1=proc.time()[[3]]
  # main1(144)
  # t2=proc.time()[[3]]
  # t2-t1=67.483
  bseq=c(36,64,100,144,196,289,400,484)
  res=sapply(bseq,main1)                     # a (2*5) by length(bseq) matrix, storing energy distances and rmse for all rep-points and all bs
  
  ### resort the outputs
  return(cbind(res,c(0,0,0,0,0,0,0,0,0,fullres)))  # a 2 by 9 matrix
  # The total computational time is about 4913.414/3600=1.36 hours
}


### Three options for reproducibility
# Option 1. Do one replicate to save time
result=Main(1)
res=result

# Option 2. Do 100 replicates as what we did 
# cl<- makeCluster(4) 
# registerDoParallel(cl) 
# result=foreach(i=1:100,
#                .combine=cbind,
#                .packages=c("geoR","fields","support","SPlit","mvtnorm",
#                             "assertthat","magrittr","purrr","rlist")) %dopar% Main(666*i+108)
# stopCluster(cl)
# res=matrix(0,10,9)
# for (i in 1:100){
#   res=res+result[,(9*i-8):(9*i)]
# }
# res=res/100


# Option 3. Load the outputs obtained by us
# result=as.matrix(read.csv(here("S2_VersusK","Setting4.csv"))[,-1])
# res=matrix(0,10,9)
# for (i in 1:100){
#   res=res+result[,(9*i-8):(9*i)]
# }
# res=res/100


### Plot Figure 3(d)
Set4Sum=data.frame(Edist=c(res[1,-9],res[3,-9],res[5,-9],res[7,-9],res[9,-9]),
                   RMSPE=c(res[2,-9],res[4,-9],res[6,-9],res[8,-9],res[10,-9]),
                   Method=as.factor(rep(c("SP","Rand","SPU","Grid","DPP"),each=8)),
                   k=rep(bseq,times=5))
# write.csv(Set4Sum,here("S2_VersusK","Set4Sum.csv"))
p=ggplot(data=Set4Sum,aes(x=k,y=RMSPE,colour=Method))+
  geom_point()+geom_line()+
  geom_hline(yintercept=res[10,9],linetype=4)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position =c(1,1),
                   #legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  #scale_linetype_manual(values=c(2,3,1,5))+
  scale_colour_manual(values=c("#0072B2","#4DAF4A","#984EA3","#E41A1C","orange"),name="Rep-points")+
  xlab("k")+ylab("RMSPE")
p









