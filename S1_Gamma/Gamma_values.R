################################################################################
# This file includes steps and outputs to reproduce Figure 2 and Table 1       #
################################################################################
# Necessary packages and functions
# library(fields)
# library(ggplot2)
# library(doParallel)
# library(foreach)
# source(here("S1_Gamma","helper-functions.R"))



################# Reproduce Figure 2(a) #################
### Evaluate \gamma for \phi=0.168639
Mainphi=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \phi=0.168639 and let \nu=0.5, 1, 1.5, 2, 3, 4.5
  Thetaset=cbind(rep(1.5,6),rep(0.168639,6),c(0.5,1,1.5,2,3,4.5))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  # t1=proc.time()[[3]]
  # gammafuc(12e3,Thetaset,1)
  # t2=proc.time()[[3]]
  # t2-t1=3920.124
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
# Repeat the above procedures for ten times and take the average. 
cl=makeCluster(4)
registerDoParallel(cl)
result=foreach(i=1:10,
               .combine = cbind,
               .packages = c("fields")) %dopar% Mainphi(2*i)
stopCluster(cl)
# write.csv(result,here("S1_Gamma","Res_Phi1.csv"))
result=as.matrix(read.csv(here("S1_Gamma","Res_Phi1.csv"))[,-1])
res=matrix(0,6,4)
for(i in 1:10){
  res=res+result[,(4*(i-1)+1):(4*i)]
}
res=res/10
# phi=0.168639, res table
#             2e3      4e3      6e3     8e3     \hat\gamma
# 0.5      1.411752 1.424425 1.430085 1.430298     1.43 
# 1.0      2.198637 2.211150 2.217124 2.215376     2.22 
# 1.5      2.889874 2.901090 2.906350 2.902640     2.90 
# 2.0      3.535389 3.544425 3.547988 3.542143     3.55 
# 3.0      4.765963 4.767902 4.765141 4.753011     4.75        
# 4.5      6.401526 5.564007 4.775728 4.158396     6.41
# Note. For most cases, the value of \hat\gamma would converege as n grows. However,
#       when \nu=4.5, the trend seems wired. Take the fitting results (0.9775+1.2316*\nu) 
#       into consideration, we choose \hat\gamma=6.41 when \nu=4.5.

### Evaluate \gamma for \phi=0.063240
Mainphi=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \phi=0.063240 and let \nu=0.5, 1, 1.5, 2, 3, 4.5
  Thetaset=cbind(rep(1.5,6),rep(0.063240,6),c(0.5,1,1.5,2,3,4.5))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3,10e3,12e3),G1))
}
# Repeat the above procedures for ten times and take the average. 
cl=makeCluster(4)
registerDoParallel(cl)
result=foreach(i=1:10,
               .combine = cbind,
               .packages = c("fields")) %dopar% Mainphi(2*i)
stopCluster(cl)
# write.csv(result,here("S1_Gamma","Res_Phi2.csv"))
result=as.matrix(read.csv(here("S1_Gamma","Res_Phi2.csv"))[,-1])
res=matrix(0,6,6)
for(i in 1:10){
  res=res+result[,(6*(i-1)+1):(6*i)]
}
res=res/10
# phi=0.063240, res1 table
#           2e3      4e3      6e3      8e3     10e3     12e3      \hat\gamma
# 0.5    1.265491 1.328137 1.355779 1.368845 1.378955 1.384615       1.38
# 1.0    1.991430 2.075371 2.112578 2.129035 2.142722 2.149258       2.15 
# 1.5    2.619812 2.724140 2.770128 2.790161 2.807289 2.815188       2.82
# 2.0    3.198972 3.323787 3.378091 3.401850 3.422296 3.431732       3.43    
# 3.0    4.284815 4.451496 4.521673 4.553100 4.579794 4.592055       4.60
# 4.5    5.836610 6.068308 6.156176 6.180876 6.171721 6.123283       6.12 

### Plot Figure 2(a)
Data=data.frame(Gam=c(1.43,2.22,2.90,3.55,4.75,6.41,1.38,2.15,2.82,3.42,4.60,6.12),
                nu=rep(c(0.5,1.0,1.5,2.0,3.0,4.5),times=2),
                phi=as.factor(rep(c("\u03C8=0.169","\u03C8=0.063"),each=6)))
lmres=lm(Gam~nu,Data[1:6,])
lmres1=lm(Gam~nu,Data[-(1:6),])
Data$Gamhat=c(lmres$fitted.values,lmres1$fitted.values)
p=ggplot(data=Data,aes(x=nu,y=Gam,colour=phi))+
  geom_point()+
  geom_line(aes(x=nu,y=Gamhat))+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=14),
                   axis.text.y = element_text(angle=90,size=14),
                   axis.title=element_text(size=14),
                   legend.justification=c(0,1),
                   legend.position =c(0,1),
                   legend.text = element_text(size=14),
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_colour_manual(values=c("#4DAF4A","#E41A1C"))+
  xlab(expression(nu))+ylab(expression(gamma))
p



################# Reproduce Figure 2(b) #################
### Evaluate \gamma for \nu=1.5
Mainnu=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \nu=1.5 and let \phi=0.168639*c(1/5,2/5,3/5,4/5,1,6/5)
  Thetaset=cbind(rep(1.5,6),0.168639*c(1/5,2/5,3/5,4/5,1,6/5),rep(1.5,6))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
# Repeat the above procedures for ten times and take the average. 
cl=makeCluster(4)
registerDoParallel(cl)
result=foreach(i=1:10,
               .combine = cbind,
               .packages = c("fields")) %dopar% Mainnu(2*i)
stopCluster(cl)
# write.csv(result,here("S1_Gamma","Res_nu1.csv"))
result=as.matrix(read.csv(here("S1_Gamma","Res_nu1.csv"))[,-1])
res=matrix(0,6,4)
for(i in 1:10){
  res=res+result[,(4*(i-1)+1):(4*i)]
}
res=res/10
# \nu=1.5, res table 
#               2e3       4e3      6e3      8e3      \hat\gamma
# 0.0337278   2.262864 2.468430 2.564422 2.615463      2.62
# 0.0674556   2.646679 2.742418 2.784463 2.802138      2.80
# 0.1011834   2.782078 2.832357 2.854164 2.859936      2.86
# 0.1349112   2.849660 2.875779 2.887259 2.887085      2.89
# 0.1686390   2.889874 2.901090 2.906350 2.902640      2.90
# 0.2023668   2.916499 2.917606 2.918715 2.912666      2.91

### Evaluate \gamma for \nu=1.0
Mainnu=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \nu=1.0 and let \phi=0.168639*c(1/5,2/5,3/5,4/5,1,6/5)
  Thetaset=cbind(rep(1.5,6),0.168639*c(1/5,2/5,3/5,4/5,1,6/5),rep(1,6))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
cl=makeCluster(4)
registerDoParallel(cl)
result=foreach(i=1:10,
               .combine = cbind,
               .packages = c("fields")) %dopar% Mainnu(2*i)
stopCluster(cl)
# write.csv(result,here("S1_Gamma","Res_nu2.csv"))
result=as.matrix(read.csv(here("S1_Gamma","Res_nu2.csv"))[,-1])
res=matrix(0,6,4)
for(i in 1:10){
  res=res+result[,(4*(i-1)+1):(4*i)]
}
res=res/10
# \nu=1.0, res1 table
#                6e3      8e3     10e3    12e3      \hat\gamma
# 0.0337278   1.712502 1.875339 1.951593 1.992286      2.00          
# 0.0674556   2.012301 2.089570 2.123714 2.138340      2.14                 
# 0.1011834   2.116839 2.158998 2.177521 2.182963      2.18
# 0.1349112   2.168355 2.192091 2.202746 2.203659      2.20
# 0.1686390   2.198637 2.211150 2.217124 2.215376      2.22
# 0.2023668   2.218456 2.223445 2.226332 2.222844      2.223

### Plot Figure 2(b)
Data=data.frame(Gam=c(2.62,2.80,2.86,2.89,2.90,2.91,
                      2.00,2.14,2.18,2.20,2.216,2.223),
                phi=c(0.0337278,0.0674556,0.1011834,0.1349112,0.1686390,0.2023668,
                      0.0337278,0.0674556,0.1011834,0.1349112,0.1686390,0.2023668),
                nu=as.factor(rep(c("\u03BD=1.5","\u03BD=1.0"),each=6)))
p=ggplot(data=Data,aes(x=phi,y=Gam,colour=nu))+
  geom_point()+geom_line()+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=14),
                   axis.text.y = element_text(angle=90,size=14),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position =c(1,0.8),
                   legend.text = element_text(size=14),
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  scale_colour_manual(values=c("#4DAF4A","#E41A1C"))+
  xlab(expression(psi))+ylab(expression(gamma))
p



################# Reproduce Table 1 #################
### Evaluate \gamma for consistent parameter settings
MainConsistent=function(setseed){
  # Set consistent parameters 
  Thetaset=rbind(c(1.5,0.168639,1.5),c(1,0.147320,1.5),c(2,0.185611,1.5))
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
# Repeat the above procedures for ten times and take the average.
cl=makeCluster(5)
registerDoParallel(cl)
result=foreach(i=1:10,
               .combine = cbind,
               .packages = c("fields")) %dopar% MainConsistent(2*i)
stopCluster(cl)
# write.csv(result,here("S1_Gamma","Res_Consistent.csv"))
result=as.matrix(read.csv(here("S1_Gamma","Res_Consistent.csv"))[,-1])
res=matrix(0,3,4)
for(i in 1:10){
  res=res+result[,(4*(i-1)+1):(4*i)]
}
res=res/10
# consistent parameters, res table
# 2.889874 2.901090 2.906350 2.902640
# 2.866635 2.886514 2.895376 2.893709
# 2.904499 2.910187 2.913171 2.908175



