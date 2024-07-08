################################################################################
# This file includes steps and outputs to reproduce Figure 2 and Table 1       #
################################################################################
# Necessary packages and functions
# library(fields)
# library(ggplot2)
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
# Repeat the above procedures for multiple (e.g.,five) times and take the average. 
res=matrix(0,6,4)
for(i in 1:5){
  res=res+Mainphi(2*i)
}
res=res/5
# phi=0.168639, res table
#             2e3      4e3      6e3     8e3     \hat\gamma
# 0.5      1.417090 1.424805 1.430408 1.430501     1.43 
# 1.0      2.207545 2.211646 2.217445 2.215701     2.22 
# 1.5      2.901197 2.901686 2.906587 2.903141     2.90 
# 2.0      3.548410 3.545074 3.548163 3.542762     3.55 
# 3.0      4.781148 4.768472 4.765344 4.753931     4.75        
# 4.5      6.407396 5.563310 4.772543 4.149624     6.41
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
res1=matrix(0,6,6)
for(i in 1:5){
  res1=res1+Mainphi(2*i)
}
res1=res1/5
# phi=0.063240, res1 table
#           2e3      4e3      6e3      8e3     10e3     12e3      \hat\gamma
# 0.5    1.270816 1.328523 1.356096 1.369033 1.379598 1.385001       1.39
# 1.0    2.000287 2.075881 2.112889 2.129332 2.143759 2.149890       2.15 
# 1.5    2.631077 2.724759 2.770352 2.790623 2.808241 2.815972       2.82
# 2.0    3.211935 3.324470 3.378250 3.402420 3.422997 3.432680       3.43    
# 3.0    4.299955 4.452119 4.521856 4.553649 4.580061 4.593407       4.60
# 4.5    5.853063 6.068613 6.157206 6.181066 6.169712 6.122045       6.12 

### Plot Figure 2(a)
Data=data.frame(Gam=c(1.43,2.22,2.90,3.55,4.75,6.41,1.39,2.15,2.82,3.42,4.60,6.12),
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
Mainphi=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \nu=1.5 and let \phi=0.168639*c(1/5,2/5,3/5,4/5,1,6/5)
  Thetaset=cbind(rep(1.5,6),0.168639*c(1/5,2/5,3/5,4/5,1,6/5),rep(1.5,6))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
# Repeat the above procedures for multiple (e.g.,five) times and take the average. 
res=matrix(0,6,4)
for(i in 1:5){
  res=res+Mainphi(2*i)
}
res=res/5
# \nu=1.5, res table 
#               2e3       4e3      6e3      8e3      \hat\gamma
# 0.0337278   2.274080 2.469064 2.564640 2.615883      2.62
# 0.0674556   2.657949 2.743036 2.784688 2.802604      2.80
# 0.1011834   2.793378 2.832965 2.854395 2.860421      2.86
# 0.1349112   2.860974 2.876380 2.887494 2.887580      2.89
# 0.1686390   2.901197 2.901686 2.906587 2.903141      2.90
# 0.2023668   2.927828 2.918198 2.918954 2.913171      2.91

### Evaluate \gamma for \nu=1.0
Mainphi=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \nu=1.0 and let \phi=0.168639*c(1/5,2/5,3/5,4/5,1,6/5)
  Thetaset=cbind(rep(1.5,6),0.168639*c(1/5,2/5,3/5,4/5,1,6/5),rep(1,6))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
res1=matrix(0,6,6)
for(i in 1:5){
  res1=res1+Mainphi(2*i)
}
res1=res1/5
# \nu=1.0, res1 table
#                6e3      8e3     10e3    12e3      \hat\gamma
# 0.0337278   1.721353 1.875860 1.951904 1.992558      2.03          
# 0.0674556   2.021162 2.090079 2.124025 2.138640      2.15                 
# 0.1011834   2.125724 2.159501 2.177837 2.183276      2.19
# 0.1349112   2.177254 2.192590 2.203065 2.203979      2.20
# 0.1686390   2.207545 2.211646 2.217445 2.215701      2.22
# 0.2023668   2.227370 2.223939 2.226654 2.223172      2.23

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
Mainphi=function(setseed){
  # Set consistent parameters 
  Thetaset=rbind(c(1.5,0.168639,1.5),c(1,0.147320,1.5),c(2,0.185611,1.5))
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
# Repeat the above procedures for multiple (e.g.,five) times and take the average. 
res=matrix(0,3,4)
for(i in 1:5){
  res=res+Mainphi(2*i)
}
res=res/5
# consistent parameters, res table
# 2.892871 2.895336 2.902760 2.905831
# 2.869637 2.880762 2.891783 2.896902
# 2.907493 2.904432 2.909582 2.911366



