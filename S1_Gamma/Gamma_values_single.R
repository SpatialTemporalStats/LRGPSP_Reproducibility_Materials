########################################################################################
# This file reproduces Figure 2 and Table 1 using results from a single replicate      #
########################################################################################    
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
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
# t1=proc.time()[[3]]
# res=Mainphi(1)
# t2=proc.time()[[3]]
# t2-t1=1814.81
# Note that you may consider to repeat the above procedures for multiple times and then 
# take the average as the final result. Please refer to the script "Gamma_values.R" for
# more details.


### Evaluate \gamma for \phi=0.063240
Mainphi=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \phi=0.063240 and let \nu=0.5, 1, 1.5, 2, 3, 4.5
  Thetaset=cbind(rep(1.5,6),rep(0.063240,6),c(0.5,1,1.5,2,3,4.5))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3,10e3,12e3),G1))
  # t1=proc.time()[[3]]
  # gammafuc(12e3,Thetaset,1)
  # t2=proc.time()[[3]]
  # t2-t1=3920.124
}
res1=Mainphi(1)


### Plot Figure 2(a)
Data=data.frame(Gam=c(res[1:5,4],res[6,1],res1[,6]),
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
res=Mainphi(1)


### Evaluate \gamma for \nu=1.0
Mainphi=function(setseed){
  # Set parameters theta=(sigma^2,phi,nu): fix \sigma^2=1.5, \nu=1.0 and let \phi=0.168639*c(1/5,2/5,3/5,4/5,1,6/5)
  Thetaset=cbind(rep(1.5,6),0.168639*c(1/5,2/5,3/5,4/5,1,6/5),rep(1,6))  
  # Evaluate \gamma under the given parameter settings with different sizes of locations
  G1=function(n){return(gammafuc(n,Thetaset,setseed))}
  return(sapply(c(2e3,4e3,6e3,8e3),G1))
}
res1=Mainphi(1)


### Plot Figure 2(b)
Data=data.frame(Gam=c(res[,4],res1[,4]),
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
res=Mainphi(1)

