############# Calculate the smoothness parameter \gamma  #############
library(fields)

# Name: gammafuc
# Purpose: Evaluate the value of \gamma for the Matern covariance function in (7) basing on the method descripted on Page 14
# Inputs: n           size of locations
#         thetaset    candidate of theta=(sigma^2,phi,nu), number of candidate \by 3
#         setseed     seed
# Output: \gamma values for all candidate theta

gammafuc=function(n,thetaset,setseed){
  set.seed(setseed)
  loc=cbind(runif(n),runif(n))
  Ga=function(theta){
    C=theta[1]*stationary.cov(loc,Covariance = "Matern",smoothness=theta[3],theta=theta[2])
    resc=svd(C)$d
    out=-lm(log(resc)~log(1:n))$coefficients[2]
  }
  return(apply(thetaset,1,Ga))
}