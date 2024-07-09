############# Evaluate Gussian process prediction \hat{f}_{\hat{c}} in (2)  #############
library(fields)

# Name: hatf.hatc.rmse
# Purpose: Calculate the Gaussian process prediction \hat{f}_{\hat{c}} in (2) with 
#          an estimated covariance function \hat{c} and full training data. Then, 
#          calculate the rooted mean squared error (RMSE) with the testing data.
# Inputs: X           n by p training location matrix, n and p are the size and dimension of the training locations, respectively
#         Y           n by 1 training measurements
#         thetaset    candidate of theta=(sigma^2,phi,tau^2,nu), number of candidate \by 4
#         Xt          nt by p testing location matrix, nt is the size of the testing locations
#         Yt          nt by 1 testing measurements
# Output: RMSE of all \hat{f}_{\hat{c}}     
hatf.hatc.rmse=function(X,Y,thetaset,Xt,Yt){
  Hatf1=function(theta){
    CH=theta[1]*stationary.cov(X,Xt,Covariance = "Matern",smoothness=theta[4],theta=theta[2])
    C=theta[1]*stationary.cov(X,Covariance = "Matern",smoothness=theta[4],theta=theta[2])
    resc=svd(C+diag(theta[3],nrow(X)))
    Yfhat=t(CH)%*%resc$u%*%diag(1/resc$d)%*%t(resc$v)%*%Y
    return(sqrt(mean((Yt-Yfhat)^2)))
  }
  return(apply(thetaset,1,Hatf1))
}



# Name: tildef.hatc.rmse
# Purpose: Calculate the predictive process prediction \tilde{f}_{\hat{c}} in (3) with 
#          an estimated covariance function \hat{c} and rep-points. Then, 
#          calculate the rooted mean squared error (RMSE) with the testing data.
# Inputs: X           n by p training location matrix, n and p are the size and dimension of the training locations, respectively
#         Y           n by 1 training measurements
#         thetaset    candidate of theta=(sigma^2,phi,tau^2,nu), number of candidate \by 4
#         Xt          nt by p testing location matrix, nt is the size of the testing locations
#         Yt          nt by 1 testing measurements
# Output: RMSE of all \hat{f}_{\hat{c}}     
tildef.hatc.rmse=function(X,Y,repX,thetaset,Xt,Yt){
  Tildef1=function(theta){
    Ch=theta[1]*stationary.cov(X,repX,Covariance = "Matern",smoothness=theta[4],theta=theta[2])
    Cth=theta[1]*stationary.cov(Xt,repX,Covariance = "Matern",smoothness=theta[4],theta=theta[2])
    Ckk=theta[1]*stationary.cov(repX,Covariance = "Matern",smoothness=theta[4],theta=theta[2])
    resck=svd(t(Ch)%*%Ch+theta[3]*Ckk)
    ythat=Cth%*%resck$u%*%diag(1/resck$d)%*%t(resck$v)%*%t(Ch)%*%Y
    return(sqrt(mean((Yt-ythat)^2)))
  }
  return(apply(thetaset,1,Tildef1))
}




# Name: energy.dist
# Purpose: Calculate the energy distance between locations and rep-points. 
# Inputs: X           n by p location matrix
#         repid       indeces of rep-points according to X
# Output: Energy distance.
energy.dist=function(X,repid){
  repX=X[repid,]        
  DXY=rdist(X,repX)
  DXX=rdist(X)
  DYY=rdist(repX)
  return(2*mean(DXY)-mean(DXX)-mean(DYY))
}




