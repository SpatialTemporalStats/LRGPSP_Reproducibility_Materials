############# Select support points (SPs)  #############
library(support)
library(SPlit)

# Purpose: Select k SPs from the given location set X.
# Inputs: k       number of SPs
#         X       n \by p location set, where n and p are the number and dimension of locations points, respectively.
#         setseed seed
# Output: Indeces of SPs w.r.t. the location set X.

sp.sample=function(k,X,setseed){
  set.seed(setseed)
  spres=sp(k,ncol(X),dist.samp=X)
  spid=subsample(X,spres$sp)
  return(spid)
}

