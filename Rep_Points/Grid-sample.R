############# Select grid points (Grids) #############
library(SPlit)

# Purpose: Select k Grids from the given location set X.
# Inputs: k       number of Grids
#         X       n \by p location set, where n and p are the number and dimension of locations points, respectively.
# Output: Indeces of Grids w.r.t. the location set X.

grid.sample=function(k,X){
  kk=round(sqrt(k))
  gr=seq(0,1,length=kk+2)[-c(1,kk+2)]
  Gr=cbind(rep(gr,each=kk),rep(gr,times=kk))
  gridid=subsample(X,Gr)
  return(gridid)
}