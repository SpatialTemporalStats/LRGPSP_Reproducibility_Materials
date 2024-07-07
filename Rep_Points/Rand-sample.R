############# Select random subsamples (Rands) #############
# Purpose: Select k Rands from the given location set X.
# Inputs: k       number of Rands.
#         X       n \by p location set, where n and p are the number and dimension of locations points, respectively.
#         setseed seed
# Output: Indeces of Rands w.r.t. the location set X.

rand.sample=function(k,X,setseed){
  set.seed(setseed)
  rid=sample(nrow(X),k,replace=FALSE)
  return(rid)
}