########################################################
## Goal: Create Functions used in GAR(1) Model Selection
## 2024-09-12
#########################################################





######################
##### GAR(1): log-likelihood & BIC given Normal with Simga^{-1}=Omega=(theta0+theta1*L)^2
###### loglikelihood 
LogLike<-function(S, theta0, theta1 = 1, L, n){
  ##S: p by p sample covariance matrix  
  ## n -- sample size
  ## theta0, theta 1: GAR(1) parameters; theta1:  for the Laplacian model, theta1=1;
  ## L: p by p Laplacian matrix (normalized or not)
  p=nrow(S)
  temp=diag(theta0, p)+theta1*L
  
  temp1=sum(diag(S%*%temp%*%temp))
  temp2=determinant(temp, logarithm = TRUE)$modulus[[1]]
  
  loglike=-n/2*temp1+n*temp2-n*p/2*log(2*pi)
  
  return(loglike)
  
}


##### BIC
BIC<-function(loglike, n, k){
  ## loglike: log-likelihood 
  ## n: sample size 
  ## k: number of parameters in the model
  bic=k*log(n)-2*loglike
  return(bic)
}



########################################
#### Functions used for GF.fit #########
########################################

###### generate data according to GAR(1) model: Simga^{-1}=Omega=(theta0+theta1*L)^2
GenData.L2<-function(n,theta0,theta1, L, rep=1){
  ##n: sample size
  ## theta0>0, theta1>0
  ##L: (normalized)  Laplacian, p by p matrix 
  ## rep: number of replicates 
  ## return: list of n by p data matrices; p by p concentration matrix: Omega; p by p covariance matrix: Sigma
  
  library(mnormt)
  data.rep=NULL
  
  p=nrow(L)
  Omega=(theta0*diag(1,p)+theta1*L)%*%(theta0*diag(1,p)+theta1*L)
  Sigma=solve(Omega)
  
  for (i in 1:rep){
    data=rmnorm(n, mean=rep(0,p), varcov=Sigma)
    data.rep[[i]]=data
  }
  
  return(list(data=data.rep, Omega=Omega, Sigma=Sigma))  
}



