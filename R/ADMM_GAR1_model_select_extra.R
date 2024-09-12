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

