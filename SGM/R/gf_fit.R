#############
##### goodness of fit by parametric bootsrap
##### Added 2024-09-27: Use Step 1 Estimates 

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

bootstrap.like <- function(L, theta0, theta1, n, lambda.v, rho.v=lambda.v, rep.boot=100, num.thread = 1, seed = 1){
  
  ## Generate bootstrap samples
  print("Samples generating...")
  set.seed(1)
  temp = GenData.L2(n, theta0, theta1, L, rep.boot)
  data.boot = temp$data # Extract datasets 
  print(paste("Fitting Bootstrap Samples..."))

  p = nrow(L)
  
  ## Refit GAR(1) with step 1 to bootstrap samples
  temp = foreach(i = 1:rep.boot, .combine = 'c', .packages = c("mnormt")) %dopar% {
    
    ## Sample Covariance Matrix
    S.c = var(data.boot[[i]])*(n-1)/n 
    theta0.i.ini = 1/sqrt(max(eigen(S.c, symmetric = T, only.values = T)$values))
    
    ## Fit GAR(1) to each replicate
    fit.rep.i = ADMM_L2(S.c, theta0.i.ini, v=rep(0,ncol(S.c)), rho.v, lambda.v, model="LN", Z_ini=matrix(0, p, p), W_ini=matrix(0, p, p), 
                        eps_thre=1e-6, eps_abs=1e-5, eps_rel=1e-3, max_iter=100000, verbose=FALSE)
    theta1.i = fit.rep.i$theta1
    L.i = fit.rep.i$L
    LogLike(S.c, theta0.i.ini, theta1.i, L.i, n)
  }
  return(temp)
}




#' Goodness of Fit Test
#' 
#' @description
#' This function provides a goodness of fit test to see whether a GAR(1) model with the normalized laplacian is applicable.  It is valid for when the signal dimension is a most the number of observations available.
#' 
#' @param S Estimate for covariance matrix, such as the MLE
#' @param nobs The number of observations used to calculate `S`
#' @param lambda.v Tuning parameter for GAR(1).  Positive number.
#' @param rho.v ADMM parameter.  Positive number.
#' @param eps_thre Small positive number.
#' @param eps_abs ADMM convergence criterion.
#' @param eps_rel ADMM convergence criterion. 
#' @param max_iter Number of iterations to run initial fit on observed data
#' @param num.threads Number of threads to use for computing.
#' @param rep.boot Number of bootstrap samples to generate for test.
#' 
#' @returns p-value for the goodness of fit test 
#' 
#' @export
GAR1_gf = function(S, nobs, lambda.v, rho.v=lambda.v, eps_thre = 1e-6, eps_abs = 1e-5, eps_rel = 1e-3, max_iter = 10000, num.thread = 1, rep.boot = 100, seed = 1){
  
  ## Reserve cores
  registerDoParallel(num.thread)
  
  ## Fit step 1 to observed data
  init.step = fit_step_0a(S, nobs)
  init.fit = fit_step_1a(init.step, lambda.v, rho.v, "LN", eps_thre, eps_abs, eps_rel, max_iter, F)
  
  print("Initial fit completed")

  ## Extract the fit results
  theta0.est = init.step$theta0
  theta1.est = init.fit[[1]]$theta1
  L.est = init.fit[[1]]$L
  
  ## Run bootstrap resamples and store log-likleihoods for each bootstrap sample
  log.like.boot = bootstrap.like(L.est, theta0.est, theta1.est, nobs, lambda.v, rho.v, rep.boot, num.thread, seed)

  
  ## Calculate observed log-likelihood
  log.like.obs = LogLike(S, theta0.est, theta1.est, L.est, nobs)
  
  ## Calculate vector of whether obs is greater than boot
  pvals = log.like.obs > log.like.boot

  
  ## Return mean (p-value)
  print("p-value:")
  return(mean(pvals))
}