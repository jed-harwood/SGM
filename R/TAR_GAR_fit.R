###############################
## TAR-GAR Fitting Functions
## Added: 2025-02-22
## TAR GAR Fitting Procedure and Necessary Functions;
## Using Version A of Fitting Procedure
###############################



###############################
### PREREQUISITE FUNCTIONS  ###
###############################

#' Calculate sample auto-covariance matrix at lag=`order`
#' 
#' @description
#' Calculates the sample autocovariance matrix at a specified lag/ordering.
#' 
#' @param data n by p data matrix
#' @param order nonnegative integer
#' 
#' @returns Returns a p-by-p symmetric matrix.
#' 
AutoCov<-function(data, order=0){
  ## data: n by p data matrix
  ## order: nonnegative integer 
  ## return: Gamma_order: p by p symmetric matrix 
  
  n= nrow(data) ## sample size 
  p = ncol(data)
  
  if(order<0 || order>n-1){
    stop("order must be an integer between {0,..., n-1}")
  }
  
  data.m=matrix(apply(data, 2, mean), n, p, byrow=TRUE) ## mean matrix
  data.c=data-data.m  ## centered data 
  
  
  data.c.pre=data.c[1: (n-order),] ##(n-order) by p 
  data.c.lag=data.c[(1+order):n, ]
  
  result=(t(data.c.pre)%*%data.c.lag)/(n-order) ## Gamma_order: when order=0, this should be the same as var(data)*(n-1)/n
  result=(result+t(result))/2 ## symmetrizing 
  
  return(result) 
  
}


#' Initial Estimate of Filter Matrix
#' 
#' @description
#' Calculates an initial estimate of the graph filter matrix R1. 
#' 
#' @param data An n by p data matrix
#' @param eps.thre Threshold for pseudoinverse of Gamma0
#' 
#' @returns A p by p symmetrix matrix

step.00 <-function(data,  eps.thre=1e-6){
  ## data: n by p data matrix
  ##  eps.thre: threshold for pseudo inverse of Gamma0
  ##  return: R1 : p by p symmetric
  
  p=ncol(data)
  ##
  Gamma0=AutoCov(data, order=0)
  Gamma0.eigen=eigen(Gamma0, symmetric=TRUE)
  
  Gamma0.lam=Gamma0.eigen$values
  Gamma0.lam=pmax(Gamma0.lam,0) ## remove small negative values due to rounding error  
  
  ## Gamma0^{-0.5}: use pseudo inverse when Gamma0 is singular 
  Gamma0.lam.neg0.5=numeric(p)
  Gamma0.lam.neg0.5[Gamma0.lam>eps.thre]=(Gamma0.lam[Gamma0.lam>eps.thre])^{-0.5}
  Gamma0.neg0.5=(Gamma0.eigen$vectors)%*%diag(Gamma0.lam.neg0.5)%*%t(Gamma0.eigen$vectors)
  
  ## R1 estimate 
  Gamma1=AutoCov(data, order=1)
  R1=Gamma0.neg0.5%*%Gamma1%*%Gamma0.neg0.5
  R1=(R1+t(R1))/2 ## remove asymmetry due to rounding error 
  
  ##
  return(R1)
}


#' Sample covariance of the innovation process  
#' 
#' @description
#' Calculates the sample covariance of the innovation process, and the corresponding theta0 estimate, given a filter matrix R1
#' 
#' @param data n by p matrix
#' @param R1 p by p matrix
#' 
#' @returns A list containing:
#' * `S`: A p by p sample covariance matrix of the innovation process
#' * `theta0`: A positive number
#' 
step.0 <-function(data, R1){
  ## data: n by p data matrix
  ## R1: p by p TAR filter 
  ## return: S: p by p sample covariance of the innovation process U.t= Y.t-R1%*%Y.{t-1}	
  ## theta0: nonnegative scalar
  
  data.pre=data[1:(n-1),] ## (n-1) by p: from t=1,...,n-1
  data.lag=data[2:n,] ##(n-1) by p: from t=2,..,n
  residual=data.lag-data.pre%*%R1 ## (n-1) by p: innovation process  
  
  S=var(residual)*(n-2)/(n-1)
  theta0=(max(eigen(S, symmetric=TRUE)$values))^{-0.5}
  
  return(list(S=S, theta0=theta0))
}



#' Step 2 of TAR-GAR Fitting Procedure
#' 
#' @description
#' Estimate the filter parameters (eta0, eta1) given (theta0, L), through a quadratic program with linear inequality constraints
#' 
#' @param data An n by p data matrix
#' @param L A p by p (normalized) Laplacian Matrix
#' @param theta0 A positive number
#' @param eps A (small) positive number
#' @importFrom quadprog solve.QP
step.2 <- function(data, L, theta0, eps=1e-6){
  ## data: n by p data matrix
  ## L: p by p Laplacian; p.s.d.
  ## theta0: positive scalar 
  ## eps: small positive values to bound |eta0+eta1*lambda_j|<= 1-eps
  ## center data:
  n=nrow(data)
  p=ncol(data)
  
  data.m=matrix(apply(data, 2, mean), n, p, byrow=TRUE) ## mean matrix
  data.c=data-data.m  ## centered data 
  
  ## eigen-decomposition of the Laplacian 
  eigen.L=eigen(L, symmetric=TRUE)
  Q=t(eigen.L$vectors)  ##  L= t(Q)%*%diag(lambda)%*%Q
  lambda=eigen.L$values
  
  ## transform data by L's eigenvectors 
  data.tilde= data.c%*%t(Q)
  
  ## Gamma0 and Gamma1 of the transformed data 
  Gamma0 = t(data.tilde[2:n,])%*%data.tilde[2:n,]/(n-1)
  Gamma1 = t(data.tilde[2:n,])%*%data.tilde[1:(n-1),]/(n-1)
  Gamma1 = (Gamma1+t(Gamma1))/2
  
  ## QP components
  ## dvec 
  d.1= 2*sum((theta0+lambda)^2*diag(Gamma1))
  d.2= 2*sum(lambda*(theta0+lambda)^2*diag(Gamma1))
  d=c(d.1,d.2)
  
  #Dmat
  d.11= 2*sum((theta0+lambda)^2*diag(Gamma0))
  d.22= 2*sum(lambda^2*(theta0+lambda)^2*diag(Gamma0))
  d.12= 2*sum(lambda*(theta0+lambda)^2*diag(Gamma0))
  D=matrix(c(d.11, d.12, d.12, d.22), 2,2)
  
  #Amat: 2 by (2p) 
  A.left=rbind(rep(1,p), lambda)
  A.right=-A.left
  A=cbind(A.left, A.right)
  
  #bvec: (2p) by 1 
  b=rep(-1,2*p)+eps
  
  # minimize -d^T eta + 1/2*eta^T D eta subject: A^T eta >= b
  result=try(solve.QP(Dmat=D, dvec=d, Amat=A, bvec=b))
  if(is(result, "try-error")){
    return(NULL)
    
  }else{
    # return eta and R1
    eta = result$solution
    R1=diag(eta[1],p)+eta[2]*L
  }
  
  #
  return(list(result=result, eta=eta, R1=R1))
}


#' Step 3 of TAR-GAR Fitting Procedure
#' 
#' @description
#' Update R1  given the estimated (eta0, eta1) and L parameters from steps 0-to-2.
#' 
#' @param eta0 A number
#' @param eta1 A number
#' @param L A p-by-p matrix
step.3 = function(eta0, eta1, L){
  p = nrow(L)
  R1 = diag(eta0, p)+eta1*L
  return(R1)
}

##########################################################################
## Extra Step: Re-estimate Final-Pass L using net.thre values, 
##             with option to jointly estimate theta0 and L, or separately.
##########################################################################

#' Thresholding Step
#' 
#' @description
#' Use net.thre values to invoke sparsity in final pass estimate of L.  Takes in resList from final pass, and produces list of zero-patterns
#' @param resList Output from passes of TAR-GAR fitting procedure
#' @param lambda.v A vector of positive numbers
#' @param net.thre A vector of (small) positive numbers
step.thre = function(resList, lambda.v, net.thre){
  A.net.e = vector(mode = "list", length = length(lambda.v))
  
  for (j in 1:length(lambda.v)){
    
    A.net.e[[j]] = vector(mode = "list", length = length(net.thre))
    
    L.lambda = resList[[j]]$L * resList[[j]]$theta1 ### L = Lhat * theta1hat
    
    for (k in 1:length(net.thre)){
      temp = abs(L.lambda) > net.thre[k]
      diag(temp) = 0
      A.net.e[[j]][[k]] = temp
    }
    
    
  }
  return(A.net.e)
  
}

#' Refitting Step of TAR-GAR Procedure
#'
#' @description
#' Given 0 pattern from net.thre, give an option to re-estimate L using GAR-strategy.  Takes in 0 pattern and results from the TAR-GAR passes, returns 3-step estimation procedure
#' 
#' @param data An n-by-p matrix
#' @param resList A list
#' @param A.net.e A p-by-p matrix
#' @param lambda.v A vector of positive numbers
#' @param net.thre A vector of (small) positive numbers
#' @param refit A positive integer
#' @param model A character
#' @param eps_thre A small positive number
#' @param eps_abs A small positive number
#' @param eps_rel A small positive number
#' @param max_iter A large integer
#' @param verbose A boolean
#' @param timevec A vector
step.lap.est = function(data, resList, A.net.e, lambda.v, net.thre, refit = 2, model = "LN", eps_thre = 1e-6, eps_abs = 1e-5, eps_rel = 1e-3, max_iter = 50000, verbose = FALSE, time.vec){
  results = vector(mode = "list", length = length(lambda.v))
  p = ncol(data) 
  n = nrow(data) - 1
  ini.mat = matrix(0, p, p)
  ini.vec = rep(0,p)
  
  for (j in 1:length(lambda.v)){
    results[[j]] = vector(mode = "list", length = length(net.thre))
    resList.j = resList[[j]] # Extract final pass results for jth lambda
    S.c = resList.j$S # Final pass S
    theta0.c = resList.j$theta0 # Final pass theta0
    
    for(k in 1:length(net.thre)){
      
      ## time for refitting (start)
      time.start = proc.time()[3]
      
      A.c = A.net.e[[j]][[k]] # Current zero-pattern
      
      ## Step 2 of GAR (refit off-diagonal entries)
      refit.2 = try(ADMM_L2_Zero(SS = S.c, theta0 = theta0.c, v = ini.vec, rho = sqrt(log(p)/n), AA = A.c, model = model, Z_ini = ini.mat, W_ini = ini.mat, eps_thre = eps_thre, 
                                 eps_abs = eps_abs, eps_rel = eps_rel, max_iter = max_iter, verbose = verbose))
      conv.2 = T
      if(inherits(refit.2, "try-error")){
        refit.2 = NULL
        conv.2=FALSE
      }
      
      L.est = refit.2$L
      v0.res = NULL
      conv.v0 = FALSE
      refit.3 = NULL
      v0.est = ini.vec
      
      
      ## Step 3.1 of GAR (get v0 estimate from previous L estimator)
      if ((refit > 1) && (model != "L") && conv.2 ==T){
        v0.res = try(ADMM.Deg.L(L.est, rho = 0.1, epsilon = sqrt(1/(2*ncol(L.est))), eps.abs = 1e-05, 
                                eps.rel = 0.001, max.iter = 50000, verbose = verbose))
        v0.est = v0.res$v
        conv.v0 = v0.res$conv
        
        if(inherits(v0.res, "try-error")){
          v0.res = NULL
          v0.est = ini.vec
          conv.v0=FALSE
        }
        
        
      }
      
      ## Step 3.2 of GAR (Put L in space of (normalized) graph Laplacians given v)
      if (conv.v0 == TRUE && (model != "L") ){
        
        refit.3 = try(ADMM_Lap_Zero(SS = S.c, V0 = v0.est, rho = sqrt(log(p)/n), AA = A.c, model = model, ZZ_ini = ini.mat, WW_ini = ini.mat, phi_ini = 1e-6, eps_thre = eps_thre, 
                                    eps_abs = eps_abs, eps_rel = eps_rel, max_iter = max_iter, Z_max_iter = max_iter, Z_conv_abs = 1e-5, Z_conv_rel = 1e-3, verbose = verbose))
        conv.3 = TRUE
        if(inherits(refit.3, "try-error")){
          refit.3 = NULL
          conv.3=FALSE
        }
        
      } else{
        refit.3 = NULL
        conv.3 = FALSE
      }
      
      ## Step 2 of TAR-GAR: Refit eta0, eta1.tilde Matrix
      if (conv.3){
        L.est = refit.3$L * refit.3$theta1 
        theta0.est = refit.3$theta0
      } else{
        warning("Step 3 of Refitting Failed.  Using Step 2 Estimate for L instead.")
        L.est = refit.2$L * refit.2$theta1
        theta0.est = theta0.c
      }
      eta.refit = step.2(data = data, L = L.est, theta0 = theta0.est, eps = eps_thre)
      
      ## Step 3 of TAR-GAR: Refit R1
      eta.0.refit = eta.refit$eta[1]
      eta.1.refit = eta.refit$eta[2]
      R1.refit = step.3(eta.0.refit, eta.1.refit, L.est) # L.est has already absorbed theta1
      
      ## time for refitting (end)
      time.end = proc.time()[3]
      
      ## time for refitting
      time.refit = time.end - time.start + time.vec[j]
      
      results[[j]][[k]] = list("S" = S.c, "theta0.ini" = theta0.c, "result.0.post" = refit.2, "conv.2" = conv.2, "theta0.0S" = refit.3$theta0, "result.0S" = refit.3, "v0.est" = v0.est, "conv.3" = conv.3, "A.net" = A.c, "R1.0S" = R1.refit, "eta0.0S" = eta.0.refit, "eta1.0S" = eta.1.refit, "time" = time.refit)
      
    }
  }
  
  return(results)
}



#' TAR-GAR Fitting Procedure
#' 
#' @include RcppExports.R
#' 
#' @description
#' Using the TAR-GAR fitting procedure, infers the (normalized) graph Laplacian matrix of a latent graph from temporally dependent data.  Fits to a sequence of tuning parameters (`lambda.v`, `net.thre`).  
#' 
#' @param data An n-by-p matrix
#' @param lambda.v A vector of positive numbers.  These values are tuning parameters for the TAR-GAR procedure.
#' @param net.thre A vector of positive numbers.  These values are tuning parameters for the TAR-GAR procedure. 
#' @param rho.v A vector of positive numbers.  The step size for the ADMM algorithm.  Set to `lambda.v` by default.
#' @param num.pass A positive integer.  Indicates how many passes to run before refitting. Three passes is recommended.
#' @param model A character.  Either `"LN"`, `"LN.noselfloop"`, or `"L"`.  Corresponds to GAR structure of the TAR-GAR innovation term.
#' @param eps.thre A small positive number.
#' @param eps_abs A small positive number.  Refers to ADMM stopping criteria. By default, set to 1e-5.
#' @param eps_rel A small positive number.  Refers to ADMM stopping criteria.  By default, set to 1e-3.
#' @param max_iter A large positive integer.  Refers to ADMM stopping criteria.  By default, set to 50000.
#' @param refit An integer.  Indicated number of steps in the GAR fitting procedure to use in the refitting step of the TAR-GAR procedure.  By default, set to 2 (indicates to use full GAR fitting procedure).
#' @param verbose A boolean.  Setting to TRUE prints out trace of the fitting procedure.
#' 
#' @returns A list containing:
#' * `result.pass`: A list corresponding to the fitting results, pre-refitting.  Length is equivalent to `length(lambda.v)`.
#' * `ini`: A list corresponding to the initial fitting results and first pass.  Length is equivalent to `length(lambda.v)`.
#' * `refit`: A list corresponding to the final estimates using the TAR-GAR procedure.  Contains a collection of lists, with inner index corresponding to `lambda index`, and outer index corresponding to `net.thre` index.
#' * `A.net`: A list of matrices containing the zero-pattern of the estimated network.   Contains a collection of matrices, with inner index corresponding to `lambda index`, and outer index corresponding to `net.thre` index.
#' * `time.total`: Total time to run the TAR-GAR procedure on the specified sequence of tuning parameters.   
#' 
#' @example man-roxygen/TARGAR_fit_example.R
#' @export
TARGAR_fit = function(data, lambda.v, net.thre, rho.v=lambda.v, num.pass = 2,  model = "LN", eps.thre = 1e-6, eps_abs = 1e-5, eps_rel = 1e-3, max_iter = 50000, refit=2, verbose = F){
  
  ## Start Time
  time.tot.beg = proc.time()[3]
  
  ## Storage objects
  resList = vector(length = length(lambda.v), mode = "list")
  iniRes = resList
  time.vec = vector(length = length(lambda.v))
  
  ## Sample size
  n = nrow(data)
  
  ## Step 00
  R1.fit = step.00(data, eps.thre)
  R1.ini = R1.fit
  print("Step 00 Complete")
  
  for (j in 1:length(lambda.v)){
    print(paste("Starting fitting for lambda =", lambda.v[j]))
    time.start = proc.time()[3]
    
    for (i in 1:num.pass){
      
      
      ## Step 0
      fit.S.theta0 = step.0(data, R1.fit)
      S.i = fit.S.theta0$S
      theta0.i = fit.S.theta0$theta0
      print("step 0 complete")
      
      ## Step 1
      p = nrow(S.i)
      
      ## v0 for model; not necessary, since ADMM_L2 automatically will change v0 depending on model
      
      if(model == "L"){
        v0 = rep(1, p)
      } else{
        v0 = rep(0,p)
      }
      fit.L.theta1 = try(ADMM_L2(s = S.i, theta0 = theta0.i, v = v0, rho = rho.v[j], lambda =lambda.v[j],
                                 model = model, Z_ini = matrix(0,p,p), W_ini = matrix(0,p,p), eps_thre = eps.thre,
                                 eps_abs = 1e-5, eps_rel = 1e-3, max_iter = 50000, verbose = F))
      if(inherits(fit.L.theta1, "try-error")){
        fit.L.theta1 = NULL
        conv=FALSE
      }
      L.i = fit.L.theta1$L*fit.L.theta1$theta1
      conv.step.1.i = fit.L.theta1$conv
      print("step 1 complete")
      
      
      ## Step 2
      #if (i == 1){ ## To match with Dr. Peng's code.  Better results for R1 if we remove the if statement. 
      fit.eta = step.2(data, L = L.i, theta0 = theta0.i, eps = eps.thre)
      print("step 2 complete")
      
      ## Step 3
      eta.0 = fit.eta$eta[1]
      eta.1 = fit.eta$eta[2]
      R1.fit = step.3(eta.0, eta.1, L.i)
      
      print("step 3 complete")
      #}
      
      
      ## For eval purposes, store 1-pass estimates 
      if (i == 1){
        iniRes[[j]] = list("R1.ini" = R1.ini, "L.ini" = fit.L.theta1$L, "theta1" = fit.L.theta1$theta1, "eta.ini" = c(eta.0,eta.1), "S.ini" = S.i, "theta0.ini" = theta0.i, "conv.step1" = conv.step.1.i)
      }
      
    }
    time.end = proc.time()[3]
    time.vec[j] = time.end - time.start
    
    resList[[j]] = list("S" = S.i, "theta0" = theta0.i, "L" = fit.L.theta1$L, "theta1" = fit.L.theta1$theta1,
                        "eta" = c(eta.0, eta.1), "R1" = R1.fit, "conv.step1" = conv.step.1.i)
    print(paste("Pass", i, "complete"))
    R1.fit = R1.ini # Reset for other lambda
    
  }
  
  if(refit != 0){
    print("Starting refitting for L")
    
    ## Extract zero-patterns for refitting
    refit.zero = step.thre(resList, lambda.v, net.thre)
    
    ## Refit L estimates
    refit.L = step.lap.est(data, resList, refit.zero, lambda.v, net.thre, refit, model, eps.thre, eps_abs, eps_rel, max_iter, verbose, time.vec) # n-1 because that is sample size for innovations
  }
  
  ## End Time
  time.tot.end = proc.time()[3]
  time.total = time.tot.end - time.tot.beg
  
  
  return(list("result.pass" = resList, "ini" = iniRes, "refit" = refit.L, "A.net" = refit.zero, "time.total" = time.total))
  
}
