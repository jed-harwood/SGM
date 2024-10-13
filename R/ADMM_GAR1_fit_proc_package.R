## Goal: Create functions to apply GAR(1) fitting procedure, and model selection
## 2024/09/12

####################
## Required Packages 
####################

#require(SGM) # R Package that has our C++ functions
#require(gmp) # R package used for handling numerical issues
#require(foreach) # For running in Parallel
#require(doParallel) # For running in Parallel




##############
## Module 1: 3-Step Estimation Procedure for GAR(1) Model
########

################
### Step 0a: Get sample covariance matrix and theta0 estimate
##########

fit_step_0a = function(S, nobs){
  
  n = nobs
  
  sigma = eigen(S, symmetric = T, only.values = T)
  theta0.e = sqrt(1/max(sigma$values)) ## initial estimated theta0; for this to be consistent we need p=o(n)
  
  temp = list("S" = S, "theta0" = theta0.e, "n"=nobs)
  return(temp)
}

###############
### Step 1a: fit L given theta_0.e from step 0a and set v0=0: Separate algorithm
#############

fit_step_1a = function(step0a, lambda.v, rho.v, model, eps_thre, eps_abs, eps_rel, max_iter_1a, verbose){
  
  S = step0a$S
  theta0.e = step0a$theta0
  
  p = ncol(S)
  Z = matrix(0,p,p)
  W = Z
  phi=0
  
  result.L2.0 = vector("list", length(lambda.v))
  
  for (j in 1:length(lambda.v)){
    result.L2.0[[j]] = ADMM_L2(S, theta0.e, rep(0,ncol(S)), rho.v[j], lambda.v[j], model, Z, W, eps_thre, eps_abs, eps_rel, max_iter_1a, verbose)
  }
  
  return(result.L2.0)
  
}

#################
### Step 2a: refit L: given the pattern from Step 1a Sep results; set v0=0: Sep-Zero  algorithm
############

## Extract 0 pattern
fit_step_2a_1 = function(step1a, lambda.v, net.thre){
  
  ## Storage objects 
  conv.0.stat=rep(NA, length(lambda.v))
  A.0.net= vector("list",  length(lambda.v))     ## store the estimated network for post-selection estimation
  net.0.size=matrix(NA, length(lambda.v), length(net.thre))  # # of edges/2 in the estimated network == number of free non-zero parameters in the Laplacian
  
  
  ## Extract zero pattern
  for (j in 1:length(lambda.v)){
    result.c = step1a[[j]]
    
    if(!is.null(result.c)&&(conv.0.stat[j]=result.c$conv)==T){
      A.0.net[[j]]=vector("list", length(net.thre))
      
      for (k in 1:length(net.thre)){
        net.e = abs(result.c$L)>net.thre[k]
        diag(net.e) = 0
        net.0.size[j,k] = sum(net.e)/2
        A.0.net[[j]][[k]] = net.e
      }
    }
  }
  return(list("A.0.net"=A.0.net, "net.0.size" = net.0.size, "conv.0.stat" = conv.0.stat))
}

## refit L given zero patterns: still fix v0 at 0
fit_step_2a_2 = function(step0a, step2a.1, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_2a, verbose){
  
  ## Extract results from step0a and step 2a.1
  S = step0a$S
  theta0.e = step0a$theta0
  p = ncol(S)
  n = step0a$n
  result.0.post = vector("list", length(lambda.v))
  Z = matrix(0, p,p)
  W = Z
  
  A.0.net = step2a.1$A.0.net
  net.0.size = step2a.1$net.0.size
  conv.0.stat = step2a.1$conv.0.stat
  
  for(j in 1:length(lambda.v)){
    if(conv.0.stat[j]==T){
      result.0.post[[j]] = vector("list", length(net.thre))
      for(k in 1:length(net.thre)){
        result.0.post[[j]][[k]] = ADMM_L2_Zero(S, theta0.e, v=rep(0,p), rho=sqrt(log(p)/n), A=A.0.net[[j]][[k]], model, Z_ini=Z, W_ini = W, 
                                               eps_thre, eps_abs, eps_rel, max_iter_2a, verbose)
      }
    }
  }
  
  return(list("result.0.post" = result.0.post, "A.0.net" = A.0.net, "net.0.size" = net.0.size))
  
}

fit_step_2a = function(step0a, step1a, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_2a, verbose){
  
  ## Extract 0 pattern
  step_2a_1 = fit_step_2a_1(step1a, lambda.v, net.thre)
  
  ## Refit non-zero elements
  step_2a_2 = fit_step_2a_2(step0a, step_2a_1, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_2a, verbose)
  
  return(step_2a_2)
}

######################
## Step 3a.0: obtain v0 estimate from Step 2a L: 2024/6/25: Made as a separate function
#############

## Estimate degree vector
fit_step_3a_0 = function(step0a, step2a, lambda.v, net.thre, eps_abs, eps_rel, verbose){
  
  ## Initialize Objects
  matrix(NA, length(lambda.v), length(net.thre)) 
  v0.0.est=vector("list",  length(lambda.v))
  conv.0.v0=matrix(NA, length(lambda.v), length(net.thre)) 
  
  ## Extract result from previous steps
  result.0.post = step2a$result.0.post
  
  
  for (j in 1:length(lambda.v)){
    v0.0.est[[j]] = vector("list", length(net.thre))
    
    for (k in 1:length(net.thre)){
      result.c = result.0.post[[j]][[k]]
      
      if(!is.null(result.c)&&result.c$conv){
        L.est = result.c$L
        temp=try(ADMM.Deg.L(L.est,rho=0.1, epsilon=sqrt(1/(2*ncol(L.est))), eps.abs=1e-5, eps.rel=1e-3, max.iter=50000, verbose=FALSE))
        
        if (inherits(temp, "try-error")){
          temp = list("v" = rep(0,p), "conv" = F)
        }
        
        v0.0.est[[j]][[k]]=temp$v
        conv.0.v0[j,k]=temp$conv
      }
    }
  }
  
  return(list("v0.0.est" = v0.0.est, "conv.0.v0" = conv.0.v0))
  
}

## Estimate theta0 and L esimultaneously with 0-pattern from step1a, and estimated v0 from step 3a.1
fit_step_3a_1 = function(step0a, step2a, step3a.0, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_3a, verbose){
  
  ## Extract results from previous steps
  S = step0a$S
  n = step0a$n
  p = ncol(S)
  
  A.0.net = step2a$A.0.net
  v0.0.est = step3a.0$v0.0.est
  conv.0.v0 = step3a.0$conv.0.v0
  
  ## Initialize matrices
  result.0S=vector("list",  length(lambda.v))
  theta0.0S=matrix(NA, length(lambda.v), length(net.thre))
  Z = matrix(0,p,p)
  W = Z
  phi = 0
  
  for(j in 1:length(lambda.v)){
    result.0S[[j]] = vector("list", length(net.thre))
    
    for(k in 1:length(net.thre)){
      if(conv.0.v0[j,k] == T){
        v0.e = v0.0.est[[j]][[k]]
        result.0S[[j]][[k]] = ADMM_Lap_Zero(S, v0.e, rho=sqrt(log(p)/n), AA=A.0.net[[j]][[k]], model=model, ZZ_ini = Z, WW_ini = W, phi_ini = phi, 
                                            eps_thre = eps_thre, eps_abs = eps_abs, eps_rel = eps_rel, max_iter = max_iter_3a, Z_max_iter = 100000, 
                                            Z_conv_abs = 1e-5, Z_conv_rel = 1e-3, verbose = verbose)
        theta0.0S[j,k] = result.0S[[j]][[k]]$theta0
      }
    }
  }
  return(list("result.0S" = result.0S, "theta0.0S" = theta0.0S, "v0.0.est" = v0.0.est, "conv.0.v0" = conv.0.v0))
}

fit_step_3a = function(step0a, step2a, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_3a, verbose){
  
  
  ## Estimate degree vector
  step_3a_0 = fit_step_3a_0(step0a, step2a, lambda.v, net.thre, eps_abs, eps_rel, verbose)
  
  ## Estimate theta0 and L simultaneously given v0 and 0-patten
  step_3a_1 = fit_step_3a_1(step0a, step2a, step_3a_0, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_3a, verbose)
  
  
  return(step_3a_1)
}


#################################
## Model Fitting function
################################
#' GAR(1) fitting procedure
#' 
#' @description
#' `GAR1_fit` performs a three-step estimation procedure, using a penalized MLE approach, to estimate graph filter parameters `theta0` and `theta`, and the normalized graph Laplacian `L`.  
#' 
#' @param S An estimate of the covariance matrix, such as the MLE.
#' @param nobs The number of samples used to calculate `S`
#' @param lambda.v Tuning parameter to control sparsity of the estimated graph
#' @param net.thre Tuning parameter to control noisy entries in estimated graph
#' @param model
#' * "LN" Fits a normalized graph Laplacian
#' * "L" Fits a graph Laplacian
#' * "LN.noselfloop" Fits a normalized graph laplacian assuming no self-loops.
#' @param step How many steps of the estimation procedure you want to run. Either 2 or 3.
#' @param rho.v ADMM parameter (typically equal to `lambda.v`)
#' @param eps_thre Small positive number
#' @param eps_abs ADMM convergence criterion
#' @param eps_rel ADMM convergence criterion
#' @param max_iter_1, max_iter_2, max_iter_3a Maximum number of iterations for algorithm
#' 
#' @returns
#' A list object
#' * `S` The p by p estimated covariance matrix
#' * `result.L2.0` A list containing the Step 1a results
#' * `result.0.post` A list containing the Step 2a results
#' * `result.0S` A list containing the Step 3 results (NULL if `step<3`)
#' * `A.0.net` A matrix that encodes the estimated graph topology
#' * `net.0.size` An integer containing the number of edges in the estimated graph
#' * `v0.0.est` A p by 1 matrix containing the estimated degree vector
#' * `theta0.0` A positive number. The estimated theta0 from Step 2
#' * `theta0.0S` A positive number. The estimated theta0 from Step 3
#' * `conv.0.v0` A matrix containing convergence results for each combination of `lambda.v` (rows) and `net.thre` (columns)
#' 
#' @example man-roxygen/GAR1_fit_example.R
#' @export
GAR1_fit = function(S, nobs, lambda.v, net.thre, model, step = 3, rho.v=lambda.v, eps_thre=1e-6, eps_abs=1e-5, eps_rel=1e-3, max_iter_1a=10000, max_iter_2a = 10000, max_iter_3a = 10000, verbose=F){
  
  ## Get sample covariance and initial theta0 estimate
  step0a = fit_step_0a(S, nobs)
  print("Step 0a complete")
  
  ## Get 0-pattern form initial theta0 estimate and sample covariance
  step1a = fit_step_1a(step0a, lambda.v, rho.v, model, eps_thre, eps_abs, eps_rel, max_iter_1a, verbose)
  print("Step 1a complete")
  
  ## Using 0-pattern, re-estimate non-zero elements in L
  step2a = fit_step_2a(step0a, step1a, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_2a, verbose)
  print("Step 2a complete")
  
  if (step == 3){
    ## Use estimated L to estimate degree vector, then simultaneously re-estimate L and theta0
    step3a = fit_step_3a(step0a, step2a, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_3a, verbose)
    print("Step 3a complete")
  }
  else{
    step3a = NULL
  }
  
  resultList = list("S" = step0a$S, "result.L2.0" = step1a, "result.0.post" = step2a$result.0.post, "result.0S" = step3a$result.0S, 
                    "A.0.net" = step2a$A.0.net, "net.0.size" = step2a$net.0.size, "v0.0.est" = step3a$v0.0.est,
                    "theta0.e" = step0a$theta0, "theta0.0S" = step3a$theta0.0S, "conv.0.v0" = step3a$conv.0.v0)
  
  return(resultList)
}



##############
## Module 2: Fitting Process for GAR(1) Model
## Goal: Select model via eBIC and provide bootstrap goodness of fit
########

####################################
## eBIC and log-likelihood selection
######
#' Select tuning parameters for GAR(1) model
#' 
#' @description
#' Given a list of models, uses eBIC criterion to select the appropriate tuning parameters for GAR(1) and conducts a goodness of fit test.
#' 
#' @param resultList A list output from `GAR1_fit`
#' @param n An integer referring to the number of observations
#' @param step 2 or 3; How many steps were used to fit the model.
#' @param model Which model to consider
#' * "LN" Normalized graph Laplacian
#' * "L" Graph Laplacian
#' * "LN.noselfloop" Normalized graph Laplacian without self-loops. 
#' 
#' @returns A list object
#' * model.selec A list containing the model with the model chosen by the eBIC criterion.
#' * A.net.e A matrix encoding the (unweighted) graph chosen by the eBIC criterion.
#' * index Index for the optimal tuning parameters (lambda, net.thre) for the eBIC-selected model.
#' * ebic ebic score for the selected model
#' 
#' @export
model_selec = function(resultList, n, step = 3, model = "LN"){
  
  ## extract p, n, and results
  S = resultList$S
  p = nrow(S)
  
  theta0.ini = resultList$theta0.e # Step 0: Initial Theta estimate   
  
  result.0 = resultList$result.L2.0 # Step 1a: Zero-patterns
  A.0.net = resultList$A.0.net # Step 1a: Estimated Networks (Step 1a)
  net.0.size = resultList$net.0.size # Step 1a: Estimated network size
  
  result.post = resultList$result.0.post # Step 2a
  
  result.0S = resultList$result.0S # Step 3a: L, phi, etc.
  theta0.0S = resultList$theta0.0S # Step 3a: theta0 est. 
  v0.0.est = resultList$v0.0.est # Step 3a: Est. v0 
  
  ## Number of lambda and net.thre (tuning parameters)
  n.lambda = dim(resultList$theta0.0S)[1]
  n.net.thr = dim(resultList$theta0.0S)[2]
  
  ## set gamma for ebic
  if(p/n>0.5){## e.g., for p=100,n=100
    gamma=1  ## eBIC parameter: p/n~1 set as 1; when p/n<0.5: set as 0.5 
  }else{
    gamma=0.5
  }
  
  ## Number of edges possible 
  P.total=p*(p-1)/2
  
  ## Convergence status 
  conv.post=array(NA, dim= dim(resultList$theta0.0S)) ## Step 2a + Step3a.0 convergence 
  conv.0S=conv.post ## Step 3a.1 convergence 
  
  ## log-likelihood and ebic
  log.post.like=conv.post  ## Step 2a
  log.0S.like=conv.post ## Stap 3a.1
  
  bic.post=conv.post ## Step 2a
  bic.0S=conv.post ## Step 3a.1
  
  ebic.post=conv.post ## Step 2a
  ebic.0S=conv.post ## Stap 3a.1
  
  ## Calculate BIC, eBIC, and log-likelihood
  for (j in 1:n.lambda){
    for (k in 1:n.net.thr){
      
      ## Extract network information
      net.size.c = net.0.size[j,k] # Extract network size
      A.0.net.c = A.0.net[[j]][[k]] # Extract Adjacency matrix
      
      ## Calculate eBIC term
      ebic.term=2*gamma*(lfactorial(P.total)-lfactorial(net.size.c)-lfactorial(P.total-net.size.c))
      
      ## If stopping at step 2a, calculate BIC, eBIC, and log-likelihood
      if(step == 2){
        result.c = result.post[[j]][[k]]
        if(!is.null(result.c) && (conv.post[j,k]=result.c$conv)==TRUE){
          ## Extract the estimates for L and theta1 
          L.est=result.c$L
          theta1.e=result.c$theta1
          
          ## Calculate log-likelihood
          log.post.like[j,k]=LogLike(S, theta0.ini, theta1.e, L.est, n)
          
          ## Calculate BIC and eBIC
          if(model == "LN"||model=="LN.noloop"){## LN models have additional p parameters accounting for the estimated v0
            bic.post[j,k]=BIC(log.post.like[j,k], n, net.size.c+1+p)
          }else{
            bic.post[j,k]=BIC(log.post.like[j,k], n, net.size.c+1)
          }
          ebic.post[j,k]=bic.post[j,k]+ebic.term
          
        }
        
      }
      
      ## If stopping at step 3a, calculate BIC, eBIC, and log-likelihood
      if(step==3){
        result.c = result.0S[[j]][[k]]
        if(!is.null(result.c)&&(conv.0S[j,k]=result.c$conv)==TRUE){
          
          ## Extract estimates for theta0, v0, L, and theta1
          theta0.e = theta0.0S[j,k]
          v0.e = v0.0.est[[j]][[k]]
          v0.e = v0.e/sqrt(sum(v0.e^2)) ## make norm=1
          L.est = result.c$L
          theta1.e = result.c$theta1
          
          ## Calculate log-likelihood
          log.0S.like[j,k]=LogLike(S, theta0.e, theta1.e, L.est, n)
          
          ## Calculate BIC and eBIC
          if(model == "LN"||model=="LN.noloop"){## LN models have additional p parameters accounting for the estimated v0, but we do not count it
            bic.0S[j,k]=BIC(log.0S.like[j,k], n, net.size.c+1)
          }else{
            bic.0S[j,k]=BIC(log.0S.like[j,k], n, net.size.c+1)
          }
          
          ebic.0S[j,k]=bic.0S[j,k]+ebic.term
        }
      }
    }
  }
  
  ######################
  ## Select the optimal lambda, net.thre, according to desired threshold
  ####################
  
  
  ## Select optimal index, depending on which steps used
  if(step == 2){
    ## Step 2a: post estimator selection results 
    index.c=which.min(ebic.post)
    index.c = arrayInd(index.c, .dim = dim(ebic.post))
    
    resultOptimal = result.post[[index.c[1]]][[index.c[2]]]
    A.0.net.opt = A.0.net[[index.c[1]]][[index.c[2]]]
    ebic.opt = ebic.post[index.c]
  }
    
  if(step == 3){
    ## Step 3a: 0S estimator selection results 
    index.c = which.min(ebic.0S)
    index.c = arrayInd(index.c, .dim = dim(ebic.0S))
    resultOptimal = result.0S[[index.c[1]]][[index.c[2]]]
    A.0.net.opt = A.0.net[[index.c[1]]][[index.c[2]]]
    ebic.opt = ebic.0S[index.c]
  }
  
  
  ############################
  ## Return the optimal lambda, net.thre, and goodness of fit test
  ###########################
  colnames(index.c) = c("lambda", "net.thre")
  rownames(index.c) = "index"
  retList = list("model.selec" = resultOptimal, "A.net.e" = A.0.net.opt, "ebic" = ebic.opt, "index" = index.c)
  return(retList)
}
