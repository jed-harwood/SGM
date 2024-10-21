###################################################
## Goal: ADMM Algorithm to Estimate Degree Vector From a Laplacian Matrix
## Date: 2024-09-12
#####################################################



#################################
## Update Functions for Iterates
##############################

##### v-minimization
v.Min<-function(Q, Lambda, lambda.s, rho, w,u, verbose=TRUE){
  ## (Q,Lambda): eigen-decomposition of S = Q^T Lambda Q
  ## lambda.s: leading eigenvalue of S 
  ## w: p by 1: v-w=0
  ## u: p by 1, dual variable 
  ## rho: ADMM parameter 
  ## return: updated v
  p=nrow(Q)
  
  tk=matrix(rho*(w-u),p,1)
  tk.tilde=Q%*%tk
  
  norm.mu<-function(mu){
    Lambda.tilde= (Lambda-lambda.s)^2+2*mu 
    Lambda.inv=1/Lambda.tilde
    val = sum((Lambda.inv*tk.tilde)^2)-1
    return(val)
  }
  
  temp=uniroot(norm.mu, interval=c(1e-6, 1e+6))
  if(verbose){
    print(paste("v step root:", temp$root, temp$f.root))
  }
  mu.s = temp$root 
  
  Lambda.tilde.s= (Lambda-lambda.s)^2+2*mu.s
  Lambda.inv.s = 1/Lambda.tilde.s
  v=t(Q)%*%matrix((Lambda.inv.s*tk.tilde), p,1)
  
  return(v)
  
  
}


##### w-minimization
w.Min<-function(v,u, epsilon=1e-2){
  ##v: p by 1, v-w=0
  ##u: p by 1, dual variable 
  ##epsilon: positive small value 
  ## return: updated w
  
  temp=v+u
  w=temp
  w[temp<epsilon] = epsilon
  
  return(w)
}






########################
## WRAPPER FUNCTION ####
########################

#' Estimate degree vector from a given Laplacian matrix
#' 
#' @description
#' `ADMM.Deg.L` returns an estimated degree vector to minimize ||Lv0||, with each element restricted to be at least epsilon.
#' 
#' @returns
#' * `v`: A p by 1 matrix
#' * `w`: A p by 1 matrix
#' * `deg`: A p by 1 matrix
#' * `conv`: A boolean indicating whether the algorithm converged (TRUE) or didn't (FALSE).  
#' 
#' @param L A p by p normalized graph Laplacian
#' @param rho ADMM parameter (a positive number)
#' @param epsilon A positive small value
#' @param eps.abs ADMM convergence criterion
#' @param eps.rel ADMM convergence criterion
#' @param max.iter Maximum number of iterations to run the algorithm
#' @param verbose Trace of fitting procedure
#' 
#'@export
ADMM.Deg.L<-function(L,rho=0.01, epsilon=sqrt(1/nrow(S)), eps.abs=1e-5, eps.rel=1e-3, max.iter=1000, verbose=TRUE){
  ## ADMM algorithm to find the 0-eigenvector of LN, proportional to square-root of the degree vector of the graph
  ## L: p by p, Laplacian 
  ## rho: ADMM parameter 
  ## epsilon: positive small value
  ## return: v0: p by 1 matrix to minimize ||Lv0||, unit norm and  each element restricted to be at least epsilon
  ## deg: v0^2, proportional to degree vector, each element to be positive 
  ## conv: convergence status 
  
  p=nrow(L)
  
  temp = eigen(L, symmetric = TRUE)
  Q = t(temp$vectors)
  Lambda=temp$values
  lambda.s = 0
  
  w= rep(1/sqrt(p),p) ## initial values 
  u=numeric(p)
  
  conv=FALSE
  iter=0
  while((!conv)&&iter<=max.iter){
    v.up = v.Min(Q, Lambda, lambda.s, rho, w,u, verbose)
    w.up = w.Min(v.up,u, epsilon)
    u.up=u+(v.up-w.up) ## dual update
    
    res.prim=sqrt(sum((v.up-w.up)^2)) ##primal and dual residuals 
    res.dual=rho*sqrt(sum((w-w.up)^2))
    
    eps.prim=sqrt(2*p)*eps.abs+eps.rel*max(sqrt(sum(v.up^2)), sqrt(sum(w.up^2)))
    eps.dual=sqrt(p)*eps.abs+eps.rel*rho*sqrt(sum((u.up)^2))
    
    if(verbose){
      print(paste("Primal:", res.prim, eps.prim))
      print(paste("Dual:", res.dual, eps.dual))
    }
    
    if(res.prim<=eps.prim&&res.dual<=eps.dual){## check convergence 
      conv = TRUE
    }
    
    iter=iter+1
    
    v=v.up  ##make updates
    w=w.up
    u=u.up
  }
  
  return(list(v=v, w=w, deg=v^2, conv=conv))
}