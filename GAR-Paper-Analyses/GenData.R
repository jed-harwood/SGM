##11/15/2018
##Jie: generating data based on spectral graphical model
## 04/2023: likelihood, BIC; TAR-GAR(1) generation 
## 2024/5/12: updated the Laplacian.Norm<-function(A)

############################
##### generate graph
###########################


##### generate random adjacency matrix:  ### 2024/5/12: add a self-loop probability argument  
##### generate edge connection according to a prob; generate weights uniformly from [min,max];
Rand.Graph<-function(p,edge.prob, self.prob=edge.prob, min=0.1, max=1, selfloop=FALSE, isolate=FALSE){
  ## p: number of nodes 
  ## edge.prob: probability of edge connection
  ## min: minimum nonzero weights, max : maximum weights 
  ## selfloop: if false, then diagonal must be zero
  ## isolate: if false, then each node has to have at least one edge 
  ## return: an adjacent matrix: A = t(A), all entries >=0
  A = matrix(0, p,p)
  A[upper.tri(A)] = runif(p*(p-1)/2, min=min, max=max)*sample(0:1,p*(p-1)/2,replace=TRUE,prob=c(1-edge.prob, edge.prob))
  A= A+t(A)
  
  if(selfloop){
    diag(A)=runif(p, min=min, max=max)*sample(0:1,p,replace=TRUE,prob=c(1-self.prob, self.prob))
  }
  
  if(!isolate){
    deg=apply(A>0, 1, sum)
    for (i in 1:(p-1)){
      if(deg[i]==0){
        j=sample((i+1):p,1)
        A[i,j]=runif(1, min=min, max=max)
        A[j,i]=A[i,j]
        deg[j]=1
      }
    }
  }
  
  return(A)
}


###### Random power-law graph:
PowerLaw.Graph<-function(p,power=1, min=0.1, max=1, edge.prob=1/p, selfloop=FALSE){
  ## p: number of nodes 
  ## power: power parameter for power-law degree distribution
  ## min: minimum nonzero weights, max : maximum weights 
  ## edge.prob: small probability for random connection added on top of the power-law net to add a bit flexibility 
  ## selfloop: if false, then diagonal must be zero
  ## return: an adjacent matrix: A = t(A), all entries >=0
  
  library(igraph)
  temp<-sample_pa(n=p, power=power, directed=FALSE)  
  temp.ends=ends(temp, E(temp))
  
  adj=matrix(0, p, p)
  for (i in 1: nrow(temp.ends)){
    adj[temp.ends[i,1], temp.ends[i,2]]<-1
    adj[temp.ends[i,2], temp.ends[i,1]]<-1
  }
  
  ##jie notes: 11/12/2018: the above graph will always have p-1 edges; and each node has at least one edge, there is no self-loop
  ## add another round of edges: so edge number does not have to be p-1
  for (i in 1:p){
    temp.new=sample(0:1,1,prob=c(1-edge.prob, edge.prob)) ## will be a new edge for node i? 
    if(temp.new==1){##if so 
      j=sample(1:p,1)
      adj[i,j]=1
      adj[j,i]=1
    }
  }
  
  if(!selfloop){##if no selfloop
    diag(adj)=0
  }
  
  ## generate adjacent matrix A 
  for (i in 1:p){
    for (j in 1:i){
      if(adj[i,j]>0){
        adj[i,j]=runif(1, min=min, max=max)
        adj[j,i]=adj[i,j]
      }
    }
  }
  
  #  
  return(adj)
}


###########################
#### Calculate Laplacian given graph 
###########################

#### graph Laplacian
Laplacian<-function(A){
  ## A: p by p adjacency matrix
  ## return: Laplacian: L = D-A 
  
  deg=apply(A, 1, sum)
  L=diag(deg)-A
  return(L)
  
}



#### Normalized graph Laplacian
### 2024/5/12: fix an error (only matters for the case with self-loop)
Laplacian.Norm<-function(A){
  ## A: p by p adjacency matrix
  ## return: normalized Laplacian: L = I-D^{-1/2}AD^{-1/2} 
  p=nrow(A)
  deg=apply(A, 1, sum)
  
  L.norm=matrix(0,p,p)
  
  for (i in 1:p){
    if(deg[i]>0){
      L.norm[i,i]=1-A[i,i]/deg[i]
    }
    
    for (j in 1:p){
      if((j!=i)&&deg[j]>0&&A[i,j]>0) L.norm[i,j]=-A[i,j]/sqrt(deg[i]*deg[j]) ## updated: 2024/5/12
    }
    
  }
  
  
  return(L.norm)
}

isNormLap<-function(LN){
  ## return: True if is a noarmlized Laplacian
  d=eigen(LN, symmetric = TRUE)$values
  return(isSymmetric(LN)&all(LN[upper.tri(LN)]<=0)&all(LN[upper.tri(LN)]>=-1)&all(diag(LN)<=1)&all(d>=0)&(min(d)<1e-6))
}


###########################
##### Generate data under GAR(1) model 
###########################

###### GAR(1) model: Sigma, Omega  given L 
GGM.GAR<-function(L, theta0, theta1){
  ## L: p by p Laplacian matrix  
  ## return: the corresponding Omega and Sigma under GAR(1) model 
  p=nrow(L)
  Omega=(theta0*diag(1,p)+theta1*L)%*%(theta0*diag(1,p)+theta1*L)
  Sigma=solve(Omega)
  return(list(Sigma=Sigma, Omega=Omega))
  
}


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



#############################
########### generate data under TAR-GAR(1) model  
#############################
###### order 0 and 1 population auto-covariance matrix of TAR-GAR(1) model given eta, theta, L
###
AutoCov01.TAR.GAR<-function(eta0,eta1, theta0, theta1, L){
  ## L: p by p Laplacian matrix; Outcome of the Laplacian() function
  ## theta0, theta1 : nonnegative scaler 
  ## eta0, eta1: filter parameters 
  p=nrow(L) ##dimension
  
  ## check whether eta0, eta1 satisfies the stationarity condition 
  temp=eigen(L, symmetric=TRUE)
  lambda=temp$values 
  
  if(!all(abs(eta0+eta1*lambda)<1)){
    stop("|eta0+eta1 lambda_j| not all less than 1: process not stationary!")
  }
  
  
  
  Q=t(temp$vectors) ## L= t(Q)%*% Lambda%*%Q
  
  ## Sigma=var(U.t): innovation covariance according to GAR(1) model 
  Omega=t(Q)%*%diag((theta0+theta1*lambda)^2)%*%Q
  Sigma=t(Q)%*%diag((theta0+theta1*lambda)^{-2})%*%Q
  
  ## R1= eta0+eta1*L: TAR filter 
  R1=diag(eta0,p)+eta1*L
  
  ## Calculate Gamma0=var(Y.t)=Sigma%*%solve(I-R1%*%R1): order zero covariance of the process 
  Gamma0= t(Q)%*%diag((theta0+theta1*lambda)^{-2}*(1-(eta0+eta1*lambda)^2)^{-1})%*%Q
  
  ## Gamma1=R1%*%Gamma0
  Gamma1= t(Q)%*%diag((theta0+theta1*lambda)^{-2}*(1-(eta0+eta1*lambda)^2)^{-1}*(eta0+eta1*lambda))%*%Q
  ##
  return(list(Gamma0=Gamma0, Gamma1=Gamma1, Omega=Omega, Sigma=Sigma, R1=R1))
}

###
GenData.TAR.GAR<-function(n,eta0,eta1, theta0,theta1, L, rep=1){
  ## n: sample size
  ## eta0, eta1: TAR filter parameters: eta0+eta1 L 
  ## theta0>0, theta1>0: GAR(1) parameters 
  ##L: (normalized)  Laplacian, p by p matrix 
  ## rep: number of replicates 
  ## return: list of n by p data matrices; p by p concentration matrix: Omega; p by p covariance matrix: Sigma
  
  library(mnormt)
  p=nrow(L) ##dimension
  data.rep=NULL
  
  ##
  temp=AutoCov01.TAR.GAR(eta0,eta1, theta0, theta1, L)
  Gamma0=temp$Gamma0  ##variance of the process: (Y.t)
  Sigma=temp$Sigma    ## variance of the innovation process: var(U.t)
  R1=temp$R1          ## TAR filter: Y.t=R1%*%Y_{t-1}+ U.t 
  
  ##
  for (i in 1:rep){##ith replicate
    
    ## generate innovation process U.t
    U.c=rmnorm(n, mean=rep(0,p), varcov=Sigma)
    
    ## generate process Y.t
    data.c=matrix(0, n,p)
    for (t in 1:n){
      if(t==1){
        ## Y.1 ~ N(0, Gamma0)  
        data.c[1,]=rmnorm(1, mean=rep(0,p), varcov=Gamma0)
      }else{
        ## Y.t= R1%*%Y.{t-1} + U.t
        data.c[t,]=data.c[t-1,]%*%R1+U.c[t,]  
      }
    }##end t-loop
    
    data.rep[[i]]=data.c
  }## end i-loop
  
  return(data.rep)  
}


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
  temp2= determinant(temp, logarithm = T)$modulus[[1]] #previously: log(det(temp))
  
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


### 
LogLike.glasso<-function(S, Omega, n){
  ##S: p by p sample covariance matrix  
  ###Omega: p by p concentration matrix, "wi" component of glasso returned value 
  p=nrow(S)
  
  temp1=sum(diag(S%*%Omega))
  temp2 = determinant(Omega, logarithm = TRUE)$modulus[[1]]
  #temp2=log(det(Omega))
  
  loglike=-n/2*temp1+n/2*temp2-n*p/2*log(2*pi)
  
  return(loglike)
  
}

### create CV folds
CV.create<-function(data, K){
  ## inputs: data: n by p date matrix 
  ## K: cv folds 
  ## return: index.train, index.test: list of indices for  training and testing data of each CV fold
  ## data.train, data.test: list of data matrices for trainign and testing of each CV fold 
  ## test indices are every K observations 
  
  index.train=vector("list",K) 
  index.test=vector("list",K)  ##every Kth 
  data.train=vector("list",K) 
  data.test=vector("list",K) 
  
  n=nrow(data)
  
  for (j in 1:K){
    index.test[[j]]=seq(from=j, to=n, by=K)
    index.train[[j]]=setdiff(1:n, index.test[[j]])
    
    data.test[[j]]=data[index.test[[j]], ]
    data.train[[j]]=data[index.train[[j]],]
  }
  
  return(list(index.train=index.train, index.test=index.test, data.train=data.train, data.test=data.test))
  
}

### calculate CV score based on likelihood of test data for GAR model 
CV.like.score<-function(data.test, theta0, theta1 = 1, L){
  ## data.test: test data set 
  ## theta0, theta 1: GAR(1) parameters; theta1:  for the Laplacian model, theta1=1;
  ## L: p by p Laplacian matrix (normalized or not)
  ## return: log-likelihood on data.test divided by sample size n  
  
  n=nrow(data.test)
  p=ncol(data.test)
  S=t(data.test)%*%data.test/n
  
  temp=diag(theta0, p)+theta1*L
  
  temp1=sum(diag(S%*%temp%*%temp))
  temp2=log(det(temp))
  
  #loglike=-n/2*temp1+n*temp2-n*p/2*log(2*pi)
  loglike=-1/2*temp1+temp2-p/2*log(2*pi)
  return(loglike)
}


### 
CV.like.glasso<-function(data.test, Omega){
  ## data.test: test data set 
  ###Omega: p by p concentration matrix, "wi" component of glasso returned value 
  ## return: log-likelihood on data.test divided by sample size n  
  
  n=nrow(data.test)
  p=ncol(data.test)
  S=t(data.test)%*%data.test/n
  
  temp1=sum(diag(S%*%Omega))
  temp2=log(det(Omega))
  
  #loglike=-n/2*temp1+n/2*temp2-n*p/2*log(2*pi)
  loglike=-1/2*temp1+1/2*temp2-p/2*log(2*pi)
  
  return(loglike)
  
}


### added on 5/12/2024
###### evaluation and summary functions
########evaluation function: error measures 
EvalErr<-function(admm.res, net.thre, L.tr, theta0.tr, theta1.tr, net.tr, err.v0,rep=100){
  
  #
  p=nrow(L.tr)
  
  #
  theta0.err=numeric(rep)+NA
  L.err=numeric(rep)+NA
  
  Omega.err=L.err
  KL.div=L.err
  
  fdr=L.err
  power=L.err
  
  conv=L.err
  
  #
  temp=GGM.GAR(L.tr, theta0.tr, theta1.tr)
  Omega.tr=temp$Omega
  Sigma.tr=temp$Sigma
  
  #
  for (i in 1:rep){
    
    if(!err.v0[i]){
      result.c=admm.res[[i]]
      
      conv[i] = FALSE
      
      if(!is.null(result.c)&&result.c$conv){
        
        conv[i]=TRUE
        
        L.est=(result.c$L) 
        theta0.e=result.c$theta0
        theta1.e=result.c$theta1
        Omega.e= GGM.GAR(L.est, theta0.e, theta1.e)$Omega
        
        ##
        theta0.err[i]=(theta0.e-theta0.tr)^2
        
        ##estimation accuracy of theta1*L
        L.err[i]=sum((theta1.e*L.est-theta1.tr*L.tr)^2)/sum((theta1.tr*L.tr)^2)
        
        #
        Omega.err[i]=sum((Omega.e-Omega.tr)^2)/sum(Omega.tr^2)
        KL.div[i]=0.5*(-log(det(Omega.e))-log(det(Sigma.tr))-p+sum(Omega.e*Sigma.tr))
        
        ##graph inference 
        #estimated latent graph
        if(net.thre==0){## for fixed pattern estimators use the pattern preserving W
          net.e=abs(result.c$W)>net.thre
        }else{## for other estimators, use the estimated L
          net.e=abs(L.est)>net.thre  
        }
        
        diag(net.e)=0
        power[i]=sum(net.e*net.tr)/sum(net.tr) ## power 
        fdr[i]=sum(net.e*(1-net.tr))/sum(net.e) ##fdr
      }## end of if 
      
    }
  }##end of i loop
  
  return(list(theta0=theta0.err, L=L.err, Omega=Omega.err, KL=KL.div, fdr=fdr, power=power, conv=conv))
}




##### summaries of error measures 

EvalSum<-function(err.list){
  
  
  conv_perc<- mean(err.list$conv, na.rm=TRUE)
  
  # theta0
  theta0_l2<-round(mean(err.list$theta0,na.rm=TRUE),4)
  
  
  # L
  L_l2=round(mean(err.list$L,na.rm=TRUE),4)
  
  # Omega and KL
  Omega_l2=round(mean(err.list$Omega,na.rm=TRUE),4)
  KL=round(mean(err.list$KL,na.rm=TRUE),3)
  
  # FDR and power 
  fdr=round(mean(err.list$fdr, na.rm=TRUE),4)
  power=round(mean(err.list$power,na.rm=TRUE),4)
  
  #
  return(c("conv_perc"=conv_perc, "theta0_l2"=theta0_l2, "L_l2"=L_l2, "Omega_l2"=Omega_l2, "KL"=KL, "fdr"=fdr, "power"=power))
}
