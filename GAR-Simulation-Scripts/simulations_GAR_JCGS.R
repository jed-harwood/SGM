##############################
## Simulation Script for GAR:
## - Dimension: p
## - Sample size: n
## - Edge probability: edge.prob
## - rep: Number of samples to simulate
## - model: Use as LN; Can also use L for the (unnormalized) graph Laplacian
## - Includes the GF Measure in the results
##############################


## Clear the environment
rm(list = ls())

## Call the necessary packages
library(SGM)
library(glasso)
library(doParallel)
library(doFuture)
library(future)
library(foreach)
library(mnormt)


## Load the Data Generation Script
source("GenData.R")


#####################################
#### Parallel Thread Setup
####################################
num.thread=25
registerDoFuture()
plan(multisession, workers = num.thread)
options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB


#####################################
#### Model Setup & Data Generation
#####################################
p=100 #
n=100 # 
rep=100 #
model="LN" ##fit the general normalized Laplacian model 
theta0=1
theta1=2
filename=paste("sep_GAR1_",model, "_p", p, "_n",n,"_seqlamnet.Rdata", sep="")

#########
#########  generate adjacent matrix and Laplacians
set.seed(1)
edge.prob=2/p
A.tr=Rand.Graph(p=p,edge.prob=edge.prob, self.prob=2*edge.prob, min=0.5, max=1, selfloop=FALSE, isolate=FALSE)
net.tr=(A.tr>0)
diag(net.tr)=0
deg=apply(A.tr,1,sum)
summary(deg)
L=Laplacian(A.tr)  ##Laplacian
LN=Laplacian.Norm(A.tr) ##normalized Laplacian

#### Specify truth given model
if(model == "LN"||model=="LN.noloop"){
  L.tr=LN
}else{
  L.tr=L
}

theta0.tr=theta0
theta1.tr=theta1

######network summaries
summary(diag(A.tr)) ## all zero: no self-loop; otherwise could have non-zero diagonals 
summary(diag(LN)) ## diagonals of LN: all 1 as if self-loop; all<=1 for general case 
summary(LN[upper.tri(LN)]) ## off-diagonals of LN: all <=0  
sum(net.tr)/2  ## p=100: 105; p=250: 307; p=500: 555

###### the true v0 vector:
if(model=="LN"||model=="LN.noloop"){ ## normalized Laplacian models 
  v0=deg^{0.5}
  v0=matrix(v0/sqrt(sum(v0^2)),p,1)
}else{ ## Laplacian model 
  v0=rep(1,p)
}
v0.tr = v0
summary(v0.tr)

###### Draw Samples from the GAR Model:
set.seed(5) # For reproducability 

######## Determine whether to use no self loops, normalized Laplacian, or unnormalized Laplacian
if(model == "LN"||model=="LN.noloop"){## normalized Laplacian models
  temp=GenData.L2(n,theta0,theta1,LN, rep)
}else{ ## Laplacian model 
  temp=GenData.L2(n,theta0,theta1,L,rep)
}

######## Extract the replicates, true covariance, and the true Gaussian Graphical Model
data=temp$data
Omega.tr=temp$Omega
Sigma.tr=temp$Sigma
A.ggm.tr = abs(Omega.tr) > 1e-6
diag(A.ggm.tr) = 0
sum(A.ggm.tr)/2



######################
## Fit the Models (GAR)
######################


##### GAR Tuning Parameters

## Lambda Tuning parameter
C.v=c(1,0.5)  # decreasing seq:
lambda.v=C.v*sqrt(log(p)/n) 

### ADMM parameter 
rho.v=pmax(lambda.v, 0.01)

## Epsilon_{thre} Tuning Parameter
if(p==100){
  C.thre=exp(seq(log(1),log(0.05), length.out=10)) ##p=100
}else if(p==250){
  C.thre=exp(seq(log(1),log(0.075), length.out=10))  ##p=250
}else{
  C.thre=exp(seq(log(1),log(0.1), length.out=10))  ##p=500
}
net.thre=C.thre*sqrt(log(p)/n)  ##threshold for the estimated adjacency matrix entries 





########## Begin Model Training for GAR

results.GAR = foreach(i = 1:rep, .maxcombine=max(rep,2))%dopar%{ 
  ## Status message
  print(paste("Fitting GAR for replicate: ", i))
  
  ## Fit the GAR model to the data
  S.i = var(data[[i]])*((n-1)/n) # Sufficent Statistic for Sigma
  GAR.i.res = SGM::GAR1_fit(S = S.i, nobs = n, lambda.v = lambda.v, net.thre = net.thre,
                            model = model, rho.v = rho.v, eps_thre = 1e-6, eps_abs = 1e-5,
                            eps_rel = 1e-3, max_iter_3a = 10000, verbose=FALSE)
  
  ## Store results important results in a list
  res.i = list("S" = S.i, "conv" = GAR.i.res$conv, "step3b" = GAR.i.res$step3b,
               "modelList" = GAR.i.res, "A.net" = GAR.i.res$A.net,
               "step1" = GAR.i.res$step1)
  
  ## Return
  res.i
}


###################################################
## Model Selection and Evaluation Metrics  (GAR) ##
###################################################

## GAR and Graph Recovery Metrics
GAR.ebic = vector(mode = "list", length = rep) # Storage list for eBIC-chosen GAR models
L.ebic.err = rep(NA, rep) # Storage vector for L Errors
theta0.ebic.err = rep(NA, rep) # Storage vector for theta_0 Errors
v0.ebic.err = rep(NA, rep) # Storage vector for v0 Errors
power.ebic.vec = rep(NA,rep) # Storage vector for Power
fdr.ebic.vec = rep(NA,rep) # Storage vector for FDR
F1.ebic.vec = rep(NA,rep) # Storage vector for F1-Score
L.s1.err = rep(NA, rep) # Vector for Step 1 estimator of L, C_lambda = 1
theta0.s1.err = rep(NA, rep) # Vector for Step 1 estimator of theta0, from Step 1

## GGM-induced-by-GAR Recovery Metrics
Sigma.gar.err = rep(NA,rep) ## Error for Covariance Matrix
Omega.gar.err = rep(NA,rep) ## Error for Precision Matrix
power.ggm.gar = rep(NA,rep) ## Power for GGM
FDR.ggm.gar = rep(NA, rep)  ## FDR for GGM
F1.ggm.gar = rep(NA,rep)    ## F1 for GGM


## For each replicate, find the eBIC selected model
for (i in 1:rep){
  ## Extract the GAR models for ith sample
  GAR.models.i = results.GAR[[i]]$modelList
  
  ## Begin the model selection via eBIC
  GAR.ebic.i = SGM::model_selec(resultList = GAR.models.i)
  GAR.ebic[[i]] = GAR.ebic.i$selected.model
  
  ## Extract the estimated parameters
  A.ebic.i = GAR.ebic.i$A.net.e # The 0-1 adjacency matrix for the eBIC selected model
  L.ebic.i = GAR.ebic.i$L * GAR.ebic.i$theta1 # L * theta_1
  theta0.ebic.i = GAR.ebic.i$theta0 # theta0 for the eBIC selected model
  v0.ebic.i = GAR.ebic.i$v0 # v0 from eBIC selected model 
  v0.ebic.i = v0.ebic.i/sqrt(sum(v0.ebic.i^2)) # Renormalize (just-in-case) 
  
  ## Record Error and Evaluation Metrics
  #### GAR Graph Metrics
  net.size.ebic = sum(A.ebic.i > 0)/2 # Network size
  L.ebic.err[i] = sum((L.ebic.i-L.tr*theta1.tr)^2)/sum((theta1.tr*L.tr)^2) # L-Error
  theta0.ebic.err[i] = abs(theta0.ebic.i - theta0.tr)^2 # theta0-error
  v0.ebic.err[i] = sum(abs(v0.ebic.i - v0.tr)^2) # v0-error
  #print(paste("The V0s are equal:", all.equal(v0.ebic.i, v0.tr)))
  power.ebic.vec[i] = sum(A.ebic.i*net.tr)/sum(net.tr) # power
  fdr.ebic.vec[i] = sum(A.ebic.i*(1-net.tr))/sum(A.ebic.i) # fdr
  F1.ebic.vec[i] = (2*(1-fdr.ebic.vec[i])*power.ebic.vec[i])/(1-fdr.ebic.vec[i] + power.ebic.vec[i]) # F1 Score
  
  
  ## Record the Step 1 Estimates for Comparison with Oracle Estimator
  L.gar.s1 = results.GAR[[i]]$step1[[1]]$L * results.GAR[[i]]$step1[[1]]$theta1
  theta0.gar.s1 = results.GAR[[i]]$step1[[1]]$theta0
  
  theta0.s1.err[i] = abs(theta0.gar.s1 - theta0.tr)^2 # Step 1 theta0 Error
  L.s1.err[i] = sum((L.gar.s1 - L.tr*theta1.tr)^2)/sum((theta1.tr*L.tr)^2) # Step 1 L error
  
  ## GGM Graph Recovery using GAR
  ### Metrics for GAR-estimated GGM
  temp.ggm = GGM.GAR(GAR.ebic.i$L, 
                     theta0.ebic.i, 
                     GAR.ebic.i$theta1)
  
  Sigma.gar.err[i] = sum( (temp.ggm$Sigma - Sigma.tr)^2 )/sum(Sigma.tr^2) # Cov
  Omega.gar.err[i] = sum( (temp.ggm$Omega - Omega.tr)^2 )/sum(Omega.tr^2) # Precision
  
  ## Use -W for the Laplacian estimate (thresholded entries of 0)
  temp.L = (-1) * GAR.ebic.i$selected.model$W
  
  ## Calculate Omega and Sigma using GAR Model
  temp.ggm = GGM.GAR(temp.L, theta0.ebic.i, GAR.ebic.i$theta1) 
  
  ## Adjacency Matrix for GGM Induced by GAR
  A.ggm = abs(temp.ggm$Omega) > 1e-6
  diag(A.ggm) = 0
  
  ## GGM Graph Recovery Metrics
  FDR.ggm.gar[i]= sum(A.ggm*(1-A.ggm.tr))/sum(A.ggm)
  power.ggm.gar[i]= sum(A.ggm*A.ggm.tr)/sum(A.ggm.tr)
  F1.ggm.gar[i] = (2*(1-FDR.ggm.gar[i])*power.ggm.gar[i])/(1-FDR.ggm.gar[i] + power.ggm.gar[i])
}

###########################################################################

##########################
# Oracle Estimator for L #
##########################

Z <- matrix(0, p, p)
W <- Z

result.orc.tr=NULL


### Fit the Oracle Estimator
print(paste("Fitting GAR-Oracle with true v0, theta0 and Graph-Zero-pattern "))
result.orc.tr = foreach(i = 1:rep, .maxcombine=max(rep,2))%dopar%{##fit the ith replicate
  
  print(paste("fit replicate ", i))
  S=var(data[[i]])*((n-1)/(n))
  
  SGM::ADMM_L2_Zero(S,theta0.tr, v0.tr,rho=sqrt(log(p)/n), A=net.tr, model=model, Z_ini=Z, W_ini=W, 
                    eps_thre=1e-6, eps_abs=1e-5, eps_rel=1e-3, max_iter=10000, verbose=FALSE)
}

### Extract the Oracle L. Error
L.oracle.err = rep(NA, rep)
for (i in 1:rep){
  
  ## Oracle Estimator of L
  L.orc.i = result.orc.tr[[i]]$L*result.orc.tr[[i]]$theta1
  
  ## Store the Error for the L Oracle Estimator
  L.oracle.err[i] = sum((L.orc.i - theta1.tr*L.tr)^2)/sum((theta1.tr*L.tr)^2)
}



###########################################################################



#########################
# Fit the models (GLASSO)
#########################

### GLASSO Tuning Parameters
C.glasso = exp(seq(log(1.5), log(.15), length.out = 100))
lambda.glasso=C.glasso*sqrt(log(p)/n) 


### Initial Fit to Data using GLASSO
result.glasso=vector("list", rep) # Storage Vector for GLASSO models

## For each sample, fit sequence of GLASSO models
result.glasso=foreach(i = 1:rep, .maxcombine=max(length(rep),2))%dopar%{ ## List of GLASSO models for samples
  
  result.glasso.i = vector("list", length(lambda.glasso))
  
  S=var(data[[i]])*((n-1)/n)
  for (j in 1:length(lambda.glasso)){ 
    print(paste("glasso replicate ", i)) 
    result.glasso.i[[j]] = glasso(s = S, nobs=n, rho = lambda.glasso[j])
  }
  result.glasso.i
}


### Extract Terms involved in eBIC, and Graph Recovery Metrics
A.glasso.net=vector("list", rep) # Storage for Graph Topology
edge.zero=vector("list", rep) # Storage for Zero-pattern
net.glasso.size=vector("list", rep) # Storage for graph size
fdr.glasso = matrix(NA, nrow = rep, ncol = length(lambda.glasso)) # Storage for FDR
power.glasso = matrix(NA, nrow = rep, ncol = length(lambda.glasso)) # Storage for Power
F1.glasso = matrix(NA, nrow = rep, ncol = length(lambda.glasso)) # Storage for F1-Score

###### Calculate terms involved in eBIC-terms
for(i in 1:rep){
  
  ## Initialize i'th element of storage vectors
  A.glasso.net.i = vector("list", length(lambda.glasso))
  net.glasso.size.i = rep(NA, length(lambda.glasso))
  edge.zero.i = vector("list", length(lambda.glasso))
  
  ## For each sample, fit sequence of glasso models
  for (j in 1:length(lambda.glasso)){
    temp=(result.glasso[[i]][[j]]$wi+t(result.glasso[[i]][[j]]$wi))/2 ## make concentration matrix symmetric due to rounding errors
    A.glasso.net.i[[j]]=abs(temp)>1e-6 ## Extract network topology
    net.glasso.size.i[j]=(sum(A.glasso.net.i[[j]])-p)/2 ## Network size
    edge.zero.i[[j]]=which(A.glasso.net.i[[j]]==FALSE, arr.ind = TRUE)
    
    ## Calculate the FDR, F1, Power metrics:
    A.glasso.c = A.glasso.net.i[[j]] # The GGM topology adjacency matrix
    diag(A.glasso.c) = 0 
    
    
    fdr.glasso[i,j] = sum((A.glasso.c*(1-A.ggm.tr)))/sum(A.glasso.c) ## FDR
    power.glasso[i,j] = sum(A.glasso.c*A.ggm.tr)/sum(A.ggm.tr) ## Power
    F1.glasso[i,j] = (2*(1-fdr.glasso[i,j])*power.glasso[i,j])/(1-fdr.glasso[i,j]+ power.glasso[i,j]) ## F1
  }
  
  ## GLASSO network size and graph topology
  A.glasso.net[[i]] = A.glasso.net.i
  net.glasso.size[[i]] = net.glasso.size.i
  edge.zero[[i]] = edge.zero.i
}

### Refit the GLASSO models (re-estimate off-diagonal elements of precision matrix)


## Storage for refitted models
result.glasso.refit=vector("list", rep)
result.glasso.refit=foreach(i = 1:rep, .maxcombine=max(length(rep),2))%dopar%{ 
  
  result.glasso.refit.i = vector("list", length(lambda.glasso))
  S = var(data[[i]])*(n-1)/n
  
  for (j in 1:length(lambda.glasso)){
    
    result.glasso.refit.i[[j]] = glasso(s = S, rho = 1e-10, nobs=n, zero=edge.zero[[i]][[j]], penalize.diagonal=FALSE)
  }
  
  result.glasso.refit.i
}


##################################################
# Model Selection for GLASSO models via eBIC    ##
##################################################

### eBIC penalty
if(p/n>0.5){## e.g., for p=100,n=100
  gamma=1  ## eBIC parameter: p/n~1 set as 1; when p/n<0.5: set as 0.5 
}else{
  gamma=0.5
}


## Store Errors for Covariance matrix
Sigma.glasso.err = matrix(NA, rep, length(lambda.glasso))
## Store Errors for Precision Matrix
Omega.glasso.err = Sigma.glasso.err

## Terms for eBIC 
log.like.glasso=matrix(NA, nrow = rep, ncol = length(lambda.glasso))
bic.glasso=matrix(Inf, nrow = rep, ncol = length(lambda.glasso))
ebic.glasso=matrix(Inf, nrow = rep, ncol = length(lambda.glasso))
P.total=p*(p-1)/2

## Storage for Error Metrics for GLASSO Selected by eBIC
Sigma.sele.glasso = rep(NA, rep)
Omega.sele.glasso = rep(NA, rep)

fdr.sele.glasso = rep(NA, rep)
power.sele.glasso = rep(NA, rep)
F1.sele.glasso = rep(NA, rep)
ebic.sele.glasso = rep(NA, rep)
size.sele.glasso = rep(NA, rep)


### Calculate eBIC for each GLASSO
for (i in 1:rep){ ## For each sample
  S=var(data[[i]])*((n-1)/n) ## Sufficient Statistic for Sample Covariance Matrix
  for (j in 1:length(lambda.glasso)){ ## For each lambda value
    
    ## Estimated Omega (Precision Matrix)
    temp=(result.glasso.refit[[i]][[j]]$wi+t(result.glasso.refit[[i]][[j]]$wi))/2
    
    ## Covariance (Sigma) and Omega (Precision) Error
    Omega.glasso.err[i,j] = sum((temp - Omega.tr)^2)/sum(Omega.tr^2) ## Covariance
    
    Sigma.e = solve(temp)
    Sigma.glasso.err[i,j] = sum((Sigma.e - Sigma.tr)^2)/sum(Sigma.tr^2) ## Precision
    
    ## loglike, BIC, eBIC
    log.like.glasso[i, j]=LogLike.glasso(S=S, Omega=temp, n=n)
    bic.glasso[i, j]=-2*log.like.glasso[i, j]+(net.glasso.size[[i]][j]+p)*log(n)
    ebic.glasso[i, j]=bic.glasso[i,j]+gamma*(lfactorial(P.total)-lfactorial(net.glasso.size[[i]][j])-lfactorial(P.total-net.glasso.size[[i]][j])) ## added 5/19/2024
  }
}


### Conduct Model Selection for eBIC
for (i in 1:rep){
  if(!all(is.infinite(bic.glasso[i, ]))){
    j.sele=which.min(ebic.glasso[i, ]) ## ebic selected glasso
    
    ## Sigma and Omega Errors
    Sigma.sele.glasso[i] = Sigma.glasso.err[i,j.sele]
    Omega.sele.glasso[i] = Omega.glasso.err[i,j.sele]
    
    ## FDR and Power for eBIC selection
    fdr.sele.glasso[i] = fdr.glasso[i,j.sele]
    power.sele.glasso[i] = power.glasso[i,j.sele]
    F1.sele.glasso[i] = F1.glasso[i, j.sele]
    
    size.sele.glasso[i]=net.glasso.size[[i]][j.sele] # Size of models selected by eBIC
    ebic.sele.glasso[i] = ebic.glasso[i, j.sele] # The values of the eBIC
  }
}

################################################################################


##############################
## Results in the Main Text ##
##############################



################################################
## Summary for the GAR Error Metrics (Table 1) #
################################################
table1.metric.names = c("theta0.err", "v0.err", "L.err", "Power", "FDR", "F1")
table1.metric = c(mean(theta0.ebic.err), mean(v0.ebic.err), mean(L.ebic.err), mean(power.ebic.vec),
                  mean(fdr.ebic.vec), mean(F1.ebic.vec))
names(table1.metric) = table1.metric.names

## Print out the Table
print(table1.metric)



#################################################
## Step 1 of GAR (C_l = 1) vs. Oracle Estimator #
#################################################

table2.metric.names = c("theta0.ini", "L.Step1.err", "L.Oracle.err")
table2.metric = matrix(NA, nrow = 1, ncol = 3)
colnames(table2.metric) = table2.metric.names
table2.metric[1,] = c(mean(theta0.s1.err), mean(L.s1.err), mean(L.oracle.err))
print(table2.metric)




################################################
## Summary for the GGM-Induced by GAR Metrics ##
################################################

table3.metric = matrix(NA, nrow = 2, ncol = 5)
colnames(table3.metric) = c("Sigma.err", "Omega.err", "FDR", "Power", "F1")
rownames(table3.metric) = c("GAR", "GLASSO")
table3.metric[1,] = c(mean(Sigma.gar.err), mean(Omega.gar.err), mean(FDR.ggm.gar),
                      mean(power.ggm.gar), mean(F1.ggm.gar))
table3.metric[2,] = c(mean(Sigma.sele.glasso), mean(Omega.sele.glasso), mean(fdr.sele.glasso),
                      mean(power.sele.glasso), mean(F1.sele.glasso))

## Print Results for Table 3: 
print(table3.metric)

