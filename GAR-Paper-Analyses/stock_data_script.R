######################################################
## Apply GAR Model to the S&P 500 Data              ##
## Requires: stock_data_processing.R, stockdata.rda ##
######################################################

## Load SGM and the Stock data
rm(list=ls())
library(SGM)

## Load pre-processed data
## for preprocessing steps, "/stock-auxiliary-scripts/stock_data_processing.R"
load("stock.data.X.Rdata")
n=nrow(X)
p=ncol(X)

#######################
## Fit the GAR Model
########
C.v=c(8,4,1)
lambda.v=C.v*sqrt(log(p)/n) 

C.thre=exp(seq(log(1),log(0.1), length.out=12)) # decreasing seq for epsilon.thre vals
net.thre=C.thre*sqrt(log(p)/n) 

rho.v=pmax(lambda.v, 0.01) # ADMM parameters
S = var(X) * (n-1)/n # Sample covariance matrix

############################
## Goodness of Fit Measure
##########
GAR1_gf(S = S, nobs = n, lambda.v = sqrt(log(p)/n), num.thread = 10, seed = 1)



#####################
## Fit the GAR model
######
gar.stock = GAR1_fit(S = S, nobs = n, lambda.v = lambda.v, net.thre = net.thre, model = "LN",
                     step = 3, rho.v = rho.v, verbose = FALSE, max_iter_s1 = 100000, max_iter_s2 = 100000, max_iter_s3 = 100000)

save(gar.stock, file = "gar_stock.rda")

###################
## Fit the GLASSO 
###############
library(glasso)
C.glasso=exp(seq(log(11), log(6), length.out=100)) ## decreasing seq.    
lambda.glasso=C.glasso*sqrt(log(p)/n) 

## Storage lists for GLASSO (with refitting)
result.glasso=vector("list", length(lambda.glasso))
result.glasso.refit = result.glasso

## glasso:
A.glasso.net=vector("list", length(lambda.glasso)) # storage for graph topology
for (j in 1:length(lambda.glasso)){
  ## Initial Fit
  result.glasso[[j]]=glasso(s = S, nobs=n, rho = lambda.glasso[j])
  
  ## Extract graph topology and refit
  ### Topology
  A.glasso.net[[j]]=abs(result.glasso[[j]]$wi)>1e-6
  edge.zero=which(A.glasso.net[[j]]==FALSE, arr.ind = TRUE)
  
  ## Refit
  result.glasso.refit[[j]]=glasso(s = S, rho = 1e-10, nobs=n, zero=edge.zero, penalize.diagonal=FALSE)
}
save(result.glasso.refit, file = "glasso_refit.rda")






############
## Model Evaluation
#########
load("gar_stock.rda")
load("glasso_refit.rda")

source("stock-auxiliary-scripts/stock_gar_res_process.R") # Script for Processing Results; Gives function stock_results()

## Get results of Model Fits
### (GAR)
GAR_res = stock_results(modelList = gar.stock, model = "GAR", p = p, n = n, sp.num = sp.num, S = S)


### (GLASSO)
GLASSO_res = stock_results(result.glasso.refit, model = "glasso", p = p, n = n, sp.num = sp.num, S = S)



#### Print out table of stock connectivity
print(GAR_res$table)




#########################################
## Plot the number of parameters vs. log-likelihood
##############################

## GAR likelihood and parameter sizes
GAR.net.size = GAR_res$size
GAR.log.like = GAR_res$loglike
GAR.selec = GAR_res$index.selec
GAR.size.selec = GAR_res$size.selec

## GLASSO likelihood and parameter sizes
glasso.net.size = GLASSO_res$size
glasso.log.like = GLASSO_res$loglike
glasso.selec = GLASSO_res$index.selec
glasso.size.selec = GLASSO_res$size.selec



## Begin Plot comapring GLASSO to GAR
pdf("stock-GAR1LN-loglike-GAR-vs-glasso.pdf")
plot(GAR.net.size+1+p, GAR.log.like, xlab='# of parameters in model', ylab="log-likelihood", 
     xlim=range(GAR.net.size+1+p, glasso.net.size+p), 
     ylim=range(GAR.log.like,glasso.log.like), col=1, type='n')

## Add the points for GAR
points(GAR.net.size+1+p, GAR.log.like,col=2, pch=2, type='p') ## GAR-LN model 
points(x=GAR.size.selec+1+p, y=GAR.log.like[GAR.selec[1], GAR.selec[2]], pch=17, col=2, type='p')
abline(v=GAR.size.selec+1+p, col=2, lty=2)

## Add the points for GLASSO
points(glasso.net.size+p, glasso.log.like, col=1, pch=1, type='p') ## Glasso model 
points(x=glasso.size.selec+p, y=glasso.log.like[glasso.selec],pch=16, col=1)
abline(v=glasso.size.selec+p, col=1, lty=2)

legend("bottomright",  c("GAR","Glasso"), pch = c(2,1), col = c(2,1), cex = 0.5) 

dev.off()
