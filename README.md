# SGM
An R package for fitting Spectral Graph Models

*** 
The R package `SGM` provides a set of functions that learn a latent graph from the data by modeling the observations as stationary signals on this graph, i.e., the covariance matrix co-diagonalizes with the (normalized) graph Laplacian.  
Currently, the package offers functions to fit a `graphical autoregressive -- GAR(1)` model through a 3-step procedure.

*** 

## Dependencies 
Please make sure that the following packages are installed before using the R package `SGM`. 

```
install.packages(c("doParallel", "foreach", "gmp"))
```

## Installation
The R Package `SGM` can be installed from source code in this GitHub repository.  The `devtools` package is required for installation.  
```
library(devtools)
install_github(repo="jed-harwood/SGM")
```

## Main Functions

`GAR1_gf`: calculate a goodness-of-fit measure to determine if `GAR(1)` is an appropriate model for the data. Valid for $n \geq p$.  Return a value between $0$ and $1$. If the returned value is close to $1$, then it means that the `GAR(1)` model is a good fit to the data.

`GAR1_fit`: fit a `GAR(1)` model for a given set of tuning parameters, using a 3-step estimation procedure based on the penalized MLE.  
* Step 0: given an initial estimate `S` for the covariance matrix (e.g., the sample covariance matrix), obtain an initial estimate for `theta0` by the reciprocal of the largest eigenvalue of `S`, square-rooted.
* Step 1: given `S`  and `theta0` from Step 0, use an ADMM algorithm to estimate the (normalized) graph Laplacian `L` of the latent graph.  
* Step 2: given the zero-pattern in `L` from Step 1, use an ADMM algorithm to refit the non-zero elements of the (normalized) graph Laplacian `L`.
* Step 3:
    * a. given the (normalized) graph Laplacian `L` from Step 2, use an ADMM algorithm to estimate the normalized degree vector `v0`.
    * b. given `v0` and the zero-pattern of `L` from Step 1, use an ADMM algorithm to simultaneously re-estimate `theta0` and the (normalized) graph Laplacian `L`.

`model_selec`: conduct model selection via the eBIC criterion.

***

## Data
To load a dataset into the working environment, run `data("<dataname>")`.  For more information on a given dataset, please run `?<dataname>`.  

`stocks`: (standardized) log-returns of $283$ S&P 500 stocks from January 1, 2007 to January 1, 2011, based on raw data obtained via the R package `quantmod`. 

`gar1`:  A list object that contains:
1. `data`: data simulated under a `GAR(1)` model.
2. `A.tr`: the true adjacency matrix for the underlying graph.
3. `LN`: the true normalized graph Laplacian matrix.
4. `theta0`, `theta1`: the true graph filter parameters.  

*** 

## Examples

For more information on the `GAR1_gf`, `GAR1_fit` and `model_selec` functions, run `?GAR1_gf`, `?GAR1_fit` and `?model_selec` in R.  

### Fit `GAR(1)` to `stocks` data.
```
### See ?stocks for more information
data("stocks")
n=nrow(stocks)
p=ncol(stocks)

### Set the model to fit: GAR(1)-normalized Laplacian model
model="LN" 

### Set tuning parameters: lambda and net.thre sequence
C.v=c(8,4,1)
lambda.v=C.v*sqrt(log(p)/n)

C.thre=exp(seq(log(1),log(0.1), length.out=12))
net.thre=C.thre*sqrt(log(p)/n)

### Set ADMM parameter
rho.v=pmax(lambda.v, 0.01)

### Calculate sample covariance
S = var(stocks)*(n-1)/n

### Determine if GAR(1) model is appropriate
GAR1_gf(S, n, lambda.v[1])

### Run GAR1_fit (up to step 3)
resList = GAR1_fit(S, n, lambda.v, net.thre, model, 3, rho.v)

### Conduct model selection via eBIC
optModel = model_sele(resList, n, step = 3, model = "LN")
```


### Fit `GAR(1)` to the simulated `gar1` data. 
```
### See ?gar1 for more information
data("gar1")

### Get data 
gar1_data = gar1$data
n = nrow(gar1_data)
p = ncol(gar1_data)

### Set model to fit
model = "LN"

### Set tuning parameters: lambda and net.thre sequence
C.v=c(1,0.5)  
lambda.v=C.v*sqrt(log(p)/n)

C.thre=exp(seq(log(1),log(0.05), length.out=10))
net.thre=C.thre*sqrt(log(p)/n)

### Set ADMM parameter 
rho.v=pmax(lambda.v, 0.01)

### Get sample covariance 
S = var(gar1_data)*(n-1)/n

### Fit GAR(1) (up to step 3)
fit = GAR1_fit(S, n, lambda.v, net.thre, model, 3, rho.v)

### Model selection via eBIC for Step 2 estimator and Step 3 estimator
fit.gar1.2 = model_selec(fit, n, step = 2, model)
gar1_fitted.2 = fit.gar1.2$model.selec

fit.gar1.3 = model_selec(fit, n, step = 3, model)
gar1_fitted.3 = fit.gar1.3$model.selec

### Evaluation: estimation errors for theta0, theta1*L, and FDR and Power for graph inference 
## Get ground truth
A.tr = gar1$A.tr # True adjacency matrix
LN = gar1$LN # True (normalized) graph Laplacian
theta0.tr = gar1$theta0 # True graph filter parameter
theta1.tr = gar1$theta1 # True graph filter parameter

## Calculate estimation errors
theta0.err.2 = abs(gar1_fitted.2$theta0 - theta0.tr)^2
theta0.err.3 = abs(gar1_fitted.3$theta0 - theta0.tr)^2

L.err.2 = sum((gar1_fitted.2$theta1 * gar1_fitted.2$L - theta1.tr*LN)^2)/sum(theta1.tr*LN^2)
L.err.3 = sum((gar1_fitted.3$theta1 * gar1_fitted.3$L - theta1.tr*LN)^2)/sum(theta1.tr*LN^2)

## Calculate FDR and Power 
FDR.2 = sum(fit.gar1.2$A.net.e*(1-A.tr))/sum(fit.gar1.2$A.net.e)
FDR.3 = sum(fit.gar1.3$A.net.e*(1-A.tr))/sum(fit.gar1.3$A.net.e)

Power.2 = sum(fit.gar1.2$A.net.e*A.tr)/sum(A.tr)
Power.3 = sum(fit.gar1.3$A.net.e*A.tr)/sum(A.tr)

## View results 
c(theta0.err.2, theta0.err.3) # theta0 errors
c(L.err.2, L.err.3) # L errors
c(Power.2, Power.3) # Power
c(FDR.2, FDR.3) # FDR

## Results:
## theta0: 0.06893404 0.00967012
## L: 0.06785057 0.03833937
## Power: 0.9187622 0.9350193
## FDR: 0.2930108 0.2944717 
```

***

## Contact
Please report any bugs to `jedharwood@ucdavis.edu`
