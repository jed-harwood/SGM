# SGM
An R package for fitting Spectral Graph Models

*** 
The R package `SGM` provides a set of functions that learn a latent graph from the data by modeling the observations as stationary signals on this graph, i.e., the covariance matrix co-diagonalizes with the (normalized) graph Laplacian.  
Currently, the package offers functions to fit a `graphical autoregressive -- GAR(1)` model through a 3-step procedure, as well as a `time-varying autoregressive -- TARGAR` model.  

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

`gar1`:  A list object that contains:
1. `data`: data simulated under a `GAR(1)` model.
2. `A.tr`: the true weighted adjacency matrix for the underlying graph.
3. `LN`: the true normalized graph Laplacian matrix.
4. `theta0`, `theta1`: the true graph filter parameters.  

*** 

## Examples

For more information on the `GAR1_gf`, `GAR1_fit` and `model_selec` functions, run `?GAR1_gf`, `?GAR1_fit` and `?model_selec` in R.  

### Fit `GAR(1)` to the simulated `gar1` data. 
```
### load libraries
library(SGM)

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

### Goodness-of-fit measure
GAR1_gf(S, n, lambda.v[1], num.thread = 2)

### Fit GAR(1) (up to step 3)
fit = GAR1_fit(S, n, lambda.v, net.thre, model, 3, rho.v)

### Model selection via eBIC for Step 2 estimator and Step 3 estimator
fit.gar1.2 = model_selec(fit, n, step = 2, model)
gar1_fitted.2 = fit.gar1.2$model.selec

fit.gar1.3 = model_selec(fit, n, step = 3, model)
gar1_fitted.3 = fit.gar1.3$model.selec

### Evaluation: estimation errors for theta0, theta1*L, and FDR and Power for graph inference 
## Get ground truth
A.tr = gar1$A.tr > 0 # True 0-1 adjacency matrix
diag(A.tr)=0
LN = gar1$LN # True (normalized) graph Laplacian
theta0.tr = gar1$theta0 # True graph filter parameter
theta1.tr = gar1$theta1 # True graph filter parameter

## Calculate estimation errors
theta0.err.2 = abs(gar1_fitted.2$theta0 - theta0.tr)^2
theta0.err.3 = abs(gar1_fitted.3$theta0 - theta0.tr)^2

L.err.2 = sum((gar1_fitted.2$theta1 * gar1_fitted.2$L - theta1.tr*LN)^2)/sum((theta1.tr*LN)^2)
L.err.3 = sum((gar1_fitted.3$theta1 * gar1_fitted.3$L - theta1.tr*LN)^2)/sum((theta1.tr*LN)^2)

## Calculate FDR and Power 
FDR.2 = sum(fit.gar1.2$A.net.e*(1-A.tr))/sum(fit.gar1.2$A.net.e)
FDR.3 = sum(fit.gar1.3$A.net.e*(1-A.tr))/sum(fit.gar1.3$A.net.e)

Power.2 = sum(fit.gar1.2$A.net.e*A.tr)/sum(A.tr)
Power.3 = sum(fit.gar1.3$A.net.e*A.tr)/sum(A.tr)

## Network sizes
net.size.2 = sum(fit.gar1.2$A.net.e)/2
net.size.3 = sum(fit.gar1.3$A.net.e)/2

## View results 
c(theta0.err.2, theta0.err.3) # theta0 errors
c(L.err.2, L.err.3) # L errors
c(net.size.2, net.size.3) # network sizes
c(Power.2, Power.3) # Power
c(FDR.2, FDR.3) # FDR

## Results:
## p=100; n=100; true network size = 105; 
## theta0: 0.06893404 0.00967012
## L: 0.03392529 0.01916969
## Sizes: 101 103
## Power: 0.9047619 0.9238095
## FDR: 0.05940594 0.05825243
```

### Fit `TARGAR` to the simulated `targar` data
```
### load libraries
library(SGM)

### See ?gar1 for more information
data("targar")

### Get data 
targar1_dat = targar$data
t = nrow(targar1_dat) # Sample Size
p = ncol(targar1_dat) # Dimensionality

### Set tuning parameters: lambda and net.thre sequence
C.v=c(1,0.5)  
lambda.v=C.v*sqrt(log(p)/(t-1))

C.thre=exp(seq(log(1),log(0.05), length.out=10))
net.thre=C.thre*sqrt(log(p)/(t-1))

### Set ADMM parameter 
rho.v=pmax(lambda.v, 0.01)

### Fit TARGAR with 3-passes
fit = TARGAR_fit(targar1_dat, lambda.v, net.thre, rho.v, num.pass = 3)


### Model selection via eBIC for TARGAR Estimates
fit.targar = model_selec(fit, t, model = "TARGAR")
targar1_fitted = fit.targar$model.selec
targar1_est = targar1_fitted$result.0S

### Evaluation: estimation errors for theta0, theta1*L, and FDR and Power for graph inference 
## Get ground truth
A.tr = targar$A.tr > 0 # True 0-1 adjacency matrix
diag(A.tr)=0
LN = targar$LN # True (normalized) graph Laplacian
theta0.tr = targar$theta0 # True graph filter parameter
theta1.tr = targar$theta1 # True graph filter parameter

## Calculate estimation errors
theta0.err = abs(targar1_est$theta0 - theta0.tr)^2
L.err = sum((targar1_est$theta1 * targar1_est$L - theta1.tr*LN)^2)/sum((theta1.tr*LN)^2)

## Calculate FDR and Power 
FDR = sum(targar1_fitted$A.net*(1-A.tr))/sum(targar1_fitted$A.net)

Power = sum(targar1_fitted$A.net*A.tr)/sum(A.tr)

## Network sizes
net.size.est = sum(targar1_fitted$A.net)/2

## View results 
(theta0.err) # theta0 errors
(L.err) # L errors
(net.size.est) # network sizes
(Power) # Power
(FDR) # FDR

## Results:
## p=100; n=100; true network size = 105; 
## theta0: 0.01596504
## L: 0.02927149
## Size: 101
## Power: 0.9428571
## FDR: 0.01980198
```


***

## Contact
Please report any bugs to `jedharwood@ucdavis.edu`
