# SGM
An R package for Estimating Spectral Graphical Models

*** 
The R package `SGM` provides a set of functions that learn a latent graph from the data.  The data are viewed as stationary signals from the inferred latent graph.  `SGM` is motivated by the need to understand complex dependence structures, and how they change over time.  Currently, the package offers code to fit a GAR(1) model (as proposed by Peng et al.) through a 3-step procedure.

*** 

## Dependencies 
Please make sure that the following packages are installed before using the R package `SGM`. 

```
install.packages(c("doParallel", "foreach", "gmp"))
```

## Installation
The R Package `SGM` can be installed from source files, in this GitHub repository.  The `devtools` package is required for installation.  
```
library(devtools)
install_github(repo="jed-harwood/SGM")
```

## Main Functions

`GAR1_gf`: run a goodness of fit test to determine if GAR(1) is appropriate. Valid for n >= p.  Return a value between 0 and 1. If the returned value is close to 1, then it means that the GAR(1) model is a good fit to the data.

`GAR1_fit`: learn a GAR(1) SGM for a given set of tuning parameters, using a 3 step estimation procedure based on a penalized MLE.  
* Step 1: given an initial estimate (S) for the covariance matrix (e.g., by the sample covariance matrix) and an initial estimate for theta0, uses an ADMM algorithm to estimate the (normalized) graph Laplacian.  The initial estimate for theta0 is the reciprocal of the largest eigenvalue of S, squared.
* Step 2: given the zero-pattern in the Laplacian from Step 1, uses an ADMM algorithm to refit the off-diagonal elements of the (normalized) graph Laplacian.
* Step 3:
    a. given the (normalized) graph Laplacian from Step 2, uses an ADMM algorithm to estimate the degree vector v0.
    b. given v0 and the zero-pattern from Step 1, uses an ADMM algorithm to simultaneously re-estimate theta0 and the (normalized) graph Laplacian.

`model_selec`: conduct model selection via the eBIC criterion

***

## Data
To load a dataset into the working environment, run `data("<dataname>")`.  For more information on a given dataset, please run `?<dataname>`.  

`stocks`: Data on the S&P 500 stock prices

`gar1`:  A list object that contains:
1. `data`: Simulated data that comes from an underlying GAR(1) model, generated to have a (latent) graph with edge probability `0.02`, graph filter parameters `theta0=1`, and `theta1=2`, with no self-loops nor isolated vertices.
2. `A.tr`: The true adjacency matrix.
3. `LN`: The true (normalized) graph Laplacian

*** 

## Examples

For more information on the `GAR1_gf`, `GAR1_fit` and `model_selec` functions, run `?GAR1_gf`, `?GAR1_fit` and `?model_selec` in R.  

### Using `stocks`
```
data("stocks")
n=nrow(stocks)
p=ncol(stocks)
model="LN"

### lambda and net.thre sequence
C.v=c(8,4,1)
lambda.v=C.v*sqrt(log(p)/n) 
rho.v=pmax(lambda.v, 0.01)

C.thre=exp(seq(log(1),log(0.1), length.out=12))
net.thre=C.thre*sqrt(log(p)/n)

### Determine if GAR(1) model is appropriate
S = var(stocks)*(n-1)/n
GAR1_gf(S, n, lambda.v[1])

### Run GAR1_fit (up to step 3)
resList = GAR1_fit(S, n, lambda.v, net.thre, model)

### Conduct model selection via eBIC
optModel = model_sele(resList, n, step = 3, model = "LN")
```


### Using `gar1`
```
### ground truth (see ?gar1 for more information)
data("gar1")
A.tr = gar1$A.tr # True adjacency matrix
LN = gar1$LN # True (normalized) graph Laplacian
gar1 = gar1$data # Simulated data
theta0.tr = 1 
theta1.tr = 2

### extract p and n
n = nrow(gar1)
p = ncol(gar1)
model = "LN"

### lambda and net.thre sequence
C.v=c(1,0.5)  
lambda.v=C.v*sqrt(log(p)/n)

C.thre=exp(seq(log(1),log(0.05), length.out=10))
net.thre=C.thre*sqrt(log(p)/n)

### ADMM parameter 
rho.v=pmax(lambda.v, 0.01)

### Fit GAR(1) (up to step 3)
S = var(gar1)*(n-1)/n
fit = GAR1_fit(S, n, lambda.v, net.thre, model, 3, rho.v)

### Model selection via eBIC (Step 2 and Step 3)
fit.gar1.2 = model_selec(fit, n, step = 2, model)
fit.gar1.3 = model_selec(fit, n, step = 3, model)

### Evaluate Estimation Errors for theta0, L, and FDR, Power

theta0.err.2 = abs(fit.gar1.2$model.selec$theta0 - theta0.tr)^2
theta0.err.3 = abs(fit.gar1.3$model.selec$theta0 - theta0.tr)^2

L.err.2 = sum((fit.gar1.2$model.selec$L*fit.gar1.2$model.selec$theta1 - theta1.tr*LN)^2)/sum(theta1.tr*LN^2)
L.err.3 = sum((fit.gar1.3$model.selec$L*fit.gar1.3$model.selec$theta1 - theta1.tr*LN)^2)/sum(theta1.tr*LN^2)

FDR.2 = sum(fit.gar1.2$A.net.e*(1-A.tr))/sum(fit.gar1.2$A.net.e)
FDR.3 = sum(fit.gar1.3$A.net.e*(1-A.tr))/sum(fit.gar1.3$A.net.e)

Power.2 = sum(fit.gar1.2$A.net.e*A.tr)/sum(A.tr)
Power.3 = sum(fit.gar1.3$A.net.e*A.tr)/sum(A.tr)


### View Errors
c(theta0.err.2, theta0.err.3) # theta0 errors
c(L.err.2, L.err.3) # L errors
c(Power.2, Power.3) # Power
c(FDR.2, FDR.3) # FDR 

```

***

## Contact
Please report any bugs to `jedharwood@ucdavis.edu`
