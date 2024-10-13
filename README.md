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

`GAR1_gf`: run a goodness of fit test to determine if GAR(1) is appropriate. Valid for n >= p.

`GAR1_fit`: learn a GAR(1) SGM for a given set of tuning parameters

`model_selec`: conduct model selection via the eBIC criterion

***

## Data
To load a dataset into the working environment, run `data("<dataname>")`.  For more information on a given dataset, please run `?<dataname>`.  

`stocks`: Data on the S&P 500 stock prices

`gar1`:  Simulated data that comes from an underlying GAR(1) model, generated to have a (latent) graph with edge probability `0.02`, graph filter parameters `theta0=1`, and `theta1=2`, with no self-loops nor isolated vertices. 

## Contact
Please report any bugs to `jedharwood@ucdavis.edu`

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
data("gar1")
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

### Model selection via eBIC
optGAR1 = model_selec(fit, n, step = 2, model)

### Summarize Errors for true theta0, theta1


```
