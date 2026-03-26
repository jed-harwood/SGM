# SGM
An R package for fitting Spectral Graph Models

*** 


## Contents
- [Overview](#Overview)
- [Installation](#Installation)
- [Datasets](#Datasets)
- [Main Functions](#Main-Functions)
- [Goodness of Fit](#Goodness-of-Fit)
- [Examples](#Examples)
- [Contact](#Contact)

*** 

## Overview

This repository contains two folders.

* `GAR-Simulation-Scripts`: Contains the R Scripts for replicating the simulation and real-data application in the GAR paper.
* `SGM`: Contains the R package "SGM," for learning undirected graphs through a Graphical Auto-Regressive (GAR) model through a 3-step fitting procedure.


## Installation 

### Dependencies
Please make sure that the following packages are installed before using the R package `SGM`. 

```
install.packages(c("doParallel", "foreach", "gmp"))
```

### Installation from the Repository
The R Package `SGM` can be installed from source code in this GitHub repository.  The `devtools` package is required for installation.  
```
library(devtools)
install_github(repo="jed-harwood/SGM", subdir="SGM")
```

## Datasets

This repository currently contains two datasets.  To load a dataset into the working environment, run `data("<dataname>")`.  


* `gar1`: A sample simulated from an underlying GAR(1) model, with p=100, n=100.  The underrlying graph was a random graph with edge probability of 0.02.
* `stocks`: A collection of (standardized) log-returns from 283 stocks on the S\&P 500.  The dataset spans January 1, 2007, to January 1, 2011, covering the global financial crisis, with 1007 closing prices per stock.  The stocks come from five GICS sectors: 58 from Information Technology, 72 from Consumer Discretionary, 32 from Consumer Staples, 59 from Financials, and 62 from Industrials.

These datasets can be accessed as the following R objects.

`gar1`:  A list object that contains:
1. `data`: data simulated under a `GAR(1)` model.
2. `A.tr`: the true weighted adjacency matrix for the underlying graph.
3. `LN`: the true normalized graph Laplacian matrix.
4. `theta0`, `theta1`: the true graph filter parameters.  


`stocks`: A 1007 by 283 matrix.


For more information on a given dataset, please run `?<dataname>`.  

***

## Main-Functions

`GAR1_fit`: fit a `GAR(1)` model for a given set of tuning parameters, using a 3-step estimation procedure based on the penalized MLE.  
* Step 0: given an initial estimate `S` for the covariance matrix (e.g., the sample covariance matrix), obtain an initial estimate for `theta0` by the reciprocal of the largest eigenvalue of `S`, square-rooted.
* Step 1: given `S` and `theta0` from Step 0, use an ADMM algorithm to estimate the (normalized) graph Laplacian `L` of the latent graph and extract the thresholded zero-pattern used in later steps.
* Step 2: given the zero-pattern from Step 1, use an ADMM algorithm to refit the non-zero elements of the (normalized) graph Laplacian `L`.
* Step 3:
    * a. given the (normalized) graph Laplacian `L` from Step 2, use an ADMM algorithm to estimate the normalized degree vector `v0`.
    * b. given `v0` and the zero-pattern from Step 1, use an ADMM algorithm to simultaneously re-estimate `theta0` and the (normalized) graph Laplacian `L`.

`model_selec`: conduct model selection via the eBIC criterion.

### Arguments

**Table: Arguments for `GAR1_fit`**

| Parameter | Data type | Default | Description |
|-----------|-----------|---------|-------------|
| S | matrix |  | Estimate of the covariance matrix (e.g., the MLE). |
| nobs | integer |  | Number of samples used to compute `S`. |
| lambda.v | numeric vector |  | Tuning parameter controlling sparsity of the estimated graph. |
| net.thre | numeric vector |  | Thresholding parameter for removing noisy edges; used in Step 2 and beyond of the GAR fitting procedure. |
| model | character | "LN" | Type of graph Laplacian: "LN" (normalized Laplacian), "L" (graph Laplacian), or "LN.noselfloop" (normalized Laplacian without self-loops). |
| step | integer | 3 | Number of steps in the estimation procedure (1, 2, or 3). |
| rho.v | numeric vector | lambda.v | ADMM penalty parameter (typically set equal to `lambda.v`). |
| eps_thre | numeric scalar | 1e-6 | Small positive threshold used for numerical stability. |
| eps_abs | numeric scalar | 1e-5 | Absolute tolerance for ADMM convergence. |
| eps_rel | numeric scalar | 1e-3 | Relative tolerance for ADMM convergence. |
| max_iter_s1 | integer | 10000 | Maximum number of iterations for Step 1. |
| max_iter_s2 | integer | 10000 | Maximum number of iterations for Step 2. |
| max_iter_s3 | integer | 10000 | Maximum number of iterations for Step 3b, the joint estimation sub-step within Step 3. |
| verbose | logical | FALSE | Logical flag indicating whether to print progress during fitting. |


**Table: Arguments for `model_selec`**

| Parameter | Data type | Default | Description |
|-----------|-----------|---------|-------------|
| resultList | list |  | A list output from `GAR1_fit`. The required metadata, including `nobs`, `step`, and `model`, are read directly from this object. |




### Value


**Table: Output of `GAR1_fit`**

| Output | Data type | Description |
|--------|-----------|-------------|
| S | matrix | The \( p \times p \) covariance matrix supplied to `GAR1_fit`. |
| nobs | integer | Number of observations used to compute `S`. |
| model | character | The fitted GAR model family. |
| step | integer | The last step requested in the fitting procedure. |
| lambda.v | numeric vector | Sparsity tuning parameter sequence. |
| rho.v | numeric vector | ADMM tuning parameter sequence. |
| net.thre | numeric vector | Threshold sequence used in Steps 2 and 3. |
| theta0.init | numeric scalar | Initial estimate of \( \theta_0 \) from Step 0, used in Steps 1 and 2. |
| theta0.s3 | numeric matrix | Matrix of Step 3b estimates for \( \theta_0 \) (NULL if `step < 3`). |
| A.net | list | Estimated graph topologies indexed by `lambda.v` and `net.thre`, created in Step 1 and used in Steps 2 and 3. |
| step1 | list | Step 1 ADMM fit objects indexed by `lambda.v`. Each element is the returned output of `ADMM_L2`, including quantities such as `L`, `Z`, `W`, `theta0`, `theta1`, and `conv`, and Step 1 also defines the thresholded zero-patterns stored in `A.net`. |
| step2 | list | Step 2 refit objects indexed by `lambda.v` and `net.thre` (NULL if `step < 2`). Each element is the returned output of `ADMM_L2_Zero`, containing the refitted `L`, `Z`, `W`, `theta0`, `theta1`, and `conv` for the corresponding zero-pattern from Step 1. |
| step3a | list | Step 3a `v0` estimation results indexed by `lambda.v` and `net.thre` (NULL if `step < 3`). Each element stores the estimated `v0` vector produced from the Step 2 Laplacian path. |
| step3b | list | Step 3b joint estimation objects indexed by `lambda.v` and `net.thre` (NULL if `step < 3`). Each element is the returned output of `ADMM_Lap_Zero`, containing the jointly refitted `L`, `theta0`, `theta1`, `Z`, `W`, `phi`, and `conv`. |
| v0.s3 | list | Step 3a `v0` estimates indexed by `lambda.v` and `net.thre` (NULL if `step < 3`). |
| conv | list | Convergence diagnostics stored as `conv$step1`, `conv$step2`, `conv$step3a`, and `conv$step3b`. |



**Table: Output of `model_selec`**

| Output | Data type | Description |
|--------|-----------|-------------|
| selected.model | list | The ADMM output for the eBIC-selected model. |
| theta0 | numeric scalar | Selected estimate of \( \theta_0 \). |
| theta1 | numeric scalar | Selected estimate of \( \theta_1 \). |
| L | matrix | Selected graph Laplacian estimate. |
| v0 | numeric vector | Selected `v0` estimate (NULL when `step = 2`). |
| A.net.e | matrix | The (unweighted) graph topology selected by the eBIC criterion. |
| index | matrix | Index corresponding to the optimal tuning parameters (`lambda`, `net.thre`). |
| lambda.v | numeric scalar | Selected value of `lambda.v`. |
| net.thre | numeric scalar | Selected value of `net.thre`. |
| ebic | numeric scalar | eBIC score of the selected model. |




***

## Examples

For more information on the `GAR1_fit` and `model_selec` functions, run `?GAR1_gf`, `?GAR1_fit` and `?model_selec` in R.  

### Fit `GAR(1)` to the simulated `gar1` data. 
```
### load libraries
library(SGM)

### See ?gar1 for more information
data("gar1")
str(gar1)


### Get data 
gar1_data = gar1$data
nobs = nrow(gar1_data)
p = ncol(gar1_data)

### Set model to fit
model = "LN"

### Set tuning parameters: lambda and net.thre sequence
C.v=c(1,0.5)  
lambda.v=C.v*sqrt(log(p)/nobs)

C.thre=exp(seq(log(1),log(0.05), length.out=10))
net.thre=C.thre*sqrt(log(p)/nobs)

### Set ADMM parameter 
rho.v=pmax(lambda.v, 0.01)

### Get sample covariance 
S = var(gar1_data)*(nobs-1)/nobs

### Fit GAR(1) (up to step 3)
fit = GAR1_fit(S, nobs, lambda.v, net.thre, model, 3, rho.v)

### Model selection via eBIC for the GAR estimator
fit.gar1.3 = model_selec(fit)

### Evaluation: estimation errors for theta0, theta1*L, and FDR and Power for graph inference 
## Get ground truth
A.tr = gar1$A.tr > 0 # True 0-1 adjacency matrix
diag(A.tr)=0
LN = gar1$LN # True (normalized) graph Laplacian
theta0.tr = gar1$theta0 # True graph filter parameter
theta1.tr = gar1$theta1 # True graph filter parameter

## Calculate estimation errors
theta0.err.3 = abs(fit.gar1.3$theta0 - theta0.tr)^2
L.err.3 = sum((fit.gar1.3$theta1 * fit.gar1.3$L - theta1.tr*LN)^2)/sum((theta1.tr*LN)^2)

## Calculate FDR and Power 
FDR.3 = sum(fit.gar1.3$A.net.e*(1-A.tr))/sum(fit.gar1.3$A.net.e)
Power.3 = sum(fit.gar1.3$A.net.e*A.tr)/sum(A.tr)

## Network sizes
net.size.3 = sum(fit.gar1.3$A.net.e)/2

## View results 
c(theta0.err.3) # theta0 errors
c(L.err.3) # L errors
c(net.size.3) # network sizes
c(Power.3) # Power
c(FDR.3) # FDR

## Results:
## p=100; n=100; true network size = 105; 
## theta0: 0.00967012
## L: 0.01916969
## Sizes: 103
## Power: 0.9238095
## FDR: 0.05825243
```



***

## Goodness-of-Fit

The packages offers a goodness of fit measure as proposed in Harwood, Paul, and Peng (2024).  

`GAR1_gf`: calculate a goodness-of-fit measure to determine if `GAR(1)` is an appropriate model for the data. Valid for $n \geq p$.  Return a value between $0$ and $1$. If the returned value is close to $1$, then it means that the `GAR(1)` model is a good fit to the data.


**Table: Arguments for `GAR1_gf`**

| Parameter    | Default   | Description |
|--------------|----------|-------------|
| S            |          | Estimate of the covariance matrix (e.g., the MLE). |
| nobs         |          | Number of observations used to compute `S`. |
| lambda.v     |          | Positive tuning parameter for the GAR(1) model. This should typically be set to `sqrt(log(p)/n)` where `p` is the dimension, and `n` is the sample size. |
| rho.v        | lambda.v | ADMM penalty parameter (positive; typically set equal to `lambda.v`). |
| eps_thre     | 1e-6     | Small positive threshold used for numerical stability. |
| eps_abs      | 1e-5     | Absolute tolerance for ADMM convergence. |
| eps_rel      | 1e-3     | Relative tolerance for ADMM convergence. |
| max_iter     | 10000    | Maximum number of iterations for fitting on observed data. |
| num.thread   | 1        | Number of threads used for computation. |
| rep.boot     | 100      | Number of bootstrap samples used for the goodness-of-fit test. |
| seed         | 1        | Random seed for reproducibility. |


**Table: Output of `GAR1_gf`**

| Output  | Description |
|---------|-------------|
| p-value | P-value for the goodness-of-fit test. |







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

# > 1

```

***

## Contact
Please report any bugs to `jedharwood@ucdavis.edu`
