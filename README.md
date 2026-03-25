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
* Step 1: given `S`  and `theta0` from Step 0, use an ADMM algorithm to estimate the (normalized) graph Laplacian `L` of the latent graph.  
* Step 2: given the zero-pattern in `L` from Step 1, use an ADMM algorithm to refit the non-zero elements of the (normalized) graph Laplacian `L`.
* Step 3:
    * a. given the (normalized) graph Laplacian `L` from Step 2, use an ADMM algorithm to estimate the normalized degree vector `v0`.
    * b. given `v0` and the zero-pattern of `L` from Step 1, use an ADMM algorithm to simultaneously re-estimate `theta0` and the (normalized) graph Laplacian `L`.

`model_selec`: conduct model selection via the eBIC criterion.

### Arguments

**Table: Arguments for `GAR1_fit`**

| Parameter    | Default   | Description |
|--------------|----------|-------------|
| S            |          | Estimate of the covariance matrix (e.g., the MLE). |
| nobs         |          | Number of samples used to compute `S`. |
| lambda.v     |          | Tuning parameter controlling sparsity of the estimated graph. |
| net.thre     |          | Thresholding parameter for removing noisy edges; used in Step 2 and beyond of the GAR fitting procedure. |
| model        | "LN"     | Type of graph Laplacian: "LN" (normalized Laplacian), "L" (graph Laplacian), or "LN.noselfloop" (normalized Laplacian without self-loops). |
| step         | 3        | Number of steps in the estimation procedure (1, 2, or 3). |
| rho.v        | lambda.v | ADMM penalty parameter (typically set equal to `lambda.v`). |
| eps_thre     | 1e-6     | Small positive threshold used for numerical stability. |
| eps_abs      | 1e-5     | Absolute tolerance for ADMM convergence. |
| eps_rel      | 1e-3     | Relative tolerance for ADMM convergence. |
| max_iter_1   | 10000    | Maximum number of iterations for Step 1. |
| max_iter_2   | 10000    | Maximum number of iterations for Step 2. |
| max_iter_3a  | 10000    | Maximum number of iterations for Step 3a. |
| verbose      | FALSE    | Logical flag indicating whether to print progress during fitting. |


**Table: Arguments for `model_selec`**

| Parameter   | Default | Description |
|------------|--------|-------------|
| resultList |        | A list output from `GAR1_fit`. |
| n          |        | An integer referring to the number of observations. |
| step       | 3      | Number of steps used to fit the model (2 or 3). Requires that at least step 2 was used in `GAR1_fit`. |
| model      | "LN"   | Specifies which model was fitted: "LN" (normalized graph Laplacian), "L" (graph Laplacian), "LN.noselfloop" (normalized graph Laplacian without self-loops), "TARGAR" (time-varying GAR). |




### Value


**Table: Output of `GAR1_fit`**

| Output         | Description |
|----------------|-------------|
| S              | Matrix: The \( p \times p \) estimated covariance matrix given as supplied in GAR1_fit |
| result.L2.0    | List: A list containing the Step 1a results. The results stored in the indices of the list correspond to the index of `lambda.v`.  The result corresponding to a given `lambda.v` contains the output of the ADMM algorithm, as discussed in Harwood, Paul and Peng (2024).  |
| result.0.post  | List: A list containing the Step 2a results. The results are indexed according to the indices of `lambda.v` and `net.thre`.  The result corresponding to a given `lambda.v` and `net.thre` contains the output of the ADMM algorithm, as discussed in Harwood, Paul and Peng (2024). | 
| result.0S      | List: A list containing the Step 3 results (NULL if `step < 3`). The results are indexed according to the indices of `lambda.v` and `net.thre`. The result corresponding to a given `lambda.v` and `net.thre` contains the output of the ADMM algorithm, as discussed in Harwood, Paul and Peng (2024). |
| A.0.net        | Matrix: A matrix encoding the estimated graph topology (NULL if `step < 2`). |
| net.0.size     | Integer: Integer giving the number of edges in the estimated graph (NULL if `step < 2`). |
| v0.0.est       | Vector: A \( p \times 1 \) vector containing the estimated degree vector (NULL if `step < 3`). |
| theta0.0       | Positive scalar: estimate of \( \theta_0 \) from Step 1 (and Step 2). |
| theta0.0S      | Positive scalar: estimate of \( \theta_0 \) from Step 3 (NULL if `step < 3`). |
| conv.0.v0      | Matrix: Matrix of convergence diagnostics across (`lambda.v`, `net.thre`) combinations (NULL if `step < 3`). |



**Table: Output of `model_selec`**

| Output       | Description |
|--------------|-------------|
| model.selec  | A list containing the model selected by the eBIC criterion. The list contains the estimated (normalized) graph Laplacian `L`, `v0`, `theta0`, and `theta1`, and their associated ADMM parameters (Z, W, phi, etc.) as discussed in Harwood, Paul and Peng (2024). |
| A.net.e      | A matrix encoding the (unweighted) graph topology selected by the eBIC criterion. |
| index        | Index corresponding to the optimal tuning parameters (`lambda`, `net.thre`) for the eBIC-selected model. |
| ebic         | eBIC score of the selected model. |




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

### Model selection via eBIC for the GAR estimator
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
theta0.err.3 = abs(gar1_fitted.3$theta0 - theta0.tr)^2
L.err.3 = sum((gar1_fitted.3$theta1 * gar1_fitted.3$L - theta1.tr*LN)^2)/sum((theta1.tr*LN)^2)

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
