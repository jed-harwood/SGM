# SGM
An R package for fitting Spectral Graph Models

*** 


## Contents
- [Overview](#Overview)
- [Installation](#Installation)
- [Repo Layout](#Repo-Layout)
- [Main Functions](#Main-Functions)
- [Datasets](#Datasets)
- [Examples](#Examples)
- [References](#References)
- [Contact](#Contact)

*** 

## Overview

This repository contains the `SGM` R package for fitting Spectral Graph Models
and two companion analysis areas:

* `GAR`: simulation and stock-market analyses from the GAR paper.
* `TARGAR`: TAR-GAR simulations and the temperature application.

Most users will work with the `SGM` package directly:

1. Compute a covariance matrix `S` from data.
2. Fit a GAR model with `GAR1_fit()`.
3. Select a final model with `model_selec()`.
4. Extract `theta0`, `theta1`, `L`, and `v0` from the selected model.

The `GAR` folder builds on that package API with reproducible simulation
studies and a stock-market application for the GAR paper.
The GAR analysis documentation cites the graph VARMA forecasting paper of
Isufi et al. (2019).

For temporally dependent signals with explicit autoregressive graph filters,
use `TARGAR_fit()` instead. The TAR-GAR workflow follows the same package
rhythm: fit a path over `lambda.v` and `net.thre`, select with
`model_selec()`, then inspect the selected Laplacian and AR filter matrices.
The `TARGAR` folder builds on that package API with reproducible simulations
and a temperature application that compares TAR-GAR, G-VAR, and graphicalVAR
on year-specific temperature data.
The graphicalVAR baseline follows Epskamp (2017), and the G-VAR baseline
follows the graph VARMA forecasting framework of Isufi et al. (2019).

## Installation 

### Dependencies
Please make sure that the following packages are installed before using the R package `SGM`. 

```
install.packages(c("doParallel", "foreach", "gmp", "mnormt", "quadprog", "nloptr"))
```

### Installation from the Repository
The R Package `SGM` can be installed from source code in this GitHub repository.  The `devtools` package is required for installation.  
```
library(devtools)
install_github(repo="jed-harwood/SGM", subdir="SGM")
```

From a local checkout of this repository, install the package with:

```sh
R CMD INSTALL SGM
```

The temperature application in `TARGAR/temp-application` also relies on the installed `SGM` package for TAR-GAR, prediction, and G-VAR functions. When running from inside that application folder, reinstall the local package with:

```sh
R CMD INSTALL ../../SGM
```

## Repo Layout

This repository has three main folders.

* `SGM`: The R package itself. It contains the user-facing fitting functions, documentation, C++ solvers, datasets, and examples.
* `GAR`: Companion analysis scripts used to reproduce the simulation studies and stock-market application from the GAR paper.
* `TARGAR`: TAR-GAR simulation scripts, reusable TAR-GAR experiment helpers, and the temperature application.

### `SGM` Package

The `SGM` package folder contains:

* user-facing fitting functions for GAR(1), TAR-GAR, model selection, forecasting, and goodness-of-fit calculations,
* package documentation and datasets,
* C++ solvers used by the ADMM routines, and
* package examples.

### `GAR` Analyses

The `GAR` folder contains:

* `GenData.R`, reusable GAR graph generation, Laplacian, and data-simulation utilities,
* `simulations_GAR_JCGS.R`, the main GAR paper simulation driver,
* `stock_data_script.R`, the stock-market application driver,
* `stock-auxiliary-scripts/`, preprocessing and post-processing utilities for the stock application,
* bundled stock data and derived stock-return data files, and
* generated GAR, GLASSO, and stock-network outputs under the `GAR` folder.

The `GAR` folder is separate from the `SGM` package itself and is intended for reproducibility and worked examples rather than the package API.

### `TARGAR` Analyses And Temperature Application

The `TARGAR` folder contains:

* `simulations_TARGAR.R`, the main TAR-GAR simulation driver,
* `simulations_TARGAR11_on_TARGAR_pq.R`, the misspecified TAR-GAR(1,1)-on-TAR-GAR(p,q) simulation driver,
* `simulation-auxiliary-scripts/`, reusable helpers for TAR-GAR graph generation, data generation, fitting, eBIC computation, and metrics,
* `temp-application/`, a year-selectable temperature application for TAR-GAR, G-VAR, and graphicalVAR analyses,
* bundled temperature data, station files, kNN adjacency matrices, and Koppen-Geiger raster assets used by the temperature application, and
* generated simulation or application outputs under local `results/` folders.

The `TARGAR` folder is separate from the `SGM` package itself. It is intended
for TAR-GAR reproducibility, worked examples, and the applied temperature
workflow rather than as the package API.

## Main Functions

`GAR1_fit`: fit a `GAR(1)` model for a given set of tuning parameters, using a 3-step estimation procedure based on the penalized MLE.  
* Step 0: given an initial estimate `S` for the covariance matrix (e.g., the sample covariance matrix), obtain an initial estimate for `theta0` by the reciprocal of the largest eigenvalue of `S`, square-rooted.
* Step 1: given `S` and `theta0` from Step 0, use an ADMM algorithm to estimate the (normalized) graph Laplacian `L` of the latent graph and extract the thresholded zero-pattern used in later steps.
* Step 2: given the zero-pattern from Step 1, use an ADMM algorithm to refit the non-zero elements of the (normalized) graph Laplacian `L`.
* Step 3:
    * a. given the (normalized) graph Laplacian `L` from Step 2, use an ADMM algorithm to estimate the normalized degree vector `v0`.
    * b. given `v0` and the zero-pattern from Step 1, use an ADMM algorithm to simultaneously re-estimate `theta0` and the (normalized) graph Laplacian `L`.

`model_selec`: conduct model selection via the eBIC criterion.

`TARGAR_fit`: fit a `TAR-GAR(p,q)` model for a given set of tuning parameters. Here `p` is the AR order, not the dimension of the data. The packaged TAR-GAR implementation supports `p = 1`, `p = 2`, and `p = 3`; `q` is the polynomial order in the graph Laplacian. The defaults match the original hard-coded pipeline: `p = 1`, `q = 1`, `num.pass = 3`, `model = "LN"`, `eps.thre = 1e-6`, `eps_abs = 1e-5`, `eps_rel = 1e-3`, `max_iter = 50000`, `deg_max_iter = 50000`, `eta_max_iter = 1000`, and `stationary = TRUE`.

`GAR1_gf`: calculate a goodness-of-fit measure to determine if `GAR(1)` is an appropriate model for the data, as proposed in Harwood, Paul, and Peng (2024). Valid when the sample size is at least the signal dimension, i.e. $n \geq d$. Return a value between $0$ and $1$. If the returned value is close to $1$, then it means that the `GAR(1)` model is a good fit to the data.

`TARGAR_pred`: forecast future observations from a selected TAR-GAR fit.

`ARp_pred`: forecast future observations with a standard AR(p) baseline.

`GVAR_fit`: fit the G-VAR baseline used by the temperature application. The compatibility wrapper `G.VAR.Fit()` is also exported for older scripts.

Typical workflow:

```r
fit <- GAR1_fit(S, nobs, lambda.v, net.thre, model = "LN")
sel <- model_selec(fit)

sel$theta0
sel$theta1
sel$L
sel$v0
```

TAR-GAR workflow:

```r
fit <- TARGAR_fit(data, lambda.v, net.thre, rho.v = rho.v,
                  p = 1, q = 1, model = "LN")
sel <- model_selec(fit)

sel$theta0
sel$theta1
sel$L
sel$R_list
sel$eta
```

The compatibility wrapper `fit_TAR_GAR()` is also exported for script-style
usage:

```r
fit <- fit_TAR_GAR(data, p = 1, q = 1,
                   lambda.v = lambda.v, net.thre = net.thre)
```

### Arguments

**Table: Arguments for `GAR1_fit`**

| Parameter | Data type | Default | Description |
|-----------|-----------|---------|-------------|
| S | matrix |  | Estimate of the covariance matrix (e.g., the MLE). |
| nobs | integer |  | Sample size used to compute `S`. |
| lambda.v | numeric vector |  | Tuning parameter controlling sparsity of the estimated graph. |
| net.thre | numeric vector |  | Thresholding parameter used in Step 1 to remove noisy edges and define the zero-patterns used in later steps. |
| model | character | "LN" | Type of graph Laplacian: "LN" (normalized Laplacian), "L" (graph Laplacian), or "LN.noselfloop" (normalized Laplacian without self-loops). |
| step | integer | 3 | Number of steps in the estimation procedure (1, 2, or 3). |
| rho.v | numeric vector | lambda.v | ADMM parameter (typically set equal to `lambda.v`). |
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
| resultList | list |  | A list output from `GAR1_fit` or `TARGAR_fit`. The required metadata are read directly from this object. |


**Table: Arguments for `TARGAR_fit`**

| Parameter | Data type | Default | Description |
|-----------|-----------|---------|-------------|
| data | matrix |  | Time-ordered data matrix with observations in rows. |
| lambda.v | numeric vector |  | Tuning parameter controlling sparsity of the estimated graph. |
| net.thre | numeric vector |  | Thresholding parameter used to remove noisy edges and define final refit zero-patterns. |
| rho.v | numeric vector | lambda.v | ADMM parameter sequence, typically set to `lambda.v` or `pmax(lambda.v, 0.01)`. |
| num.pass | integer | 3 | Number of alternating TAR-GAR passes before final refitting. |
| model | character | "LN" | Innovation graph model: "LN", "L", or "LN.noselfloop". |
| p | integer | 1 | AR order; must be 1, 2, or 3. |
| q | integer | 1 | Polynomial order in the graph Laplacian; must be positive. |
| eps.thre | numeric scalar | 1e-6 | Small positive threshold used for numerical stability and stationarity constraints. |
| eps_abs | numeric scalar | 1e-5 | Absolute tolerance for ADMM convergence. |
| eps_rel | numeric scalar | 1e-3 | Relative tolerance for ADMM convergence. |
| max_iter | integer | 50000 | Maximum number of ADMM iterations. |
| deg_max_iter | integer | 50000 | Maximum number of iterations for the degree-vector ADMM refit. |
| lap_z_max_iter | integer | max_iter | Maximum number of inner Z updates in the Laplacian refit. |
| eta_max_iter | integer | 1000 | Maximum number of eta optimization iterations for the `p = 3` nonlinear program. |
| stationary | logical | TRUE | Whether to impose stationarity constraints while fitting TAR filters. |
| verbose | logical | FALSE | Logical flag indicating whether to print progress during fitting. |


**Table: Arguments for `GAR1_gf`**

| Parameter    | Default   | Description |
|--------------|----------|-------------|
| S            |          | Estimate of the covariance matrix (e.g., the MLE). |
| nobs         |          | Number of observations used to compute `S`. |
| lambda.v     |          | Positive tuning parameter for the GAR(1) model. This should typically be set to `sqrt(log(d)/nobs)` where `d` is the signal dimension and `nobs` is the sample size. |
| rho.v        | lambda.v | ADMM penalty parameter (positive; typically set equal to `lambda.v`). |
| eps_thre     | 1e-6     | Small positive threshold used for numerical stability. |
| eps_abs      | 1e-5     | Absolute tolerance for ADMM convergence. |
| eps_rel      | 1e-3     | Relative tolerance for ADMM convergence. |
| max_iter     | 10000    | Maximum number of iterations for fitting on observed data. |
| num.thread   | 1        | Number of threads used for computation. |
| rep.boot     | 100      | Number of bootstrap samples used for the goodness-of-fit test. |
| seed         | 1        | Random seed for reproducibility. |




### Value


**Table: Output of `GAR1_fit`**

| Output | Data type | Description |
|--------|-----------|-------------|
| S | matrix | The \( d \times d \) covariance matrix supplied to `GAR1_fit`, where `d` is the signal dimension. |
| nobs | integer | Number of observations used to compute `S`. |
| model | character | The fitted GAR model family. |
| step | integer | The last step requested in the fitting procedure. |
| lambda.v | numeric vector | Sparsity tuning parameter sequence. |
| rho.v | numeric vector | ADMM tuning parameter sequence. |
| net.thre | numeric vector | Threshold sequence used in Steps 2 and 3. |
| theta0.init | numeric scalar | Initial estimate of \( \theta_0 \) from Step 0, used in Steps 1 and 2. |
| theta0.s3 | numeric matrix | Matrix of Step 3b estimates for \( \theta_0 \) (NULL if `step < 3`). |
| A.net | list | Estimated graph topologies indexed by `lambda.v` and `net.thre`, created in Step 1 and used in Steps 2 and 3. |
| step1 | list | Step 1 ADMM fit objects indexed by `lambda.v`. Each element is the returned output of `ADMM_L2`, including `L`, `Z`, `W`, `theta0`, `theta1`, and `conv`. Step 1 also defines the thresholded zero-patterns stored in `A.net`. |
| step2 | list | Step 2 refit objects indexed by `lambda.v` and `net.thre` (NULL if `step < 2`). Each element is the returned output of `ADMM_L2_Zero`, containing the refitted `L`, `Z`, `W`, `theta0`, `theta1`, and `conv` for the corresponding zero-pattern from Step 1. |
| step3a | list | Step 3a `v0` estimation results indexed by `lambda.v` and `net.thre` (NULL if `step < 3`). Each element stores the estimated `v0` vector computed from the Step 2 Laplacian estimate. |
| step3b | list | Step 3b joint estimation objects indexed by `lambda.v` and `net.thre` (NULL if `step < 3`). Each element is the returned output of `ADMM_Lap_Zero`, containing the jointly refitted `L`, `theta0`, `theta1`, `Z`, `W`, `phi`, and `conv`. |
| v0.s3 | list | Step 3a `v0` estimates indexed by `lambda.v` and `net.thre` (NULL if `step < 3`). |
| conv | list | Convergence diagnostics stored as `conv$step1`, `conv$step2`, `conv$step3a`, and `conv$step3b`. |

For most applications, you do not need to inspect `step1`, `step2`, `step3a`, or `step3b` directly. The usual workflow is to fit with `GAR1_fit()`, select with `model_selec()`, and then work with the top-level `theta0`, `theta1`, `L`, and `v0` returned by `model_selec()`.



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
| R_list | list | For TAR-GAR fits, selected AR filter matrices `R1`, ..., `R_p`, one for each AR lag. |
| eta | numeric vector | For TAR-GAR fits, selected polynomial filter coefficients. |
| p, q | integer | For TAR-GAR fits, selected model AR order and polynomial order metadata. |


**Table: Output of `GAR1_gf`**

| Output  | Description |
|---------|-------------|
| p-value | P-value for the goodness-of-fit test. |


## Datasets

This repository contains package datasets and application data. Package
datasets live in `SGM` and can be loaded into the working environment with
`data("<dataname>")`.


* `gar1`: A sample simulated from an underlying GAR(1) model, with `d = 100` nodes and `n = 100` observations. The underlying graph was a random graph with edge probability 0.02.
* `stocks`: A collection of (standardized) log-returns from 283 stocks on the S\&P 500.  The dataset spans January 1, 2007, to January 1, 2011, covering the global financial crisis, with 1007 closing prices per stock.  The stocks come from five GICS sectors: 58 from Information Technology, 72 from Consumer Discretionary, 32 from Consumer Staples, 59 from Financials, and 62 from Industrials.
* `targar`: A sample simulated from an underlying TAR-GAR model, with `d = 100` nodes and `n = 100` observations.

These datasets can be accessed as the following R objects.

`gar1`:  A list object that contains:
1. `data`: data simulated under a `GAR(1)` model.
2. `A.tr`: the true weighted adjacency matrix for the underlying graph.
3. `LN`: the true normalized graph Laplacian matrix.
4. `theta0`, `theta1`: the true graph filter parameters.  


`stocks`: A 1007 by 283 matrix.

`targar`: A list object that contains:
1. `data`: data simulated under a TAR-GAR model.
2. `A.tr`: the true weighted adjacency matrix for the underlying graph.
3. `LN`: the true normalized graph Laplacian matrix.
4. `theta0`, `theta1`: the true innovation graph filter parameters.
5. `eta0`, `eta1`: the true TAR filter parameters for the included example.


For more information on a given dataset, please run `?<dataname>`.  

The TAR-GAR temperature application also includes year-specific temperature
application data under `TARGAR/temp-application/data`. These files are not
installed as package datasets; they are read by the temperature application
backend using paths inside the repository. The temperature records were
retrieved from NOAA and focus on California stations for years 2011 through
2020.

| Application data | Location | Notes |
|---|---|---|
| NOAA temperature data, 2011-2015 | `TARGAR/temp-application/data/temperature/temp_<year>.Rda` and `temp_CA_2010_2015.csv` | Year-specific processed matrices are bundled for 2011-2015; the combined CSV is retained as source data for the earlier years. |
| NOAA temperature data, 2016-2019 | `TARGAR/temp-application/data/temperature/temp_<year>.csv` and, where available, `temp_<year>.Rda` | Raw station-date-temperature CSVs plus processed year matrices for several years. |
| NOAA GSOD temperature data, 2020 | `TARGAR/temp-application/data/temperature/dailytemp_gsod.csv` | The 2020 run uses this source CSV directly; station coordinates are embedded in the file. |
| Station coordinates | `TARGAR/temp-application/data/stations/latlong_<yy>.csv` | Latitude/longitude files for the 2011-2019 station sets. |
| kNN adjacency matrices | `TARGAR/temp-application/data/knn_adj_matrices/` | Stored sparse adjacency matrices used by the G-VAR and graphicalVAR temperature workflows. |
| Koppen-Geiger rasters | `TARGAR/temp-application/data/koppen_geiger_tif/` | Climate-zone rasters used for temperature-application summaries and plots. |



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
d = ncol(gar1_data)

### Set model to fit
model = "LN"

### Set tuning parameters: lambda and net.thre sequence
C.v=c(1,0.5)  
lambda.v=C.v*sqrt(log(d)/nobs)

C.thre=exp(seq(log(1),log(0.05), length.out=10))
net.thre=C.thre*sqrt(log(d)/nobs)

### Set ADMM parameter 
rho.v=pmax(lambda.v, 0.01)

### Get sample covariance 
S = var(gar1_data)*(nobs-1)/nobs

### Fit GAR(1) (up to step 3)
fit = GAR1_fit(S, nobs, lambda.v, net.thre, model, 3, rho.v)

### Model selection via eBIC for the GAR estimator
sel = model_selec(fit)

### Evaluation: estimation errors for theta0, theta1*L, and FDR and Power for graph inference 
## Get ground truth
A.tr = gar1$A.tr > 0 # True 0-1 adjacency matrix
diag(A.tr)=0
LN = gar1$LN # True (normalized) graph Laplacian
theta0.tr = gar1$theta0 # True graph filter parameter
theta1.tr = gar1$theta1 # True graph filter parameter

## Calculate estimation errors
theta0.err = abs(sel$theta0 - theta0.tr)^2
L.err = sum((sel$theta1 * sel$L - theta1.tr*LN)^2)/sum((theta1.tr*LN)^2)

## Calculate FDR and Power 
FDR = sum(sel$A.net.e*(1-A.tr))/sum(sel$A.net.e)
Power = sum(sel$A.net.e*A.tr)/sum(A.tr)

## Network sizes
net.size = sum(sel$A.net.e)/2

## View results 
c(theta0.err) # theta0 errors
c(L.err) # L errors
c(net.size) # network sizes
c(Power) # Power
c(FDR) # FDR

## Results:
## d=100; n=100; true network size = 105; 
## theta0: 0.00967012
## L: 0.01916969
## Sizes: 103
## Power: 0.9238095
## FDR: 0.05825243
```


### Fit `TAR-GAR(p,q)` to the simulated `targar` data.

```
library(SGM)

data("targar")

targar_data = targar$data
nobs = nrow(targar_data)
d = ncol(targar_data)

model = "LN"
ar.order = 1
poly.order = 1

C.v = c(1, 0.5)
lambda.v = C.v * sqrt(log(d) / (nobs - ar.order))
rho.v = pmax(lambda.v, 0.01)

C.thre = exp(seq(log(1), log(0.1), length.out = 5))
net.thre = C.thre * sqrt(log(d) / (nobs - ar.order))

fit = TARGAR_fit(
  targar_data,
  lambda.v = lambda.v,
  net.thre = net.thre,
  rho.v = rho.v,
  p = ar.order,
  q = poly.order,
  num.pass = 3,
  model = model
)

sel = model_selec(fit)

sel$theta0
sel$theta1
sel$L
sel$R_list
sel$eta
```

For higher-order TAR-GAR models, set `p = 2` or `p = 3`. The `p = 3` path
uses the nonlinear eta optimizer from `nloptr`; the `p = 1` and `p = 2`
paths use quadratic programs from `quadprog`.



### Compute the `GAR(1)` goodness-of-fit measure.

```
### load libraries
library(SGM)

### See ?gar1 for more information
data("gar1")

### Get data 
gar1_data = gar1$data
nobs = nrow(gar1_data)
d = ncol(gar1_data)

### Set model to fit
model = "LN"

### Set tuning parameters: lambda and net.thre sequence
C.v=c(1,0.5)  
lambda.v=C.v*sqrt(log(d)/nobs)

C.thre=exp(seq(log(1),log(0.05), length.out=10))
net.thre=C.thre*sqrt(log(d)/nobs)

### Set ADMM parameter 
rho.v=pmax(lambda.v, 0.01)

### Get sample covariance 
S = var(gar1_data)*(nobs-1)/nobs

### Goodness-of-fit measure
GAR1_gf(S, nobs, lambda.v[1], num.thread = 2)

# > 1

```

***

## References

- Epskamp, S. (2017). graphicalVAR: Graphical VAR for experience sampling
  data. CRAN package documentation.
  https://cran.r-project.org/web/packages/graphicalVAR/graphicalVAR.pdf
- Isufi, E., Loukas, A., Perraudin, N., and Leus, G. (2019). Forecasting Time
  Series With VARMA Recursions on Graphs. IEEE Transactions on Signal
  Processing, vol. 67, no. 18, pp. 4870-4885, 15 Sept. 2019.
  https://doi.org/10.1109/TSP.2019.2929930

***

## Contact
Please report any bugs to `jedharwood@ucdavis.edu`
