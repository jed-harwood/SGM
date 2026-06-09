# SGM Package

The `SGM` directory contains the R package for fitting Spectral Graph Models.
It provides GAR(1), TAR-GAR, model selection, forecasting, G-VAR helper, and
goodness-of-fit functions, along with the package datasets and documentation.
The G-VAR helper follows the graph VARMA forecasting framework of Isufi et al.
(2019).

## Installation

Install the package dependencies first:

```r
install.packages(c("doParallel", "foreach", "gmp", "mnormt", "quadprog", "nloptr"))
```

From the repository root, install the local package with:

```sh
R CMD INSTALL SGM
```

Or install from GitHub with:

```r
library(devtools)
install_github(repo = "jed-harwood/SGM", subdir = "SGM")
```

## Datasets

Load package datasets with `data("<dataname>")`.

| Dataset | Description |
|---|---|
| `gar1` | A sample simulated from an underlying GAR(1) model, with `d = 100` nodes and `n = 100` observations. |
| `stocks` | Standardized log-returns from 283 stocks on the S&P 500, spanning January 1, 2007, to January 1, 2011. |
| `targar` | A sample simulated from an underlying TAR-GAR model, with `d = 100` nodes and `n = 100` observations. |

For more information on a dataset, run `?<dataname>` in R.

## Main Functions

| Function | Purpose |
|---|---|
| `GAR1_fit()` | Fit a GAR(1) model over `lambda.v` and `net.thre` tuning grids. |
| `model_selec()` | Select a fitted GAR or TAR-GAR model via eBIC. |
| `TARGAR_fit()` | Fit a TAR-GAR(p,q) model for `p = 1`, `p = 2`, or `p = 3`. |
| `fit_TAR_GAR()` | Compatibility wrapper for script-style TAR-GAR usage. |
| `GAR1_gf()` | Compute the GAR(1) goodness-of-fit measure. |
| `TARGAR_pred()` | Forecast with a selected TAR-GAR model. |
| `ARp_pred()` | Forecast with a standard AR(p) baseline. |
| `GVAR_fit()` | Fit the G-VAR baseline used by the temperature application. |
| `G.VAR.Fit()` | Compatibility alias for `GVAR_fit()`. |

## Examples

```r
library(SGM)

data("gar1")
S <- var(gar1$data) * (nrow(gar1$data) - 1) / nrow(gar1$data)
nobs <- nrow(gar1$data)
d <- ncol(gar1$data)

lambda.v <- c(1, 0.5) * sqrt(log(d) / nobs)
net.thre <- exp(seq(log(1), log(0.05), length.out = 10)) * sqrt(log(d) / nobs)
rho.v <- pmax(lambda.v, 0.01)

fit <- GAR1_fit(S, nobs, lambda.v, net.thre, model = "LN", step = 3, rho.v = rho.v)
sel <- model_selec(fit)

sel$theta0
sel$theta1
sel$L
sel$v0
```

The repository-level `README.md` contains the full argument and output tables
for these functions.

## References

- Isufi, E., Loukas, A., Perraudin, N., and Leus, G. (2019). Forecasting Time
  Series With VARMA Recursions on Graphs. IEEE Transactions on Signal
  Processing, vol. 67, no. 18, pp. 4870-4885, 15 Sept. 2019.
  https://doi.org/10.1109/TSP.2019.2929930
