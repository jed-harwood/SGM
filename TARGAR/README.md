# TARGAR

This folder is the TAR-GAR analysis area for the repository. It contains the
TAR-GAR simulation suite, reusable simulation helpers, and the temperature
application for comparing TAR-GAR, G-VAR, and graphicalVAR on year-specific
temperature data.
The graphicalVAR baseline follows Epskamp (2017), and the G-VAR baseline
follows Isufi et al. (2019).

The simulation scripts are organized around main user-facing drivers at the
top level and reusable helper functions in `simulation-auxiliary-scripts`.
The temperature application is organized separately under `temp-application`
with one user-facing runner and shared backend code for paths, defaults,
summaries, and model dispatch.

Notation in this folder uses `d` for the signal dimension, i.e. the number of
nodes, stations, or variables. In TAR-GAR(p,q), `p` is the autoregressive
order and `q` is the polynomial order in the graph Laplacian.

## Contents

| File or folder | Description |
|---|---|
| `simulations_TARGAR.R` | Main simulation script. It sets the TAR-GAR order/configuration, generates data, fits TAR-GAR models over tuning grids, computes 0S/eBIC metrics, prints summaries, and saves results. |
| `simulations_TARGAR11_on_TARGAR_pq.R` | Separate simulation script for fitting TAR-GAR(p = 1, q = 1) to data generated from the TAR-GAR(p, q) models. |
| `simulation-auxiliary-scripts/simulations_TARGAR_wrappers.R` | Helper functions for graph generation, TAR-GAR data generation, fitting replicates, computing eBIC, and extracting selected 0S/eBIC metrics. |
| `temp-application/` | Temperature application for running TAR-GAR, G-VAR, and graphicalVAR on a selected year. |
| `results/` | Default output folder created by the simulation scripts. |

## Installation

The scripts in this folder rely on the local `SGM` package for TAR-GAR fitting,
prediction, and G-VAR functionality. From the repository root, install or
reinstall the package before running either the simulations or the temperature
application:

```sh
R CMD INSTALL SGM
```

The package can also be installed from GitHub:

```r
library(devtools)
install_github(repo = "jed-harwood/SGM", subdir = "SGM")
```

## Typical Workflow

Run from the repository root:

```r
source("TARGAR/simulations_TARGAR.R")
```

Or run from inside the `TARGAR` folder:

```r
source("simulations_TARGAR.R")
```

The simulation script:

1. sets the TAR-GAR order and sample-size configuration,
2. generates the random graph, Laplacian, true TAR filters, and data replicates,
3. fits TAR-GAR models over `lambda.v` and `net.thre`,
4. computes the full BIC/eBIC paths, and
5. records the 0S/eBIC-selected metrics.

For the misspecified TAR-GAR(1,1)-on-TAR-GAR(p,q) simulations, use:

```r
source("TARGAR/simulations_TARGAR11_on_TARGAR_pq.R")
```

## Temperature Application

The temperature application lives in `TARGAR/temp-application`. It has one
user-facing script and a separate backend:

| File or folder | Description |
|---|---|
| `run_temp_application.R` | User-facing entrypoint. Pick `year` and `model` in the setup block, or pass them as command-line arguments. |
| `R/temp_application_backend.R` | Shared backend for paths, year configuration, defaults, summaries, and calls into the model pipelines. |
| `R/temp_pipeline.R` | TAR-GAR and G-VAR temperature pipeline. |
| `R/graphicalvar_temp_pipeline.R` | graphicalVAR temperature pipeline. |
| `scripts/` | Legacy one-model launchers that call the shared backend. |
| `data/temperature/` | NOAA temperature files for years 2011-2020 and combined source CSVs. |
| `data/stations/` | Station latitude/longitude files. |
| `data/knn_adj_matrices/` | Year-specific stored kNN adjacency matrices used by the G-VAR baseline. |
| `data/koppen_geiger_tif/` | Koppen-Geiger raster files and legend. |
| `results/` | Generated `.RData`, clustering summaries, and plot PDFs. |

Run the application from the repository root with:

```sh
Rscript TARGAR/temp-application/run_temp_application.R --year=2011 --model=both
```

Or run from inside `TARGAR/temp-application`:

```sh
R CMD INSTALL ../../SGM
Rscript run_temp_application.R --year=2011 --model=both
```

Available models are `both`, `targar`, and `graphicalvar`. The bundled
temperature years are 2011 through 2020. The runner validates years from files
named `data/temperature/temp_<year>.Rda`, `data/temperature/temp_<year>.csv`,
and the 2020 source file `data/temperature/dailytemp_gsod.csv`.

For TAR-GAR, the temperature application still fits the full `p = 1,2,3` by
`q = 1,2,3` grid and uses `num.pass = 3` by default.

The application expects all required files to live inside
`TARGAR/temp-application`:

| Required files | Expected location |
|---|---|
| NOAA temperature data, 2011-2019 | `data/temperature/temp_<year>.Rda` or `data/temperature/temp_<year>.csv` |
| NOAA GSOD temperature data, 2020 | `data/temperature/dailytemp_gsod.csv` |
| Station coordinates, 2011-2019 | `data/stations/latlong_<yy>.csv` or the matching file inside the year-specific kNN folder |
| Station coordinates, 2020 | Embedded in `data/temperature/dailytemp_gsod.csv` |
| G-VAR kNN adjacency matrices, 2011-2019 | `data/knn_adj_matrices/tr<year>/adjacency_matrix_k_<k>.mtx` |
| G-VAR kNN adjacency matrices, 2020 | `data/knn_adj_matrices/adjacency_matrix_k_<k>.mtx` |
| Koppen-Geiger raster | `data/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif` |

The temperature data are repository-bundled application data, not package
datasets loaded with `data()`. The NOAA files follow a station/date/temperature
schema with station identifiers, names, latitude, longitude, elevation, date,
and temperature. The backend detrends and imputes these records before fitting
TAR-GAR, G-VAR, and graphicalVAR models.

The application uses package-owned `SGM::fit_TAR_GAR()`,
`SGM::TARGAR_pred()`, and `SGM::GVAR_fit()` rather than sourcing separate local
TAR-GAR or G-VAR implementation files.

## Selecting The Simulation

Edit the user-facing setup block in `simulations_TARGAR.R`:

```r
ar.order = 1
d = 100
n = 500
q = 1
n.rep = 100
model = "LN"
theta0 = 1
theta1 = 2
edge.prob = 2 / d
num.pass = 3
graph.seed = 1
data.seed = 5 * d + n
lambda.C = c(1.5, 1, 0.5, 0.25, 0.1)
if (d == 100) {
  C.thre = exp(seq(log(1), log(0.05), length.out = 10))
} else if (d == 250) {
  C.thre = exp(seq(log(1), log(0.075), length.out = 10))
} else {
  C.thre = exp(seq(log(1), log(0.1), length.out = 10))
}
```

Use:

| `ar.order` | Default setting |
|---|---|
| `1` | TAR-GAR order 1, `q = 1`, `n = 500`, `num.pass = 3`, `stationary = TRUE`. |
| `2` | TAR-GAR order 2, `q = 3`, `n = 500`, `num.pass = 3`. |
| `3` | TAR-GAR order 3, `q = 3`, `n = 100`, `num.pass = 5`. |

All three defaults use `d = 100`, `model = "LN"`, `theta0 = 1`,
`theta1 = 2`, `edge.prob = 2 / d`, the `lambda.v` grid
`lambda.C * sqrt(log(d) / (n - fit_ar_order))`, and the `net.thre`
grid `C.thre * sqrt(log(d) / (n - fit_ar_order))`.

The original `C.thre` defaults are kept:

```r
if (d == 100) {
  C.thre = exp(seq(log(1), log(0.05), length.out = 10))
} else if (d == 250) {
  C.thre = exp(seq(log(1), log(0.075), length.out = 10))
} else {
  C.thre = exp(seq(log(1), log(0.1), length.out = 10))
}
```

The simulation scripts build the graph first so `lambda.max` is available
before the true TAR-GAR filter is specified. Edit `eta` directly using
`lambda.max`, matching the original scripts:

```r
graph.setup = prepare_targar_graph(graph.config)
lambda.max = graph.setup$lambda.max
eta = c(0.2, 0.7 / lambda.max)
```

For example, with `ar.order = 3, q = 3`, set:

```r
eta = c(-0.1, 0.1 / lambda.max, 0.1 / lambda.max^2,
        0.1 / lambda.max^3,
        -0.1, 0.1 / lambda.max, 0.12 / lambda.max^2,
        0.12 / lambda.max^3,
        -0.2, 0.1 / lambda.max, 0.2 / lambda.max^2,
        0.2 / lambda.max^3)
```

`lambda.max` is also returned as `sim.setup$lambda.max` after
`prepare_targar_simulation()`.

The data generator uses the model-specific initialization directly:

| AR order | Initialization |
|---|---|
| 1 | `Y1` is drawn from `Gamma0`. |
| 2 | The first two observations are drawn from `Gamma0.tilde`. |
| 3 | The first three observations are drawn from `Gamma0.tilde`. |

There is no burn-in. This is intentional: the modular simulations use the
model-specific initialization directly and then fit the generated `n`
observations.

For exact reproduction of the original matched TAR-GAR(1,1) forecasting
script, the setup includes `generation.n = n + 100` and then fits rows `1:n`.
Those extra rows are not used as burn-in; they only preserve the old random
number stream because the original generator drew the full innovation matrix
before drawing the initial state. For the other TAR-GAR simulations,
`generation.n = n`.

The original random seeds are kept: `graph.seed = 1` matches the original
`set.seed(1)` before graph generation, and `data.seed = 5 * d + n` matches the
original `set.seed(5*p+n)` before data generation.

## TAR-GAR(1,1) On TAR-GAR(p,q)

Use `simulations_TARGAR11_on_TARGAR_pq.R` for the simulations where the
fitted model is fixed at TAR-GAR(p = 1, q = 1), but the data-generating model
changes. In that script, edit the same manual setup variables:

```r
data.ar.order = 1
d = 100
n = 100
data.q = 3
n.rep = 100
model = "LN"
theta0 = 1
theta1 = 2
edge.prob = 2 / d
fit.ar.order = 1
fit.q = 1
graph.seed = 1
data.seed = 5 * d + n
lambda.C = c(1.5, 1, 0.5, 0.25, 0.1)
if (d == 100) {
  C.thre = exp(seq(log(1), log(0.05), length.out = 10))
} else if (d == 250) {
  C.thre = exp(seq(log(1), log(0.075), length.out = 10))
} else {
  C.thre = exp(seq(log(1), log(0.1), length.out = 10))
}
```

Use:

| `data.ar.order` | Data-generating model | Fitted model |
|---|---|---|
| `1` | TAR-GAR(p = 1, q = 3), `n = 100`. | TAR-GAR(p = 1, q = 1). |
| `2` | TAR-GAR(p = 2, q = 3), `n = 500`. | TAR-GAR(p = 1, q = 1). |
| `3` | TAR-GAR(p = 3, q = 3), `n = 500`. | TAR-GAR(p = 1, q = 1). |

These runs use the same graph, theta, eta, lambda, and threshold defaults from
the corresponding original scripts, with no burn-in.

## Metrics

The modular script focuses on the 0S/eBIC metrics from the original
simulations and stores the full eBIC path. Important outputs include:

| Output | Description |
|---|---|
| `targar.results$L.0S.ebic.err` | Relative squared error for the eBIC-selected 0S Laplacian estimate. |
| `targar.results$theta.0S.ebic.err` | Squared error for the eBIC-selected 0S `theta0` estimate. |
| `targar.results$R1.0S.ebic.err`, etc. | Relative squared error for each data-generating TAR filter `R_1, ..., R_p`. In TAR-GAR(1,1)-on-TAR-GAR(p,q) runs, unestimated higher lags are treated as zero fitted matrices. |
| `targar.results$eta.0S.ebic.err` | Squared error for selected eta coefficients. |
| `targar.results$power.0S.ebic` | Power for selected graph recovery. |
| `targar.results$fdr.0S.ebic` | FDR for selected graph recovery. |
| `targar.results$F1.0S.ebic` | F1 score for selected graph recovery. |
| `targar.results$v0.0S.ebic` | Selected `v0` estimation error. |
| `targar.results$ebic.0S` | Full 0S eBIC path over replicates, `lambda.v`, and `net.thre`. |
| `targar.results$bic.0S` | Full 0S BIC path. |
| `targar.results$ebic.0S.selec` | Selected 0S eBIC value for each replicate. |

The eBIC number of parameters is based on the fitted TAR-GAR model:

```r
fit_ar_order * (fit_q + 1) + 1 + d + net.size
```

For `simulations_TARGAR.R`, `fit_ar_order = ar_order` and `fit_q = q`.
For `simulations_TARGAR11_on_TARGAR_pq.R`, `fit_ar_order = 1` and
`fit_q = 1`. This is centralized in `targar_ebic_parameter_count()` in the
wrapper file.

## Generated Outputs

By default, both simulation scripts save:

| Output file | Description |
|---|---|
| `results/<case>_d<d>_n<n>_model<model>_rep<n.rep>_modular.RData` | Saved list containing `setup`, `results.TARGAR`, `targar.results`, and `summary.tables`. |

Set `save.results = FALSE` in `simulations_TARGAR.R` to run without writing an
output file. Set `keep.fits = FALSE` to save only setup, metrics, and summaries.

## Notes

- These scripts assume the local `SGM` package is installed. TAR-GAR,
  prediction, and G-VAR helpers are loaded from the package exports.
- Several runs are computationally expensive with the original defaults
  (`d = 100`, `n.rep = 100`). For quick checks, temporarily reduce `d`,
  `n`, `n.rep`, and `num.thread` in the setup block.
- The scripts clear the workspace with `rm(list = ls())`, so it is best to run
  them in a fresh R session.

## References

- Epskamp, S. (2017). graphicalVAR: Graphical VAR for experience sampling
  data. CRAN package documentation.
  https://cran.r-project.org/web/packages/graphicalVAR/graphicalVAR.pdf
- Isufi, E., Loukas, A., Perraudin, N., and Leus, G. (2019). Forecasting Time
  Series With VARMA Recursions on Graphs. IEEE Transactions on Signal
  Processing, vol. 67, no. 18, pp. 4870-4885, 15 Sept. 2019.
  https://doi.org/10.1109/TSP.2019.2929930
