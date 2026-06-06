# TARGAR

This folder contains modular simulation scripts for the TAR-GAR experiments.
The scripts mirror the user-facing structure of `GAR-Paper-Analyses`: a main
simulation driver at the top level and reusable helper functions in an
auxiliary script folder.

## Contents

| File or folder | Description |
|---|---|
| `simulations_TARGAR.R` | Main simulation script. It sets the TAR-GAR order/configuration, generates data, fits TAR-GAR models over tuning grids, computes 0S/eBIC metrics, prints summaries, and saves results. |
| `simulation-auxiliary-scripts/simulations_TARGAR_wrappers.R` | Helper functions for graph generation, TAR-GAR data generation, fitting replicates, computing eBIC, and extracting selected 0S/eBIC metrics. |
| `results/` | Default output folder created by `simulations_TARGAR.R`. |

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

## Selecting The Simulation

Edit the user-facing setup block in `simulations_TARGAR.R`:

```r
ar.order = 1
config = targar_default_config(ar.order)
```

Use:

| `ar.order` | Default setting |
|---|---|
| `1` | TAR-GAR order 1, `q = 1`, `n = 500`, `num.pass = 3`, `stationary = FALSE`. |
| `2` | TAR-GAR order 2, `q = 3`, `n = 500`, `num.pass = 3`. |
| `3` | TAR-GAR order 3, `q = 3`, `n = 100`, `num.pass = 5`. |

All three defaults use `d = 100`, `model = "LN"`, `theta0 = 1`,
`theta1 = 2`, `edge.prob = 2 / d`, the `lambda.v` grid
`c(1.5, 1, 0.5, 0.25, 0.1) * sqrt(log(d) / (n - ar.order))`, and the
same `net.thre` constant grids used in the original simulation scripts.

The data generator uses the model-specific initialization directly:

| AR order | Initialization |
|---|---|
| 1 | `Y1` is drawn from `Gamma0`. |
| 2 | The first two observations are drawn from `Gamma0.tilde`. |
| 3 | The first three observations are drawn from `Gamma0.tilde`. |

There is no burn-in. This is intentional: the modular simulations use the
model-specific initialization directly and then fit the generated `n`
observations.

## Metrics

The modular script focuses on the 0S/eBIC metrics from the original
simulations and stores the full eBIC path. Important outputs include:

| Output | Description |
|---|---|
| `targar.results$L.0S.ebic.err` | Relative squared error for the eBIC-selected 0S Laplacian estimate. |
| `targar.results$theta.0S.ebic.err` | Squared error for the eBIC-selected 0S `theta0` estimate. |
| `targar.results$R1.0S.ebic.err`, etc. | Relative squared error for selected TAR filters. Higher-order runs include `R2.0S.ebic.err` and `R3.0S.ebic.err` as appropriate. |
| `targar.results$eta.0S.ebic.err` | Squared error for selected eta coefficients. |
| `targar.results$power.0S.ebic` | Power for selected graph recovery. |
| `targar.results$fdr.0S.ebic` | FDR for selected graph recovery. |
| `targar.results$F1.0S.ebic` | F1 score for selected graph recovery. |
| `targar.results$v0.0S.ebic` | Selected `v0` estimation error. |
| `targar.results$ebic.0S` | Full 0S eBIC path over replicates, `lambda.v`, and `net.thre`. |
| `targar.results$bic.0S` | Full 0S BIC path. |
| `targar.results$ebic.0S.selec` | Selected 0S eBIC value for each replicate. |

The eBIC number of parameters is:

```r
ar_order * (q + 1) + 1 + d + net.size
```

This is centralized in `targar_ebic_parameter_count()` in the wrapper file.

## Generated Outputs

By default, `simulations_TARGAR.R` saves:

| Output file | Description |
|---|---|
| `results/<case>_d<d>_n<n>_model<model>_rep<n.rep>_modular.RData` | Saved list containing `setup`, `results.TARGAR`, `targar.results`, and `summary.tables`. |

Set `save.results = FALSE` in `simulations_TARGAR.R` to run without writing an
output file. Set `keep.fits = FALSE` to save only setup, metrics, and summaries.

## Notes

- These scripts assume the local `SGM` package is installed and that the TAR-GAR
  package functions are available.
- Several runs are computationally expensive with the original defaults
  (`d = 100`, `n.rep = 100`). For quick checks, temporarily reduce `d`,
  `n`, `n.rep`, and `num.thread` in the setup block.
- The scripts clear the workspace with `rm(list = ls())`, so it is best to run
  them in a fresh R session.
