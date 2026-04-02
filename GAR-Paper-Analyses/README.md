# GAR-Paper-Analyses

This folder contains the companion analysis scripts used to reproduce the simulation studies and stock-market application from the GAR paper. These scripts are separate from the `SGM` package itself and are intended for reproducibility, experimentation, and worked examples.

## Contents

| File or folder | Description |
|---|---|
| `GenData.R` | Utility functions for generating random graphs, computing Laplacians, and simulating data from GAR models. This file is shared by the simulation and stock-analysis workflows. |
| `simulations_GAR_JCGS.R` | Main simulation script for the GAR paper. It generates synthetic data, fits GAR models over tuning grids, performs model selection, and computes evaluation metrics such as estimation error, power, FDR, and F1. |
| `stock_data_script.R` | Main script for the stock-market application. It preprocesses the stock data, fits GAR and GLASSO models, performs post-processing, and creates comparison plots. |
| `stockdata.rda` | Raw stock dataset used in the stock-market application. |
| `stock.data.X.Rdata` | Preprocessed stock return data created by `stock_data_processing.R`. |
| `stock-auxiliary-scripts/` | Helper scripts used by the main stock application workflow. |

## Stock Auxiliary Scripts

| File | Description |
|---|---|
| `stock-auxiliary-scripts/stock_data_processing.R` | Preprocesses the raw stock dataset, constructs the log-return matrix `X`, stores sector counts in `sp.num`, saves `stock.data.X.Rdata`, and produces exploratory plots and autocorrelation summaries. |
| `stock-auxiliary-scripts/stock_gar_res_process.R` | Post-processes fitted GAR and GLASSO models for the stock application. It computes model-size and log-likelihood summaries, performs GAR model selection, creates sector-level connectivity summaries, and writes the sector-network plot. |

## Typical Workflows

### Simulation study

Run:

```r
source("GenData.R")
source("simulations_GAR_JCGS.R")
```

The simulation script:

1. sets the graph and sample-size configuration,
2. generates GAR data replicates,
3. fits GAR models over a tuning grid,
4. applies `model_selec()` to each replicate, and
5. records recovery and estimation metrics.

### Stock application

Run:

```r
source("stock_data_script.R")
```

The stock script:

1. sources `stock-auxiliary-scripts/stock_data_processing.R`,
2. loads or creates the preprocessed stock return data,
3. fits GAR and GLASSO models,
4. sources `stock-auxiliary-scripts/stock_gar_res_process.R`, and
5. summarizes the selected models and creates comparison plots.

## Generated Outputs

Depending on which scripts you run, this folder may contain generated files such as:

| Output file | Produced by | Description |
|---|---|---|
| `stock.data.X.Rdata` | `stock_data_processing.R` | Preprocessed stock return matrix and related objects for the stock application. |
| `gar_stock.rda` | `stock_data_script.R` | Saved GAR fit object for the stock data. |
| `glasso_refit.rda` | `stock_data_script.R` | Saved refitted GLASSO results for the stock data. |
| `stock-GAR1LN-loglike-GAR-vs-glasso.pdf` | `stock_data_script.R` | Plot comparing model size and log-likelihood for GAR and GLASSO. |
| `stock-GAR1LN-sector-wise-weighted-network.pdf` | `stock_gar_res_process.R` | Sector-level connectivity plot for the selected GAR model. |

## Notes

- These scripts assume the `SGM` package and any additional analysis dependencies are installed.
- Several scripts clear the workspace with `rm(list = ...)`, so it is best to run them in a fresh R session.
- The stock workflow uses path-aware `source()` and file-loading logic so it can be run from the repository root or from inside this folder.
