# Temperature Application

This directory contains the temperature experiments for TAR-GAR, G-VAR, and graphicalVAR.

## Layout

- `run_temp_application.R`: user-facing entrypoint. Pick `year` and `model` here, or pass them as command-line arguments.
- `R/`: backend and model pipeline code.
- `scripts/`: legacy one-model launchers that call the shared backend.
- `data/temperature/`: yearly temperature files and combined source CSVs.
- `data/stations/`: station latitude/longitude files.
- `data/knn_adj_matrices/`: year-specific stored kNN adjacency matrices.
- `data/koppen_geiger_tif/`: Koppen-Geiger raster files and legend.
- `results/`: generated `.RData`, clustering summaries, and plot PDFs.

## Run

Install or reinstall the local `SGM` package first so the application uses
package-owned TAR-GAR, prediction, and G-VAR functions:

```sh
R CMD INSTALL ../../SGM
```

```sh
Rscript run_temp_application.R --year=2011 --model=both
```

Available models are `both`, `targar`, and `graphicalvar`.
