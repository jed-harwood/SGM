# Temperature Application

This directory contains the temperature experiments for TAR-GAR, G-VAR, and
graphicalVAR. The application is organized around one user-facing script and a
shared backend so that users only need to choose the year and model family.

## Layout

- `run_temp_application.R`: user-facing entrypoint. Pick `year` and `model` here, or pass them as command-line arguments.
- `R/temp_application_backend.R`: shared backend for paths, year-specific files, defaults, summaries, and model dispatch.
- `R/temp_pipeline.R`: TAR-GAR and G-VAR temperature pipeline.
- `R/graphicalvar_temp_pipeline.R`: graphicalVAR temperature pipeline.
- `scripts/`: legacy one-model launchers that call the shared backend.
- `data/temperature/`: yearly temperature files and combined source CSVs.
- `data/stations/`: station latitude/longitude files.
- `data/knn_adj_matrices/`: year-specific stored kNN adjacency matrices.
- `data/koppen_geiger_tif/`: Koppen-Geiger raster files and legend.
- `results/`: generated `.RData`, clustering summaries, and plot PDFs.

## Installation

Install or reinstall the local `SGM` package first so the application uses
package-owned TAR-GAR, prediction, and G-VAR functions:

```sh
R CMD INSTALL ../../SGM
```

From the repository root, the equivalent command is:

```sh
R CMD INSTALL SGM
```

## Run

```sh
Rscript run_temp_application.R --year=2011 --model=both
```

Available models are `both`, `targar`, and `graphicalvar`.

You can also edit the setup block in `run_temp_application.R`:

```r
year <- 2011
model <- "both"
save_results <- TRUE
```

Command-line arguments override those values. The bundled temperature years are
2011 through 2020, detected from `data/temperature/temp_<year>.Rda`,
`data/temperature/temp_<year>.csv`, and the 2020 source file
`data/temperature/dailytemp_gsod.csv`.

## Required Files

All paths are relative to this directory.

| File type | Path pattern |
|---|---|
| Temperature data, 2011-2019 | `data/temperature/temp_<year>.Rda` or `data/temperature/temp_<year>.csv` |
| Temperature data, 2020 | `data/temperature/dailytemp_gsod.csv` |
| Combined 2010-2015 source data | `data/temperature/temp_CA_2010_2015.csv` |
| Station coordinates, 2011-2019 | `data/stations/latlong_<yy>.csv` or `data/knn_adj_matrices/tr<year>/latlong_<yy>.csv` |
| Station coordinates, 2020 | Embedded in `data/temperature/dailytemp_gsod.csv` |
| kNN adjacency matrices, 2011-2019 | `data/knn_adj_matrices/tr<year>/adjacency_matrix_k_<k>.mtx` |
| kNN adjacency matrices, 2020 | `data/knn_adj_matrices/adjacency_matrix_k_<k>.mtx` |
| Koppen-Geiger raster | `data/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif` |

The backend validates these paths before loading year-specific data. Generated
outputs are written to `results/` by default.
