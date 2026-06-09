# Temperature Application

This directory contains the temperature experiments for TAR-GAR, G-VAR, and
graphicalVAR. The application is organized around one user-facing script and a
shared backend so that users only need to choose the year and model family.
The graphicalVAR baseline follows Epskamp (2017), and the G-VAR baseline
follows Isufi et al. (2019).

Notation in this application uses `d` for the number of stations. In
TAR-GAR(p,q), `p` is the autoregressive order and `q` is the polynomial order
in the graph Laplacian.

## Layout

- `run_temp_application.R`: user-facing entrypoint. Pick `year` and `model` here, or pass them as command-line arguments.
- `R/temp_application_backend.R`: shared backend for paths, year-specific files, defaults, summaries, and model dispatch.
- `R/temp_pipeline.R`: TAR-GAR and G-VAR temperature pipeline.
- `R/graphicalvar_temp_pipeline.R`: graphicalVAR temperature pipeline.
- `scripts/`: legacy one-model launchers that call the shared backend.
- `data/temperature/`: NOAA temperature files for years 2011-2020 and combined source CSVs.
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

The TAR-GAR temperature pipeline still fits the full `p = 1,2,3` by `q = 1,2,3`
grid and uses `num.pass = 3` by default.

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

## Datasets And Required Files

All paths are relative to this directory.

| File type | Path pattern |
|---|---|
| NOAA temperature data, 2011-2019 | `data/temperature/temp_<year>.Rda` or `data/temperature/temp_<year>.csv` |
| NOAA GSOD temperature data, 2020 | `data/temperature/dailytemp_gsod.csv` |
| Combined 2010-2015 source data | `data/temperature/temp_CA_2010_2015.csv` |
| Station coordinates, 2011-2019 | `data/stations/latlong_<yy>.csv` or `data/knn_adj_matrices/tr<year>/latlong_<yy>.csv` |
| Station coordinates, 2020 | Embedded in `data/temperature/dailytemp_gsod.csv` |
| kNN adjacency matrices, 2011-2019 | `data/knn_adj_matrices/tr<year>/adjacency_matrix_k_<k>.mtx` |
| kNN adjacency matrices, 2020 | `data/knn_adj_matrices/adjacency_matrix_k_<k>.mtx` |
| Koppen-Geiger raster | `data/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif` |

The temperature files are repository-bundled application data rather than
package datasets loaded with `data()`. The NOAA records focus on California
stations for years 2011 through 2020. Raw CSV files use station identifiers,
station names, latitude, longitude, elevation, date, and temperature columns;
processed `.Rda` files store year-specific station-by-day temperature matrices
used by the model pipelines.

The backend validates these paths before loading year-specific data. Generated
outputs are written to `results/` by default.

## References

- Epskamp, S. (2017). graphicalVAR: Graphical VAR for experience sampling
  data. CRAN package documentation.
  https://cran.r-project.org/web/packages/graphicalVAR/graphicalVAR.pdf
- Isufi, E., Loukas, A., Perraudin, N., and Leus, G. (2019). Forecasting Time
  Series With VARMA Recursions on Graphs. IEEE Transactions on Signal
  Processing, vol. 67, no. 18, pp. 4870-4885, 15 Sept. 2019.
  https://doi.org/10.1109/TSP.2019.2929930
