#################################################
## TEMP APPLICATION: BACKEND AND CONFIGURATION ##
#################################################

## This file is the central backend for the temperature application.
## It keeps paths, year-specific file choices, defaults, summaries, and
## calls into the TAR-GAR/G-VAR and graphicalVAR pipelines in one place.

.temp_application_backend_file <- tryCatch(
  normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE),
  error = function(e) NA_character_
)

.temp_application_dir <- if (!is.na(.temp_application_backend_file)) {
  backend_dir <- dirname(.temp_application_backend_file)
  if (basename(backend_dir) == "R") dirname(backend_dir) else backend_dir
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

temp_app_dir <- function() {
  .temp_application_dir
}

temp_app_path <- function(...) {
  file.path(temp_app_dir(), ...)
}

temp_code_path <- function(...) {
  temp_app_path("R", ...)
}

temp_data_path <- function(...) {
  temp_app_path("data", ...)
}

temp_temperature_path <- function(...) {
  temp_data_path("temperature", ...)
}

temp_station_path <- function(...) {
  temp_data_path("stations", ...)
}

with_temp_app_dir <- function(expr) {
  old_wd <- getwd()
  setwd(temp_app_dir())
  on.exit(setwd(old_wd), add = TRUE)
  force(expr)
}

temp_source_once <- function(script, marker, force = FALSE) {
  if (!force && exists(marker, mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }

  source(temp_code_path(script), local = .GlobalEnv, chdir = TRUE)
  invisible(TRUE)
}

ensure_sgm_package <- function() {
  if (!requireNamespace("SGM", quietly = TRUE)) {
    stop(
      "The SGM package is required. Install it from the local package directory before running this application.",
      call. = FALSE
    )
  }

  required_exports <- c("fit_TAR_GAR", "TARGAR_pred", "GVAR_fit")
  missing_exports <- setdiff(required_exports, getNamespaceExports("SGM"))
  missing_internals <- c()
  for (fn in c("LogLike", "BIC")) {
    if (!exists(fn, envir = asNamespace("SGM"), mode = "function")) {
      missing_internals <- c(missing_internals, fn)
    }
  }

  if (length(missing_exports) > 0 || length(missing_internals) > 0) {
    stop(
      "The installed SGM package is missing required temperature-application functions: ",
      paste(c(missing_exports, missing_internals), collapse = ", "),
      ". Reinstall SGM from this repository before running the temperature application.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

ensure_temp_targar_backend <- function(force = FALSE) {
  ensure_sgm_package()
  temp_source_once("temp_pipeline.R", "TARGAR_on_temp", force = force)
}

ensure_temp_graphicalvar_backend <- function(force = FALSE) {
  temp_source_once("graphicalvar_temp_pipeline.R", "graphicalVAR_lag_sweep_on_temp", force = force)
}

temp_available_years <- function() {
  files <- basename(c(
    Sys.glob(temp_temperature_path("temp_*.Rda")),
    Sys.glob(temp_temperature_path("temp_*.csv"))
  ))
  files <- files[grepl("^temp_[0-9]{4}\\.(Rda|csv)$", files)]
  sort(unique(as.integer(sub("^temp_([0-9]{4})\\..*$", "\\1", files))))
}

temp_validate_year <- function(year) {
  year <- as.integer(year)
  years <- temp_available_years()

  if (length(year) != 1 || is.na(year) || !(year %in% years)) {
    stop(
      "Unsupported year: ", paste(year, collapse = ", "),
      ". Available years are: ", paste(years, collapse = ", "),
      call. = FALSE
    )
  }

  year
}

temp_year_config <- function(year, prefer_rda = TRUE) {
  year <- temp_validate_year(year)

  rda_path <- temp_temperature_path(sprintf("temp_%s.Rda", year))
  csv_path <- temp_temperature_path(sprintf("temp_%s.csv", year))
  has_rda <- file.exists(rda_path)
  has_csv <- file.exists(csv_path)

  if (prefer_rda && has_rda) {
    data_path <- rda_path
    run_prep <- FALSE
  } else if (has_csv) {
    data_path <- csv_path
    run_prep <- TRUE
  } else if (has_rda) {
    data_path <- rda_path
    run_prep <- FALSE
  } else {
    stop("No temperature data file found for ", year, call. = FALSE)
  }

  station_data_path <- if (year >= 2011 && year <= 2015) {
    temp_temperature_path("temp_CA_2010_2015.csv")
  } else if (has_csv) {
    csv_path
  } else {
    NA_character_
  }

  short_year <- substr(as.character(year), 3, 4)
  knn_dir <- if (year == 2020) {
    temp_data_path("knn_adj_matrices")
  } else {
    temp_data_path("knn_adj_matrices", paste0("tr", year))
  }

  latlong_candidates <- c(
    temp_station_path(sprintf("latlong_%s.csv", short_year)),
    file.path(knn_dir, sprintf("latlong_%s.csv", short_year))
  )
  latlong_path <- latlong_candidates[file.exists(latlong_candidates)][1]
  if (is.na(latlong_path)) {
    latlong_path <- latlong_candidates[1]
  }

  list(
    year = year,
    path_train = data_path,
    path_test = data_path,
    data_path = data_path,
    data_type = if (run_prep) "csv" else "rda",
    run_prep = run_prep,
    station_data_path = station_data_path,
    latlong_path = latlong_path,
    knn_dir = knn_dir,
    koppen_raster_path = temp_koppen_raster_path()
  )
}

temp_knn_matrix_path <- function(year, k) {
  year <- as.integer(year)
  base_dir <- if (identical(year, 2020L)) {
    temp_data_path("knn_adj_matrices")
  } else {
    temp_data_path("knn_adj_matrices", paste0("tr", year))
  }
  path <- file.path(base_dir, paste0("adjacency_matrix_k_", k, ".mtx"))

  if (!file.exists(path)) {
    stop("Missing kNN adjacency matrix: ", path, call. = FALSE)
  }

  path
}

temp_load_knn_adjacencies <- function(year, k = c(3, 7, 10), as_matrix = TRUE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required to load kNN adjacency matrices.", call. = FALSE)
  }

  out <- lapply(k, function(k_current) {
    A <- Matrix::readMM(temp_knn_matrix_path(year, k_current))
    if (as_matrix) {
      A <- as.matrix(A)
    }
    A
  })
  names(out) <- paste0(k, "NN")
  out
}

temp_koppen_raster_path <- function(period = "1991_2020",
                                    resolution = "0p00833333") {
  path <- temp_data_path(
    "koppen_geiger_tif",
    period,
    paste0("koppen_geiger_", resolution, ".tif")
  )

  if (!file.exists(path)) {
    stop("Missing Koppen-Geiger raster: ", path, call. = FALSE)
  }

  path
}

temp_results_dir <- function() {
  temp_app_path("results")
}

temp_result_file <- function(model = c("targar", "graphicalvar"),
                             year,
                             output_dir = temp_results_dir()) {
  model <- match.arg(model)
  filename <- switch(
    model,
    targar = paste0("TARGAR_list_tr_", year, "_t_", year, ".RData"),
    graphicalvar = paste0("graphicalVAR_list_tr_", year, "_t_", year, ".RData")
  )
  file.path(output_dir, filename)
}

temp_targar_defaults <- function() {
  n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)

  list(
    C.v = c(8, 4, 1.5, 1, 0.5, 0.25, 0.1),
    C.thre = exp(seq(log(0.75), log(0.075), length.out = 12)),
    plot_temp = TRUE,
    heatmap = TRUE,
    n.cores = min(9, n_cores),
    verbose = TRUE,
    run_prep = NULL,
    num.cluster = 7,
    standardize = FALSE,
    use_base_model = FALSE,
    forecast = FALSE,
    gvar = TRUE,
    stationary = TRUE
  )
}

temp_graphicalvar_defaults <- function() {
  list(
    plot_temp = TRUE,
    run_prep = NULL,
    num.cluster = 7,
    standardize = FALSE,
    forecast = FALSE,
    lag_specs = list("lag1" = 1, "lag2" = 1:2, "lag3" = 1:3),
    nLambda = 20,
    lambda_beta = c(10, 8, 4, 1, 0.5),
    lambda_kappa = exp(seq(log(10), log(0.05), length.out = 10)),
    gamma = 0.5,
    centerWithin = FALSE,
    penalize.diagonal = FALSE,
    net_size_tar = NULL,
    refit = FALSE,
    use_refit_forecast = FALSE,
    graph_type = "contemporaneous",
    include_knn = TRUE,
    save_plots = TRUE,
    n.cores.lag = 3,
    n.cores.beta = 2,
    beta_block_size = 1,
    verbose = TRUE
  )
}

temp_merge_args <- function(defaults, overrides = list()) {
  if (is.null(overrides)) {
    return(defaults)
  }
  if (!is.list(overrides)) {
    stop("Argument overrides must be provided as a list.", call. = FALSE)
  }
  utils::modifyList(defaults, overrides)
}

temp_selected_targar_net_size <- function(results_tar) {
  ebic_model <- apply(results_tar$ebic, c(1, 2), min, na.rm = TRUE)
  selected <- which(ebic_model == min(ebic_model, na.rm = TRUE), arr.ind = TRUE)
  if (nrow(selected) == 0) {
    return(NULL)
  }

  as.integer(results_tar$netsize[selected[1, 1], selected[1, 2]])
}

run_temp_targar <- function(year,
                            save_results = TRUE,
                            output_dir = temp_results_dir(),
                            args = list()) {
  config <- temp_year_config(year)
  defaults <- temp_targar_defaults()
  defaults$run_prep <- config$run_prep

  call_args <- c(
    list(
      path_train = config$path_train,
      path_test = config$path_test,
      year = config$year,
      output_dir = output_dir
    ),
    defaults
  )
  call_args <- temp_merge_args(call_args, args)

  ensure_temp_targar_backend()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results_tar <- with_temp_app_dir({
    do.call(TARGAR_on_temp, call_args)
  })

  if (isTRUE(save_results)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    save(results_tar, file = temp_result_file("targar", config$year, output_dir))
  }

  results_tar
}

run_temp_graphicalvar <- function(year,
                                  save_results = TRUE,
                                  output_dir = temp_results_dir(),
                                  net_size_tar = NULL,
                                  args = list()) {
  config <- temp_year_config(year)
  defaults <- temp_graphicalvar_defaults()
  defaults$run_prep <- config$run_prep
  if (!is.null(net_size_tar)) {
    defaults$net_size_tar <- net_size_tar
  }

  call_args <- c(
    list(
      path_train = config$path_train,
      path_test = config$path_test,
      year = config$year,
      output_dir = output_dir
    ),
    defaults
  )
  call_args <- temp_merge_args(call_args, args)

  ensure_temp_graphicalvar_backend()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results_gvar <- with_temp_app_dir({
    do.call(graphicalVAR_lag_sweep_on_temp, call_args)
  })

  if (isTRUE(save_results)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    save(results_gvar, file = temp_result_file("graphicalvar", config$year, output_dir))
  }

  results_gvar
}

run_temp_analysis <- function(year,
                              model = c("both", "targar", "graphicalvar"),
                              save_results = TRUE,
                              output_dir = temp_results_dir(),
                              targar_args = list(),
                              graphicalvar_args = list(),
                              net_size_tar = NULL,
                              print_summary = TRUE) {
  model <- match.arg(model)
  config <- temp_year_config(year)
  results <- list(year = config$year, config = config)

  message("Temperature application year: ", config$year)
  message("Training file: ", config$path_train)
  message("Data mode: ", config$data_type, " (run_prep = ", config$run_prep, ")")

  if (model %in% c("both", "targar")) {
    results$targar <- run_temp_targar(
      year = config$year,
      save_results = save_results,
      output_dir = output_dir,
      args = targar_args
    )

    if (isTRUE(print_summary)) {
      print_targar_summary(results$targar)
    }
  }

  if (model %in% c("both", "graphicalvar")) {
    if (is.null(net_size_tar) && !is.null(results$targar)) {
      net_size_tar <- temp_selected_targar_net_size(results$targar)
      message("Using selected TAR-GAR network size for graphicalVAR matching: ", net_size_tar)
    }

    results$graphicalvar <- run_temp_graphicalvar(
      year = config$year,
      save_results = save_results,
      output_dir = output_dir,
      net_size_tar = net_size_tar,
      args = graphicalvar_args
    )

    if (isTRUE(print_summary)) {
      print_graphicalvar_summary(results$graphicalvar)
    }
  }

  invisible(results)
}

print_targar_summary <- function(results_tar) {
  print("Network Sizes:")
  print(results_tar$netsize)

  print("Overlap Indices:")
  print(results_tar$overlap)

  print("eBIC")
  print(apply(results_tar$ebic, c(1, 2), min))

  print("Convergence:")
  print(results_tar$conv.prop.mat)

  print("Forecasts:")
  print(results_tar$forecasts)

  print("Purity")
  print(results_tar$purity)

  invisible(results_tar)
}

print_graphicalvar_summary <- function(results_gvar) {
  print("Network Size:")
  print(sum(results_gvar$adjacency) / 2)

  print("eBIC by Lag")
  print(results_gvar$ebic_by_lag)

  print("eBIC by Lag and Beta Block")
  print(results_gvar$ebic_by_block)

  print("Selected Lag")
  print(results_gvar$selected_lag)

  print("Selected eBIC")
  print(results_gvar$fit$EBIC)

  print("Tuning Path:")
  print(results_gvar$fit$path)

  print("Forecasts:")
  print(results_gvar$forecasts)

  print("Purity")
  print(results_gvar$purity)

  if (!is.null(results_gvar$graphicalVAR_sim_size)) {
    print("Size-matched graphicalVAR Selection:")
    print(results_gvar$graphicalVAR_sim_size$selection)

    print("Size-matched graphicalVAR Network Size:")
    print(sum(results_gvar$graphicalVAR_sim_size$adjacency) / 2)

    print("Size-matched graphicalVAR Forecasts:")
    print(results_gvar$graphicalVAR_sim_size$forecasts)

    print("Size-matched graphicalVAR Purity:")
    print(results_gvar$graphicalVAR_sim_size$purity)
  }

  invisible(results_gvar)
}

temp_parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE),
                                default_year = 2011,
                                default_model = "both",
                                default_save_results = TRUE) {
  out <- list(
    year = as.integer(default_year),
    model = default_model,
    save_results = default_save_results,
    help = FALSE
  )

  if (length(args) == 0) {
    return(out)
  }

  for (arg in args) {
    if (arg %in% c("-h", "--help")) {
      out$help <- TRUE
    } else if (arg == "--no-save") {
      out$save_results <- FALSE
    } else if (arg == "--save") {
      out$save_results <- TRUE
    } else if (grepl("^--year=", arg)) {
      out$year <- as.integer(sub("^--year=", "", arg))
    } else if (grepl("^--model=", arg)) {
      out$model <- sub("^--model=", "", arg)
    } else if (grepl("^[0-9]{4}$", arg)) {
      out$year <- as.integer(arg)
    } else if (arg %in% c("both", "targar", "graphicalvar")) {
      out$model <- arg
    } else {
      warning("Ignoring unrecognized argument: ", arg, call. = FALSE)
    }
  }

  out
}

temp_print_cli_help <- function() {
  cat(
    "Usage:\n",
    "  Rscript run_temp_application.R --year=2011 --model=both\n",
    "\n",
    "Models:\n",
    "  both, targar, graphicalvar\n",
    "\n",
    "Options:\n",
    "  --year=YYYY     Temperature year to analyze. Available: ",
    paste(temp_available_years(), collapse = ", "),
    "\n",
    "  --model=MODEL   Model family to run.\n",
    "  --no-save       Run without saving the result .RData file.\n",
    sep = ""
  )
}
