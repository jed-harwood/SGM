############################################
## TEMP APPLICATION: GRAPHICALVAR PIPELINE ##
############################################

## This script intentionally reuses the temp_detrend function from
## temp_pipeline.R without evaluating the rest of that TAR-GAR pipeline.

library(dplyr)
library(readr)
library(ggplot2)
library(zoo)
library(tidyr)
library(stringr)
library(forcats)
library(igraph)
library(raster)
library(RSpectra)
library(Matrix)
library(sf)
library(graphicalVAR)
library(foreach)
library(doParallel)


load_temp_detrend <- function(path = "temp_pipeline.R") {
  if (!file.exists(path) && exists("temp_code_path", mode = "function")) {
    path <- temp_code_path(path)
  } else if (!file.exists(path) && exists("temp_app_path", mode = "function")) {
    path <- temp_app_path(path)
  }

  exprs <- parse(path)
  is_temp_detrend <- vapply(
    exprs,
    function(expr) {
      is.call(expr) &&
        identical(expr[[1]], as.name("=")) &&
        identical(expr[[2]], as.name("temp_detrend"))
    },
    logical(1)
  )
  idx <- which(is_temp_detrend)
  if (length(idx) != 1) {
    stop("Could not uniquely locate temp_detrend in ", path, call. = FALSE)
  }
  eval(exprs[[idx]], envir = .GlobalEnv)
  invisible(temp_detrend)
}

load_temp_detrend()


make_graphicalvar_regularize_mat_beta <- function(p, lags, penalize.diagonal = TRUE) {
  regularize_mat_beta <- matrix(TRUE, nrow = p, ncol = p * length(lags))

  if (!penalize.diagonal) {
    for (lag_idx in seq_along(lags)) {
      cols <- ((lag_idx - 1) * p + 1):(lag_idx * p)
      regularize_mat_beta[cbind(seq_len(p), cols)] <- FALSE
    }
  }

  regularize_mat_beta
}


make_graphicalvar_lambda_grid <- function(train,
                                          lags = 1:3,
                                          nLambda = 20,
                                          scale = FALSE,
                                          centerWithin = FALSE,
                                          lambda_min_beta = 0.05,
                                          lambda_min_kappa = 0.05,
                                          penalize.diagonal = TRUE,
                                          lambda_beta = NULL,
                                          lambda_kappa = NULL) {
  if (!is.null(lambda_beta) && !is.null(lambda_kappa)) {
    return(list(lambda_beta = lambda_beta, lambda_kappa = lambda_kappa))
  }

  ts_data <- graphicalVAR:::tsData(
    as.data.frame(train),
    scale = scale,
    centerWithin = centerWithin,
    lags = lags
  )
  lams <- graphicalVAR:::generate_lambdas(
    as.matrix(ts_data$data_l),
    as.matrix(ts_data$data_c),
    nLambda,
    nLambda,
    lambda_min_kappa = lambda_min_kappa,
    lambda_min_beta = lambda_min_beta,
    penalize.diagonal = penalize.diagonal
  )

  if (!is.null(lambda_beta)) {
    lams$lambda_beta <- lambda_beta
  }
  if (!is.null(lambda_kappa)) {
    lams$lambda_kappa <- lambda_kappa
  }

  lams
}


fit_graphicalvar_temp <- function(train,
                                  lags = 1:3,
                                  nLambda = 20,
                                  gamma = 0.5,
                                  scale = FALSE,
                                  centerWithin = FALSE,
                                  lambda_beta = NULL,
                                  lambda_kappa = NULL,
                                  lambda_min_beta = 0.05,
                                  lambda_min_kappa = 0.05,
                                  regularize_mat_beta = NULL,
                                  regularize_mat_kappa = NULL,
                                  penalize.diagonal = TRUE,
                                  refit = FALSE,
                                  verbose = TRUE,
                                  ...) {
  train <- as.matrix(train)
  storage.mode(train) <- "double"
  p <- ncol(train)

  if (is.null(regularize_mat_beta)) {
    regularize_mat_beta <- make_graphicalvar_regularize_mat_beta(
      p = p,
      lags = lags,
      penalize.diagonal = penalize.diagonal
    )
  }

  gvar_args <- list(
    data = train,
    lags = lags,
    nLambda = nLambda,
    gamma = gamma,
    scale = scale,
    centerWithin = centerWithin,
    lambda_min_beta = lambda_min_beta,
    lambda_min_kappa = lambda_min_kappa,
    regularize_mat_beta = regularize_mat_beta,
    penalize.diagonal = penalize.diagonal,
    verbose = verbose
  )
  if (!is.null(regularize_mat_kappa)) {
    gvar_args$regularize_mat_kappa <- regularize_mat_kappa
  }
  if (!is.null(lambda_beta)) {
    gvar_args$lambda_beta <- lambda_beta
  }
  if (!is.null(lambda_kappa)) {
    gvar_args$lambda_kappa <- lambda_kappa
  }
  gvar_args <- c(gvar_args, list(...))

  fit <- do.call(graphicalVAR::graphicalVAR, gvar_args)

  refit_fit <- NULL
  if (refit) {
    refit_fit <- try(refit_graphicalvar_beta(fit, train, lags = lags), silent = TRUE)
    if (inherits(refit_fit, "try-error")) {
      warning("Refitting non-zero temporal entries failed; returning penalized fit only.")
      refit_fit <- NULL
    }
  }

  list(
    fit = fit,
    refit = refit_fit,
    lags = lags,
    tuning = list(
      nLambda = nLambda,
      gamma = gamma,
      centerWithin = centerWithin,
      lambda_beta = lambda_beta,
      lambda_kappa = lambda_kappa,
      lambda_min_beta = lambda_min_beta,
      lambda_min_kappa = lambda_min_kappa,
      regularize_mat_beta = regularize_mat_beta,
      regularize_mat_kappa = regularize_mat_kappa,
      penalize.diagonal = penalize.diagonal
    )
  )
}


split_graphicalvar_beta <- function(beta, p, lags) {
  lag_mats <- vector("list", length(lags))
  names(lag_mats) <- paste0("lag", lags)

  for (i in seq_along(lags)) {
    cols <- 1 + ((i - 1) * p + 1):(i * p)
    lag_mats[[i]] <- as.matrix(beta[, cols, drop = FALSE])
    rownames(lag_mats[[i]]) <- rownames(beta)
    colnames(lag_mats[[i]]) <- rownames(beta)
  }

  lag_mats
}


make_graphicalvar_design <- function(data, lags) {
  data <- as.matrix(data)
  max_lag <- max(lags)
  n <- nrow(data)
  p <- ncol(data)

  Y <- data[(max_lag + 1):n, , drop = FALSE]
  X <- matrix(1, nrow = n - max_lag, ncol = 1 + p * length(lags))
  colnames(X) <- c(
    "intercept",
    unlist(lapply(lags, function(lag) paste0(colnames(data), "_lag", lag)))
  )

  for (i in seq_along(lags)) {
    lag <- lags[i]
    cols <- 1 + ((i - 1) * p + 1):(i * p)
    X[, cols] <- data[(max_lag + 1 - lag):(n - lag), , drop = FALSE]
  }

  list(X = X, Y = Y)
}


refit_graphicalvar_beta <- function(fit, train, lags = fit$data$lags) {
  train <- as.matrix(train)
  colnames(train) <- fit$labels
  p <- ncol(train)
  design <- make_graphicalvar_design(train, lags = lags)
  beta_pen <- as.matrix(fit$beta)
  beta_refit <- matrix(0, nrow = p, ncol = ncol(beta_pen), dimnames = dimnames(beta_pen))
  beta_refit[, 1] <- beta_pen[, 1]

  for (j in seq_len(p)) {
    support <- which(beta_pen[j, -1] != 0) + 1
    cols <- c(1, support)
    if (length(cols) == 1) {
      beta_refit[j, 1] <- mean(design$Y[, j])
    } else {
      coef_j <- stats::lm.fit(x = design$X[, cols, drop = FALSE], y = design$Y[, j])$coefficients
      coef_j[is.na(coef_j)] <- 0
      beta_refit[j, cols] <- coef_j
    }
  }

  resid <- design$Y - design$X %*% t(beta_refit)
  sigma <- stats::cov(resid)
  kappa_refit <- try(solve(sigma), silent = TRUE)
  if (inherits(kappa_refit, "try-error")) {
    kappa_refit <- MASS::ginv(sigma)
  }
  dimnames(kappa_refit) <- list(fit$labels, fit$labels)

  list(
    beta = beta_refit,
    beta_lags = split_graphicalvar_beta(beta_refit, p = p, lags = lags),
    residuals = resid,
    sigma = sigma,
    kappa = kappa_refit,
    support = beta_pen != 0
  )
}


graphicalvar_adjacency <- function(model, use_refit = FALSE, type = c("contemporaneous", "temporal")) {
  type <- match.arg(type)
  fit <- model$fit
  p <- length(fit$labels)

  if (type == "contemporaneous") {
    A <- 1 * (fit$PCC != 0)
    diag(A) <- 0
    A <- 1 * ((A + t(A)) > 0)
    dimnames(A) <- list(fit$labels, fit$labels)
    return(A)
  }

  beta <- fit$beta
  if (use_refit && !is.null(model$refit)) {
    beta <- model$refit$beta
  }
  lag_mats <- split_graphicalvar_beta(beta, p = p, lags = model$lags)
  A_dir <- Reduce("+", lapply(lag_mats, function(B) 1 * (B != 0)))
  diag(A_dir) <- 0
  A <- 1 * ((A_dir + t(A_dir)) > 0)
  dimnames(A) <- list(fit$labels, fit$labels)
  A
}


graphicalvar_net_size <- function(model, type = c("contemporaneous", "temporal")) {
  A <- graphicalvar_adjacency(model, use_refit = FALSE, type = type)
  sum(A[upper.tri(A)] != 0)
}


extract_graphicalvar_candidate <- function(block_model, row_index) {
  fit <- block_model$fit
  candidate <- fit$allResults[[row_index]]
  fit_candidate <- fit

  fit_candidate$beta <- candidate$beta
  fit_candidate$kappa <- candidate$kappa
  fit_candidate$EBIC <- candidate$EBIC

  colnames(fit_candidate$beta) <- colnames(fit$beta)
  rownames(fit_candidate$beta) <- rownames(fit$beta)
  colnames(fit_candidate$kappa) <- colnames(fit$kappa)
  rownames(fit_candidate$kappa) <- rownames(fit$kappa)

  fit_candidate$PCC <- graphicalVAR:::computePCC(fit_candidate$kappa)
  colnames(fit_candidate$PCC) <- rownames(fit_candidate$PCC) <- fit$labels

  fit_candidate$PDC <- fit$PDC
  if (1 %in% block_model$lags) {
    pdc_candidate <- try(
      graphicalVAR:::computePDC(
        fit_candidate$beta[, c("1", paste0(fit$labels, "_lag1")), drop = FALSE],
        fit_candidate$kappa
      ),
      silent = TRUE
    )
    if (!inherits(pdc_candidate, "try-error")) {
      fit_candidate$PDC <- pdc_candidate
      colnames(fit_candidate$PDC) <- rownames(fit_candidate$PDC) <- fit$labels
    }
  }

  out <- block_model
  out$fit <- fit_candidate
  out$refit <- NULL
  out
}


select_graphicalvar_size_match <- function(beta_block_results,
                                           selected_lag_name,
                                           net_size_tar,
                                           graph_type = c("contemporaneous", "temporal")) {
  graph_type <- match.arg(graph_type)

  if (is.null(net_size_tar)) {
    return(NULL)
  }
  if (length(net_size_tar) != 1 || is.na(net_size_tar) || net_size_tar < 0) {
    stop("net_size_tar must be a single non-negative integer.", call. = FALSE)
  }
  if (net_size_tar != as.integer(net_size_tar)) {
    stop("net_size_tar must be an integer.", call. = FALSE)
  }
  net_size_tar <- as.integer(net_size_tar)

  selected_lag_blocks <- beta_block_results[vapply(
    beta_block_results,
    function(res) res$lag_name == selected_lag_name,
    logical(1)
  )]

  candidates <- list()
  candidate_index <- 1
  for (block_id in seq_along(selected_lag_blocks)) {
    block_result <- selected_lag_blocks[[block_id]]
    path_current <- block_result$model$fit$path

    for (row_index in seq_len(nrow(path_current))) {
      candidate_model <- try(
        extract_graphicalvar_candidate(block_result$model, row_index),
        silent = TRUE
      )
      if (inherits(candidate_model, "try-error")) {
        next
      }

      size_current <- try(
        graphicalvar_net_size(candidate_model, type = graph_type),
        silent = TRUE
      )
      if (inherits(size_current, "try-error") || is.na(size_current)) {
        next
      }

      candidates[[candidate_index]] <- data.frame(
        lag_name = block_result$lag_name,
        lags = paste(block_result$model$lags, collapse = ","),
        block_id = block_result$block_id,
        row_index = row_index,
        EBIC = path_current$ebic[row_index],
        lambda_beta = path_current$beta[row_index],
        lambda_kappa = path_current$kappa[row_index],
        net_size = size_current,
        net_size_tar = net_size_tar,
        abs_diff = abs(size_current - net_size_tar),
        stringsAsFactors = FALSE
      )
      candidate_index <- candidate_index + 1
    }
  }

  if (length(candidates) == 0) {
    warning("No valid graphicalVAR candidates were available for net_size_tar selection.")
    return(NULL)
  }

  candidate_table <- do.call(rbind, candidates)
  candidate_table <- candidate_table[is.finite(candidate_table$EBIC), , drop = FALSE]
  if (nrow(candidate_table) == 0) {
    warning("No finite-EBIC graphicalVAR candidates were available for net_size_tar selection.")
    return(NULL)
  }

  candidate_table <- candidate_table[order(candidate_table$abs_diff, candidate_table$EBIC), , drop = FALSE]
  best <- candidate_table[1, , drop = FALSE]
  best_block <- selected_lag_blocks[[which(vapply(
    selected_lag_blocks,
    function(res) res$block_id == best$block_id,
    logical(1)
  ))]]
  best_model <- extract_graphicalvar_candidate(best_block$model, best$row_index)

  list(
    model = best_model,
    selection = best,
    candidates = candidate_table
  )
}


predict_graphicalvar_one <- function(history, beta, lags) {
  history <- as.matrix(history)
  p <- ncol(history)
  x <- c(1, unlist(lapply(lags, function(lag) history[nrow(history) + 1 - lag, ])))
  as.numeric(beta %*% x)
}


graphicalvar_pred <- function(model, testdata, n.ahead = 1, use_refit = FALSE) {
  testdata <- as.matrix(testdata)
  colnames(testdata) <- model$fit$labels
  lags <- model$lags
  max_lag <- max(lags)

  beta <- model$fit$beta
  if (use_refit && !is.null(model$refit)) {
    beta <- model$refit$beta
  }

  n <- nrow(testdata)
  p <- ncol(testdata)
  num_forecasts <- n - max_lag - n.ahead + 1
  pred <- matrix(NA, nrow = num_forecasts, ncol = p, dimnames = list(NULL, colnames(testdata)))
  truth <- matrix(NA, nrow = num_forecasts, ncol = p, dimnames = list(NULL, colnames(testdata)))

  for (i in seq_len(num_forecasts)) {
    hist <- testdata[i:(i + max_lag - 1), , drop = FALSE]
    for (h in seq_len(n.ahead)) {
      next_pred <- predict_graphicalvar_one(hist, beta = beta, lags = lags)
      hist <- rbind(hist, next_pred)
    }
    pred[i, ] <- next_pred
    truth[i, ] <- testdata[i + max_lag + n.ahead - 1, ]
  }

  err <- truth - pred
  rmse_node <- sqrt(colMeans(err^2, na.rm = TRUE))
  rmae_node <- colMeans(abs(err), na.rm = TRUE)
  node_corr <- suppressWarnings(vapply(
    seq_len(p),
    function(j) stats::cor(truth[, j], pred[, j], use = "pairwise.complete.obs"),
    numeric(1)
  ))

  list(
    pred = pred,
    truth = truth,
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    rmse.node = summary(rmse_node),
    rmae.node = summary(rmae_node),
    node.corr = summary(node_corr)
  )
}


make_forecast_graphicalvar <- function(model, train, test, use_refit = FALSE) {
  list(
    "1-step" = graphicalvar_pred(model, test, n.ahead = 1, use_refit = use_refit)$rmse.node[4],
    "2-step" = graphicalvar_pred(model, test, n.ahead = 2, use_refit = use_refit)$rmse.node[4],
    "3-step" = graphicalvar_pred(model, test, n.ahead = 3, use_refit = use_refit)$rmse.node[4],
    "1-step.train" = graphicalvar_pred(model, train, n.ahead = 1, use_refit = use_refit)$rmse.node[4]
  )
}


spectral_cluster <- function(A, K, normalized = TRUE, nstart = 10, seed = 1) {
  stopifnot(nrow(A) == ncol(A))
  A <- as.matrix(A)
  diag(A) <- 0
  A <- 1 * ((A + t(A)) > 0)

  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  L <- laplacian_matrix(g, normalized = normalized, sparse = TRUE)

  if (requireNamespace("RSpectra", quietly = TRUE)) {
    es <- RSpectra::eigs_sym(L, k = K, which = "SM")
    U <- es$vectors
    evals <- es$values
  } else {
    E <- eigen(as.matrix(L), symmetric = TRUE)
    idx <- order(E$values, decreasing = FALSE)[1:K]
    U <- E$vectors[, idx, drop = FALSE]
    evals <- E$values[idx]
  }

  if (normalized) {
    rn <- sqrt(rowSums(U^2))
    rn[rn == 0] <- 1
    Tm <- U / rn
  } else {
    Tm <- U
  }

  set.seed(seed)
  km <- stats::kmeans(Tm, centers = K, nstart = nstart)
  mem <- km$cluster
  cl <- make_clusters(g, membership = mem)

  list(
    membership = mem,
    communities = cl,
    centers = km$centers,
    eigvecs = U,
    eigvals = evals,
    graph = g
  )
}


read_temp_station_data <- function(path_train, year = 2020) {
  if (str_detect(path_train, ".csv")) {
    return(read_csv(path_train, col_types = list("f", "f", "d", "d", "d", "D", "d", "d")))
  }

  if (year >= 2011 && year <= 2015) {
    path_station <- if (exists("temp_temperature_path", mode = "function")) {
      temp_temperature_path("temp_CA_2010_2015.csv")
    } else {
      "data/temperature/temp_CA_2010_2015.csv"
    }
    return(read_csv(path_station, col_types = list("f", "f", "d", "d", "d", "D", "d", "d")))
  }

  path_train_fix <- str_replace(path_train, "Rda", "csv")
  read_csv(path_train_fix, col_types = list("f", "f", "d", "d", "d", "D", "d", "d"))
}


get_koppen_labels <- function(path_train, train, year = 2020) {
  temp <- read_temp_station_data(path_train, year = year)

  koppen_path <- if (exists("temp_koppen_raster_path", mode = "function")) {
    temp_koppen_raster_path()
  } else {
    "data/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"
  }
  koppen_raster <- raster(koppen_path)

  latlong <- temp %>%
    dplyr::filter(STATION %in% colnames(train)) %>%
    dplyr::select(STATION, LATITUDE, LONGITUDE) %>%
    distinct() %>%
    dplyr::rename(name = STATION, x = LONGITUDE, y = LATITUDE)

  stations_sf <- st_as_sf(latlong, coords = c("x", "y"), crs = 4326)
  latlong$koppen_code <- raster::extract(koppen_raster, stations_sf)
  latlong$koppen_code <- fct_recode(as.factor(latlong$koppen_code), BWh = "4", BWk = "5", BSh = "6", BSk = "7", Csa = "8", Csb = "9", Dsb = "18")
  latlong$koppen_code <- fct_na_value_to_level(latlong$koppen_code, "uncategorized")

  stations_manual <- list(
    "99401899999" = "Csa",
    "72493723289" = "Csb",
    "99402899999" = "Csb",
    "99404199999" = "Csb",
    "72494023234" = "Csb",
    "72594024213" = "Csb",
    "72292800369" = "BSk"
  )
  for (i in seq_len(nrow(latlong))) {
    if (latlong$name[i] %in% names(stations_manual)) {
      latlong$koppen_code[i] <- stations_manual[[as.character(latlong$name[i])]]
    }
  }

  latlong$koppen_code <- droplevels(latlong$koppen_code)
  latlong
}


load_knn_adjacencies <- function(year = 2020) {
  if (exists("temp_load_knn_adjacencies", mode = "function")) {
    return(temp_load_knn_adjacencies(year = year, k = c(3, 7, 10), as_matrix = TRUE))
  }

  if ((year == 2020)) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_10.mtx")
  } else if (year == 2019) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_10.mtx")
  } else if (year == 2018) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_10.mtx")
  } else if (year == 2017) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_10.mtx")
  } else if (year == 2016) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_10.mtx")
  } else if (year == 2015) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_10.mtx")
  } else if (year == 2014) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_10.mtx")
  } else if (year == 2013) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_10.mtx")
  } else if (year == 2012) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_10.mtx")
  } else if (year == 2011) {
    A.3nn <- Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_3.mtx")
    A.7nn <- Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_7.mtx")
    A.10nn <- Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_10.mtx")
  } else {
    stop("No kNN adjacency paths configured for year ", year, call. = FALSE)
  }

  list("3NN" = as.matrix(A.3nn), "7NN" = as.matrix(A.7nn), "10NN" = as.matrix(A.10nn))
}


cluster_purity <- function(membership, labels, K) {
  majority <- numeric(K)
  for (i in seq_len(K)) {
    tab <- table(labels[which(membership == i)])
    majority[i] <- if (length(tab) == 0) 0 else max(tab)
  }
  sum(majority) / length(membership)
}


plot_cluster_pdf <- function(spec, labels, title, file, layout = NULL) {
  if (is.null(layout)) {
    layout <- layout_with_fr(spec$graph)
  }
  lab_levels <- unique(labels)
  cols <- rainbow(length(lab_levels))
  names(cols) <- lab_levels
  vertex_cols <- unname(cols[as.character(labels)])
  freq <- table(labels)[lab_levels]

  pdf(file = file, width = 8.92, height = 4.86)
  plot(
    spec$communities,
    spec$graph,
    col = vertex_cols,
    cex = 0.5,
    vertex.label = "",
    layout = layout,
    vertex.size = 7,
    main = title
  )
  legend(
    "topright",
    legend = paste(lab_levels, paste0("(", freq, ")")),
    pch = 21,
    pt.bg = cols,
    title = "True Koppen Zone (Freq)",
    col = "black"
  )
  dev.off()
}


spec_cluster_graphicalvar <- function(graphicalvar_model,
                                      path_train,
                                      train,
                                      K = 7,
                                      year = 2020,
                                      include_knn = TRUE,
                                      use_refit = FALSE,
                                      graph_type = c("contemporaneous", "temporal"),
                                      save_plots = TRUE,
                                      output_prefix = "graphicalvar",
                                      output_dir = "results") {
  graph_type <- match.arg(graph_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (graph_type != "contemporaneous") {
    warning("Purity analysis is intended for the graphicalVAR contemporaneous undirected network.")
  }

  latlong <- get_koppen_labels(path_train = path_train, train = train, year = year)

  ### Now create the graph objects
  ##### graphicalVAR contemporaneous network
  A.gvar <- graphicalvar_adjacency(graphicalvar_model, use_refit = use_refit, type = graph_type)
  graph.gvar <- graph_from_adjacency_matrix(A.gvar,
                                            add.colnames = NULL,
                                            diag = FALSE,
                                            mode = "undirected")
  V(graph.gvar)$koppen <- latlong$koppen_code

  ###### kNN Graphs
  if (include_knn) {
    knn <- load_knn_adjacencies(year = year)
    A.3nn <- knn[["3NN"]]
    A.7nn <- knn[["7NN"]]
    A.10nn <- knn[["10NN"]]
  } else {
    A.3nn <- A.7nn <- A.10nn <- matrix(0, nrow = nrow(A.gvar), ncol = ncol(A.gvar))
  }

  A.3nn[A.3nn > 1e-6] <- 1
  A.7nn[A.7nn > 1e-6] <- 1
  A.10nn[A.10nn > 1e-6] <- 1

  graph_3nn <- graph_from_adjacency_matrix(A.3nn, mode = "undirected", diag = FALSE)
  graph_7nn <- graph_from_adjacency_matrix(A.7nn, mode = "undirected", diag = FALSE)
  graph_10nn <- graph_from_adjacency_matrix(A.10nn, mode = "undirected", diag = FALSE)

  V(graph_3nn)$koppen <- latlong$koppen_code
  V(graph_7nn)$koppen <- latlong$koppen_code
  V(graph_10nn)$koppen <- latlong$koppen_code

  ## Plot settings
  ##########################################
  ##### Layouts
  lay_3nn <- layout_with_fr(graph_3nn)
  lay_7nn <- layout_with_fr(graph_7nn)
  lay_10nn <- layout_with_fr(graph_10nn)
  lay_gvar <- layout_with_fr(graph.gvar)

  #### Colors and frequencies
  col_vec <- rainbow(8)
  vertex_col_list <- c()
  freq_vec <- c()
  j <- 1
  for (i in unique(V(graph.gvar)$koppen)) {
    vertex_col_list[[i]] <- col_vec[j]
    freq_vec[j] <- table(V(graph.gvar)$koppen)[i]
    j <- j + 1
  }

  col_vec <- c()
  for (i in V(graph.gvar)$koppen) {
    col_vec <- c(col_vec, vertex_col_list[[i]])
  }
  V(graph.gvar)$col <- col_vec
  ##########################################

  spec_3nn <- spectral_cluster(A = graph_3nn, K = K, normalized = TRUE, nstart = 25)
  spec_7nn <- spectral_cluster(A = graph_7nn, K = K, normalized = TRUE, nstart = 25)
  spec_10nn <- spectral_cluster(A = graph_10nn, K = K, normalized = TRUE, nstart = 25)
  spec_gvar <- spectral_cluster(A = graph.gvar, K = K, normalized = TRUE, nstart = 25)

  ## Calculate Purity and summary statistics from clusters
  majority_3nn <- c()
  majority_7nn <- c()
  majority_10nn <- c()
  majority_gvar <- c()

  majority2_3nn <- c()
  majority2_7nn <- c()
  majority2_10nn <- c()
  majority2_gvar <- c()

  cluster_info <- vector(mode = "list", length = K)
  for (i in 1:K) {
    ## Obtain indices of nodes in ith cluster
    idx_cluster_3nn <- which(spec_3nn$membership == i)
    idx_cluster_7nn <- which(spec_7nn$membership == i)
    idx_cluster_10nn <- which(spec_10nn$membership == i)
    idx_cluster_gvar <- which(spec_gvar$membership == i)

    ## Calculate length of cluster i for each graph
    cluster_3nn_i_len <- length(idx_cluster_3nn)
    cluster_7nn_i_len <- length(idx_cluster_7nn)
    cluster_10nn_i_len <- length(idx_cluster_10nn)
    cluster_gvar_i_len <- length(idx_cluster_gvar)

    tab_3nn <- sort(table(V(graph_3nn)$koppen[idx_cluster_3nn]), decreasing = TRUE)
    tab_7nn <- sort(table(V(graph_7nn)$koppen[idx_cluster_7nn]), decreasing = TRUE)
    tab_10nn <- sort(table(V(graph_10nn)$koppen[idx_cluster_10nn]), decreasing = TRUE)
    tab_gvar <- sort(table(V(graph.gvar)$koppen[idx_cluster_gvar]), decreasing = TRUE)

    ## Calculate most prevalent class in ith cluster
    majority_3nn[i] <- if (length(tab_3nn) >= 1) tab_3nn[1] else 0
    majority_7nn[i] <- if (length(tab_7nn) >= 1) tab_7nn[1] else 0
    majority_10nn[i] <- if (length(tab_10nn) >= 1) tab_10nn[1] else 0
    majority_gvar[i] <- if (length(tab_gvar) >= 1) tab_gvar[1] else 0

    ## Calculate the second most prevalent class in ith cluster
    majority2_3nn[i] <- if (length(tab_3nn) >= 2) tab_3nn[2] else NA
    majority2_7nn[i] <- if (length(tab_7nn) >= 2) tab_7nn[2] else NA
    majority2_10nn[i] <- if (length(tab_10nn) >= 2) tab_10nn[2] else NA
    majority2_gvar[i] <- if (length(tab_gvar) >= 2) tab_gvar[2] else NA

    ## Store Majority Size and Corresponding Label for Each cluster
    majority_3nn_zone <- if (length(tab_3nn) >= 1) names(tab_3nn[1]) else NA
    majority_7nn_zone <- if (length(tab_7nn) >= 1) names(tab_7nn[1]) else NA
    majority_10nn_zone <- if (length(tab_10nn) >= 1) names(tab_10nn[1]) else NA
    majority_gvar_zone <- if (length(tab_gvar) >= 1) names(tab_gvar[1]) else NA

    ### Second Majority Class
    majority_3nn_zone_2 <- if (length(tab_3nn) >= 2) names(tab_3nn[2]) else NA
    majority_7nn_zone_2 <- if (length(tab_7nn) >= 2) names(tab_7nn[2]) else NA
    majority_10nn_zone_2 <- if (length(tab_10nn) >= 2) names(tab_10nn[2]) else NA
    majority_gvar_zone_2 <- if (length(tab_gvar) >= 2) names(tab_gvar[2]) else NA

    ## Store summary statistics in nested list: cluster_info >> graph
    info_3nn <- list("size" = cluster_3nn_i_len, "majority size" = majority_3nn[i], "majority zone" = majority_3nn_zone, "majority zone 2" = majority_3nn_zone_2, "second count" = majority2_3nn[i], "membership" = idx_cluster_3nn, "graph" = graph_3nn)
    info_7nn <- list("size" = cluster_7nn_i_len, "majority size" = majority_7nn[i], "majority zone" = majority_7nn_zone, "majority zone 2" = majority_7nn_zone_2, "second count" = majority2_7nn[i], "membership" = idx_cluster_7nn, "graph" = graph_7nn)
    info_10nn <- list("size" = cluster_10nn_i_len, "majority size" = majority_10nn[i], "majority zone" = majority_10nn_zone, "majority zone 2" = majority_10nn_zone_2, "second count" = majority2_10nn[i], "membership" = idx_cluster_10nn, "graph" = graph_10nn)
    info_gvar <- list("size" = cluster_gvar_i_len, "majority size" = majority_gvar[i], "majority zone" = majority_gvar_zone, "majority zone 2" = majority_gvar_zone_2, "second count" = majority2_gvar[i], "membership" = idx_cluster_gvar, "graph" = graph.gvar)

    cluster_info[[i]] <- list("3nn" = info_3nn,
                              "7nn" = info_7nn,
                              "10nn" = info_10nn,
                              "graphicalVAR" = info_gvar)
  }

  purity_3nn <- sum(majority_3nn) / length(spec_3nn$membership)
  purity_7nn <- sum(majority_7nn) / length(spec_7nn$membership)
  purity_10nn <- sum(majority_10nn) / length(spec_10nn$membership)
  purity_gvar <- sum(majority_gvar) / length(spec_gvar$membership)

  save(cluster_info, file = file.path(output_dir, paste0(output_prefix, "_", graph_type, "_cluster_info_temp_", year, ".Rda")))
  print("Cluster information saved.")

  if (save_plots) {
    ## 3nn Plot
    pdf(file = file.path(output_dir, paste0("3nn_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
    plot(spec_3nn$communities, spec_3nn$graph,
         col = col_vec, cex = 0.5, vertex.label = "", layout = lay_3nn,
         vertex.size = 7, main = paste0("Clustering of 3NN Graph (Purity = ", round(purity_3nn, 3), ")"))
    legend("topright", legend = paste(unique(V(graph.gvar)$koppen), paste0("(", freq_vec, ")")),
           pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)")
    dev.off()

    ## 7nn Plot
    pdf(file = file.path(output_dir, paste0("7nn_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
    plot(spec_7nn$communities, spec_7nn$graph,
         col = col_vec, cex = 0.5, vertex.label = "", layout = lay_7nn,
         vertex.size = 7, main = paste0("Clustering of 7NN Graph (Purity = ", round(purity_7nn, 3), ")"))
    legend("topright", legend = paste(unique(V(graph.gvar)$koppen), paste0("(", freq_vec, ")")),
           pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)")
    dev.off()

    ## 10nn Plot
    pdf(file = file.path(output_dir, paste0("10nn_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
    plot(spec_10nn$communities, spec_10nn$graph,
         col = col_vec, cex = 0.5, vertex.label = "", layout = lay_10nn,
         vertex.size = 7, main = paste0("Clustering of 10NN Graph (Purity = ", round(purity_10nn, 3), ")"))
    legend("topright", legend = paste(unique(V(graph.gvar)$koppen), paste0("(", freq_vec, ")")),
           pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)")
    dev.off()

    ## graphicalVAR plot
    model_name <- paste0("graphicalVAR ", graph_type)
    pdf(file = file.path(output_dir, paste0(output_prefix, "_", graph_type, "_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
    plot(spec_gvar$communities, spec_gvar$graph,
         col = col_vec, cex = 0.5, vertex.label = "", layout = lay_gvar,
         vertex.size = 8, main = paste0("Clustering of ", model_name, " Graph (Purity = ", round(purity_gvar, 3), ")"))
    legend("topright", legend = paste(unique(V(graph.gvar)$koppen), paste0("(", freq_vec, ")")),
           pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)", col = "black")
    dev.off()
  }

  purity_vec <- c(purity_gvar, purity_3nn, purity_7nn, purity_10nn)
  names(purity_vec) <- c("graphicalVAR", "3NN", "7NN", "10NN")
  return(purity_vec)
}


graphicalVAR_on_temp <- function(path_train,
                                 path_test,
                                 plot_temp = FALSE,
                                 run_prep = TRUE,
                                 num.cluster = 7,
                                 year = 2020,
                                 standardize = FALSE,
                                 forecast = FALSE,
                                 lags = 1:3,
                                 nLambda = 20,
                                 gamma = 0.5,
                                 scale = FALSE,
                                 centerWithin = FALSE,
                                 lambda_beta = NULL,
                                 lambda_kappa = NULL,
                                 lambda_min_beta = 0.05,
                                 lambda_min_kappa = 0.05,
                                 regularize_mat_beta = NULL,
                                 regularize_mat_kappa = NULL,
                                 penalize.diagonal = TRUE,
                                 refit = FALSE,
                                 use_refit_forecast = FALSE,
                                 graph_type = c("contemporaneous", "temporal"),
                                 include_knn = TRUE,
                                 save_plots = TRUE,
                                 output_dir = "results",
                                 verbose = FALSE,
                                 ...) {
  graph_type <- match.arg(graph_type)

  temp_list <- temp_detrend(
    path_train,
    path_test,
    plot_temp,
    run_prep = run_prep,
    standardize = standardize,
    forecast = forecast
  )
  temp.train <- temp_list[["train"]]
  temp.test <- temp_list[["test"]]
  colnames(temp.test) <- colnames(temp.train)

  if (verbose) {
    print("Preprocessing completed...")
  }

  gvar_fit <- fit_graphicalvar_temp(
    train = temp.train,
    lags = lags,
    nLambda = nLambda,
    gamma = gamma,
    scale = scale,
    centerWithin = centerWithin,
    lambda_beta = lambda_beta,
    lambda_kappa = lambda_kappa,
    lambda_min_beta = lambda_min_beta,
    lambda_min_kappa = lambda_min_kappa,
    regularize_mat_beta = regularize_mat_beta,
    regularize_mat_kappa = regularize_mat_kappa,
    penalize.diagonal = penalize.diagonal,
    refit = refit,
    verbose = verbose,
    ...
  )

  if (verbose) {
    print("graphicalVAR training completed...")
  }

  forecasts <- make_forecast_graphicalvar(
    model = gvar_fit,
    train = temp.train,
    test = temp.test,
    use_refit = use_refit_forecast
  )

  if (verbose) {
    print("Forecasts made...")
  }

  purity <- spec_cluster_graphicalvar(
    graphicalvar_model = gvar_fit,
    path_train = path_train,
    train = temp.train,
    K = num.cluster,
    year = year,
    include_knn = include_knn,
    use_refit = use_refit_forecast,
    graph_type = graph_type,
    save_plots = save_plots,
    output_dir = output_dir
  )

  list(
    fit = gvar_fit$fit,
    refit = gvar_fit$refit,
    lags = lags,
    tuning = gvar_fit$tuning,
    adjacency = graphicalvar_adjacency(gvar_fit, use_refit = use_refit_forecast, type = graph_type),
    forecasts = forecasts,
    train = temp.train,
    test = temp.test,
    purity = purity
  )
}


graphicalVAR_lag_sweep_on_temp <- function(path_train,
                                           path_test,
                                           plot_temp = FALSE,
                                           run_prep = TRUE,
                                           num.cluster = 7,
                                           year = 2020,
                                           standardize = FALSE,
                                           forecast = FALSE,
                                           lag_specs = list("lag1" = 1, "lag2" = 1:2, "lag3" = 1:3),
                                           nLambda = 20,
                                           gamma = 0.5,
                                           scale = FALSE,
                                           centerWithin = FALSE,
                                           lambda_beta = NULL,
                                           lambda_kappa = NULL,
                                           lambda_min_beta = 0.05,
                                           lambda_min_kappa = 0.05,
                                           regularize_mat_kappa = NULL,
                                           penalize.diagonal = TRUE,
                                           net_size_tar = NULL,
                                           refit = FALSE,
                                           use_refit_forecast = FALSE,
                                           graph_type = c("contemporaneous", "temporal"),
                                           include_knn = TRUE,
                                           save_plots = TRUE,
                                           n.cores.lag = 1,
                                           n.cores.beta = 1,
                                           beta_block_size = 1,
                                           output_dir = "results",
                                           verbose = FALSE,
                                           ...) {
  graph_type <- match.arg(graph_type)

  temp_list <- temp_detrend(
    path_train,
    path_test,
    plot_temp,
    run_prep = run_prep,
    standardize = standardize,
    forecast = forecast
  )
  temp.train <- temp_list[["train"]]
  temp.test <- temp_list[["test"]]
  colnames(temp.test) <- colnames(temp.train)

  if (verbose) {
    print("Preprocessing completed...")
  }

  lag_models <- vector("list", length(lag_specs))
  names(lag_models) <- names(lag_specs)
  ebic_by_lag <- data.frame(
    lag_name = names(lag_specs),
    lags = vapply(lag_specs, function(x) paste(x, collapse = ","), character(1)),
    EBIC = NA_real_,
    lambda_beta = NA_real_,
    lambda_kappa = NA_real_,
    stringsAsFactors = FALSE
  )
  forecasts_by_lag <- vector("list", length(lag_specs))
  names(forecasts_by_lag) <- names(lag_specs)
  extra_args <- list(...)

  lambda_grid_by_lag <- lapply(lag_specs, function(lags_current) {
    make_graphicalvar_lambda_grid(
      train = temp.train,
      lags = lags_current,
      nLambda = nLambda,
      scale = scale,
      centerWithin = centerWithin,
      lambda_min_beta = lambda_min_beta,
      lambda_min_kappa = lambda_min_kappa,
      penalize.diagonal = penalize.diagonal,
      lambda_beta = lambda_beta,
      lambda_kappa = lambda_kappa
    )
  })
  names(lambda_grid_by_lag) <- names(lag_specs)

  fit_tasks <- do.call(rbind, lapply(names(lag_specs), function(lag_name) {
    lags_current <- lag_specs[[lag_name]]
    beta_values <- lambda_grid_by_lag[[lag_name]]$lambda_beta
    beta_blocks <- split(beta_values, ceiling(seq_along(beta_values) / beta_block_size))
    data.frame(
      lag_name = lag_name,
      block_id = seq_along(beta_blocks),
      stringsAsFactors = FALSE
    )
  }))

  fit_one_beta_block <- function(lag_name, block_id) {
    lags_current <- lag_specs[[lag_name]]
    lambda_grid <- lambda_grid_by_lag[[lag_name]]
    beta_blocks <- split(
      lambda_grid$lambda_beta,
      ceiling(seq_along(lambda_grid$lambda_beta) / beta_block_size)
    )
    lambda_beta_current <- beta_blocks[[block_id]]

    if (verbose) {
      print(paste0(
        "Fitting graphicalVAR for lags = ",
        paste(lags_current, collapse = ","),
        ", beta block ",
        block_id,
        "..."
      ))
    }

    gvar_fit <- do.call(
      fit_graphicalvar_temp,
      c(
        list(
          train = temp.train,
          lags = lags_current,
          nLambda = nLambda,
          gamma = gamma,
          scale = scale,
          centerWithin = centerWithin,
          lambda_beta = lambda_beta_current,
          lambda_kappa = lambda_grid$lambda_kappa,
          lambda_min_beta = lambda_min_beta,
          lambda_min_kappa = lambda_min_kappa,
          regularize_mat_beta = NULL,
          regularize_mat_kappa = regularize_mat_kappa,
          penalize.diagonal = penalize.diagonal,
          refit = refit,
          verbose = verbose
        ),
        extra_args
      )
    )

    path_current <- gvar_fit$fit$path
    selected_idx <- which.min(path_current$ebic)
    ebic_current <- data.frame(
      lag_name = lag_name,
      lags = paste(lags_current, collapse = ","),
      block_id = block_id,
      EBIC = path_current$ebic[selected_idx],
      lambda_beta = path_current$beta[selected_idx],
      lambda_kappa = path_current$kappa[selected_idx],
      stringsAsFactors = FALSE
    )

    list(
      lag_name = lag_name,
      block_id = block_id,
      model = gvar_fit,
      ebic = ebic_current
    )
  }

  n.requested.workers <- max(1, n.cores.lag) * max(1, n.cores.beta)
  if (n.requested.workers > 1 && nrow(fit_tasks) > 1) {
    n.workers <- min(n.requested.workers, nrow(fit_tasks))
    cluster_type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"
    cl <- try(parallel::makeCluster(n.workers, type = cluster_type), silent = TRUE)

    if (inherits(cl, "try-error")) {
      warning("Could not create parallel cluster; falling back to sequential lag/beta sweep.")
      beta_block_results <- lapply(seq_len(nrow(fit_tasks)), function(task_idx) {
        fit_one_beta_block(fit_tasks$lag_name[task_idx], fit_tasks$block_id[task_idx])
      })
    } else {
      doParallel::registerDoParallel(cl)
      on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      }, add = TRUE)

      beta_block_results <- foreach(
        task_idx = seq_len(nrow(fit_tasks)),
        .packages = c("graphicalVAR"),
        .export = c(
          "fit_graphicalvar_temp",
          "make_graphicalvar_regularize_mat_beta",
          "refit_graphicalvar_beta",
          "make_graphicalvar_design",
          "split_graphicalvar_beta",
          "make_forecast_graphicalvar",
          "graphicalvar_pred",
          "predict_graphicalvar_one"
        )
      ) %dopar% {
        fit_one_beta_block(fit_tasks$lag_name[task_idx], fit_tasks$block_id[task_idx])
      }
    }
  } else {
    beta_block_results <- lapply(seq_len(nrow(fit_tasks)), function(task_idx) {
      fit_one_beta_block(fit_tasks$lag_name[task_idx], fit_tasks$block_id[task_idx])
    })
  }

  ebic_by_block <- do.call(rbind, lapply(beta_block_results, `[[`, "ebic"))

  for (lag_name in names(lag_specs)) {
    lag_block_results <- beta_block_results[vapply(beta_block_results, function(res) res$lag_name == lag_name, logical(1))]
    lag_block_ebic <- do.call(rbind, lapply(lag_block_results, `[[`, "ebic"))
    best_block_idx <- which.min(lag_block_ebic$EBIC)
    best_block_result <- lag_block_results[[best_block_idx]]
    combined_path <- do.call(rbind, lapply(lag_block_results, function(res) res$model$fit$path))
    combined_path$id <- seq_len(nrow(combined_path))
    best_model <- best_block_result$model
    best_model$fit$path <- combined_path

    lag_models[[lag_name]] <- best_model
    forecasts_by_lag[[lag_name]] <- make_forecast_graphicalvar(
      model = best_model,
      train = temp.train,
      test = temp.test,
      use_refit = use_refit_forecast
    )
    ebic_by_lag[ebic_by_lag$lag_name == lag_name, ] <- best_block_result$ebic[names(ebic_by_lag)]
  }

  selected_lag_name <- ebic_by_lag$lag_name[which.min(ebic_by_lag$EBIC)]
  selected_model <- lag_models[[selected_lag_name]]

  if (verbose) {
    print("Lag sweep completed...")
    print(ebic_by_lag)
    print(paste0("Selected lag specification: ", selected_lag_name))
  }

  forecasts <- forecasts_by_lag[[selected_lag_name]]

  if (verbose) {
    print("Final forecasts selected by global EBIC...")
  }

  purity <- spec_cluster_graphicalvar(
    graphicalvar_model = selected_model,
    path_train = path_train,
    train = temp.train,
    K = num.cluster,
    year = year,
    include_knn = include_knn,
    use_refit = use_refit_forecast,
    graph_type = graph_type,
    save_plots = save_plots,
    output_dir = output_dir
  )

  size_match <- select_graphicalvar_size_match(
    beta_block_results = beta_block_results,
    selected_lag_name = selected_lag_name,
    net_size_tar = net_size_tar,
    graph_type = graph_type
  )

  graphicalVAR_sim_size <- NULL
  if (!is.null(size_match)) {
    size_match_forecasts <- make_forecast_graphicalvar(
      model = size_match$model,
      train = temp.train,
      test = temp.test,
      use_refit = use_refit_forecast
    )

    size_match_purity <- spec_cluster_graphicalvar(
      graphicalvar_model = size_match$model,
      path_train = path_train,
      train = temp.train,
      K = num.cluster,
      year = year,
      include_knn = include_knn,
      use_refit = use_refit_forecast,
      graph_type = graph_type,
      save_plots = save_plots,
      output_prefix = "graphicalVAR_sim_size",
      output_dir = output_dir
    )

    graphicalVAR_sim_size <- list(
      fit = size_match$model$fit,
      refit = size_match$model$refit,
      lags = size_match$model$lags,
      selected_lag = selected_lag_name,
      selection = size_match$selection,
      candidates = size_match$candidates,
      tuning = size_match$model$tuning,
      adjacency = graphicalvar_adjacency(size_match$model, use_refit = use_refit_forecast, type = graph_type),
      forecasts = size_match_forecasts,
      purity = size_match_purity
    )

    if (verbose) {
      print("Size-matched graphicalVAR selected...")
      print(size_match$selection)
    }
  }

  list(
    fit = selected_model$fit,
    refit = selected_model$refit,
    lags = selected_model$lags,
    selected_lag = selected_lag_name,
    ebic_by_lag = ebic_by_lag,
    ebic_by_block = ebic_by_block,
    lag_models = lag_models,
    forecasts_by_lag = forecasts_by_lag,
    tuning = selected_model$tuning,
    adjacency = graphicalvar_adjacency(selected_model, use_refit = use_refit_forecast, type = graph_type),
    forecasts = forecasts,
    train = temp.train,
    test = temp.test,
    purity = purity,
    net_size_tar = net_size_tar,
    graphicalVAR_sim_size = graphicalVAR_sim_size
  )
}
