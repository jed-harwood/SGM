pred_n_ahead <- function(model, order = 1, obs, n.ahead = 1) {
  if (order == 1) {
    R1.est <- model[[1]]
    pred <- R1.est %*% obs

    if (n.ahead > 1) {
      for (k in 2:n.ahead) {
        pred <- R1.est %*% pred
      }
    }
  }

  if (order == 2) {
    R1.est <- model[[1]]
    R2.est <- model[[2]]

    roll.window <- obs
    pred <- R1.est %*% c(obs[1, ]) + R2.est %*% c(obs[2, ])
    roll.window <- rbind(c(pred), c(obs[1, ]))

    if (n.ahead > 1) {
      for (k in 2:n.ahead) {
        pred <- R1.est %*% c(roll.window[1, ]) + R2.est %*% roll.window[2, ]
        roll.window <- rbind(c(pred), c(roll.window[1, ]))
      }
    }
  }

  if (order == 3) {
    R1.est <- model[[1]]
    R2.est <- model[[2]]
    R3.est <- model[[3]]

    roll.window <- obs
    pred <- R1.est %*% c(obs[1, ]) + R2.est %*% c(obs[2, ]) +
      R3.est %*% c(obs[3, ])
    roll.window <- rbind(c(pred), c(obs[1, ]), c(obs[2, ]))

    if (n.ahead > 1) {
      for (k in 2:n.ahead) {
        pred <- R1.est %*% c(roll.window[1, ]) +
          R2.est %*% roll.window[2, ] +
          R3.est %*% c(roll.window[3, ])
        roll.window <- rbind(c(pred), c(roll.window[1, ]), c(roll.window[2, ]))
      }
    }
  }

  pred
}

#' Forecast using TAR-GAR filter matrices
#'
#' @param model A list containing TAR-GAR filter matrices.
#' @param testdata A numeric matrix of test observations.
#' @param order The autoregressive order of the model.
#' @param n.ahead Number of steps ahead to forecast.
#'
#' @returns A list containing predictions and forecasting error summaries.
#' @export
TARGAR_pred <- function(model, testdata, order = 1, n.ahead = 1) {
  n <- nrow(testdata)
  p <- nrow(model[[1]])

  n.valid <- n - order - n.ahead + 1
  predMat <- matrix(NA, nrow = n.valid, ncol = p)
  rownames(predMat) <- paste0("t", (order + n.ahead):n)

  MSPE.vec <- rep(NA, n.valid)
  counter <- 1
  mse.mat <- predMat
  denom.node.rmse <- NULL

  for (i in order:(n - n.ahead)) {
    predMat[counter, ] <- pred_n_ahead(
      model = model,
      obs = testdata[i:(i - order + 1), ],
      n.ahead = n.ahead,
      order = order
    )

    obs.true <- testdata[i + n.ahead, ]
    denom.node.rmse <- rbind(denom.node.rmse, obs.true)
    MSPE.vec[counter] <- sum(abs(predMat[counter, ] - obs.true)^2) /
      sum(abs(obs.true)^2)
    mse.mat[counter, ] <- predMat[counter, ] - obs.true
    counter <- counter + 1
  }

  MSPE <- mean(MSPE.vec, na.rm = TRUE)
  RMSE <- sqrt(MSPE)
  node.rmse <- sqrt(colMeans(mse.mat^2, na.rm = TRUE))
  node.rmae <- colSums(abs(mse.mat), na.rm = TRUE) /
    colSums(abs(denom.node.rmse), na.rm = TRUE)

  index.plt <- sample(seq_len(p), size = min(10, p), replace = FALSE)
  val_vec <- c()
  label_vec <- c()
  true_mat <- testdata[(order + n.ahead):n, ]
  true_mat_pred_times <- testdata[(order + n.ahead):n, index.plt, drop = FALSE]
  for (k in seq_along(index.plt)) {
    val_vec <- c(val_vec, predMat[, index.plt[k]])
    label_vec <- c(label_vec, rep(k, length(predMat[, index.plt[k]])))
  }
  plot_df <- data.frame(
    pred = val_vec,
    station = label_vec,
    true = c(true_mat_pred_times)
  )

  node.corr <- cor(predMat, true_mat, use = "pairwise.complete.obs")

  list(
    pred = predMat,
    mspe = MSPE,
    rmse = RMSE,
    rmse.node = summary(node.rmse),
    plot.df = plot_df,
    rmae.node = summary(node.rmae),
    node.corr = summary(diag(node.corr)),
    full.corr = node.corr
  )
}

#' Forecast using univariate AR(p) fits
#'
#' @param arList A list of fitted AR(p) objects, such as outputs from `ar.yw()`.
#' @param testdata A numeric matrix of test observations.
#' @param n.ahead Number of steps ahead to forecast.
#'
#' @returns A list containing predictions and forecasting error summaries.
#' @export
ARp_pred <- function(arList, testdata, n.ahead = 1) {
  p <- length(arList)
  n <- nrow(testdata)
  predMat <- matrix(NA, nrow = n, ncol = p)
  num.lags.vec <- vapply(arList, function(obj) length(obj$ar), numeric(1))

  for (i in seq_len(p)) {
    coefList <- arList[[i]]$ar
    num.lags <- num.lags.vec[i]

    for (l in seq_len(n)) {
      if (l + num.lags + n.ahead - 1 > n) {
        message(paste(
          "Note:", n - l + 1,
          "observations cannot be validated for component", i
        ))
        break
      }

      tempVec <- testdata[seq(l + num.lags - 1, l, by = -1), i]
      pred.i <- as.numeric(coefList %*% tempVec)

      if (n.ahead >= 2) {
        for (j in 2:n.ahead) {
          tempVec <- c(pred.i, tempVec[-num.lags])
          pred.i <- as.numeric(coefList %*% tempVec)
        }
      }

      predMat[l, i] <- pred.i
    }
  }

  MSPE.vec <- c()
  rmse.node <- matrix(NA, nrow = n, ncol = p)

  for (l in seq_len(n)) {
    if (any(is.na(predMat[l, ]))) {
      next
    }

    valid <- TRUE
    mspe.l <- 0

    for (i in seq_len(p)) {
      num.lags <- num.lags.vec[i]
      idx.true <- l + num.lags + n.ahead - 1

      if (idx.true > n) {
        valid <- FALSE
        break
      }

      obs.pred <- predMat[l, i]
      obs.true <- testdata[idx.true, i]
      mspe.l <- mspe.l + (obs.pred - obs.true)^2
      if (i == p) {
        mspe.l <- mspe.l / sum(abs(testdata[idx.true, ])^2)
      }

      rmse.node[l, i] <- (obs.true - obs.pred)^2
    }

    if (valid) {
      MSPE.vec <- c(MSPE.vec, mspe.l)
    }
  }

  MSPE <- mean(MSPE.vec)
  RMSE <- sqrt(MSPE)

  list(
    pred = predMat,
    mspe = MSPE,
    rmse = RMSE,
    num.lags = num.lags.vec,
    rmse.node = rmse.node
  )
}

gvar_autocov <- function(data, order = 0) {
  n <- nrow(data)
  p <- ncol(data)

  if (order < 0 || order > n - 1) {
    stop("`order` must be an integer between 0 and n - 1.")
  }

  data.m <- matrix(apply(data, 2, mean), n, p, byrow = TRUE)
  data.c <- data - data.m

  data.c.pre <- data.c[1:(n - order), ]
  data.c.lag <- data.c[(1 + order):n, ]
  (t(data.c.pre) %*% data.c.lag) / (n - order)
}

#' Fit an Isufi-style graph VAR model
#'
#' @param data An n by p numeric data matrix.
#' @param A A p by p adjacency matrix.
#' @param q The autoregressive order.
#' @param L_q_vec Polynomial orders in the graph Laplacian for each lag.
#'
#' @returns A list with graph filter coefficients, filter matrices, and AIC.
#' @export
GVAR_fit <- function(data, A, q, L_q_vec) {
  p <- nrow(A)
  n <- nrow(data)

  L <- diag(apply(A, 1, sum)) - A
  L.list <- list()
  temp <- eigen(L)
  Li_and_Lj <- sort(L_q_vec, decreasing = TRUE)[1:2]

  L_m <- 2 * Li_and_Lj[1]
  for (j in 1:(L_m + 1)) {
    L.list[[j]] <- temp$vectors %*% diag(temp$values^(j - 1)) %*%
      t(temp$vectors)
  }

  AutoCov.List <- list()
  for (i in 0:q) {
    AutoCov.List[[i + 1]] <- gvar_autocov(data, order = i)
  }

  D <- matrix(NA, nrow = sum(L_q_vec) + q, ncol = sum(L_q_vec) + q)
  d <- rep(NA, sum(L_q_vec) + q)

  for (i in 1:q) {
    for (j in 1:(L_q_vec[i] + 1)) {
      offset <- if (i == 1) 0 else sum(L_q_vec[1:(i - 1)] + 1)
      d[offset + j] <- 2 * sum(diag(AutoCov.List[[i + 1]] %*% L.list[[j]]))
    }
  }

  for (i in 1:q) {
    for (j in 1:q) {
      lag <- abs(j - i)
      R_ij <- AutoCov.List[[lag + 1]]
      D_ij <- matrix(NA, nrow = L_q_vec[i] + 1, ncol = L_q_vec[j] + 1)

      for (k in 1:(L_q_vec[i] + 1)) {
        for (l in 1:(L_q_vec[j] + 1)) {
          D_ij[k, l] <- 2 * sum(diag(R_ij %*% L.list[[k + l - 1]]))
        }
      }

      offset.row <- if (i == 1) 0 else sum(L_q_vec[1:(i - 1)] + 1)
      offset.col <- if (j == 1) 0 else sum(L_q_vec[1:(j - 1)] + 1)

      idx.row <- (offset.row + 1):(offset.row + L_q_vec[i] + 1)
      idx.col <- (offset.col + 1):(offset.col + L_q_vec[j] + 1)
      D[idx.row, idx.col] <- D_ij
    }
  }

  A.qp <- diag(1, ncol(D))
  b <- rep(-1e10, ncol(D))
  result <- try(quadprog::solve.QP(Dmat = D, dvec = d, Amat = A.qp, bvec = b))
  if (inherits(result, "try-error")) {
    return(NULL)
  }

  psi <- result$unconstrained.solution
  filterList <- vector(mode = "list", length = q)
  for (j in 1:q) {
    offset <- if (j == 1) 0 else sum(L_q_vec[1:(j - 1)] + 1)
    temp <- matrix(0, nrow = p, ncol = p)
    for (l in 1:(L_q_vec[j] + 1)) {
      temp <- temp + psi[l + offset] * L.list[[l]]
    }
    filterList[[j]] <- temp
  }

  predMat <- matrix(NA, nrow = n, ncol = p)
  for (i in (q + 1):n) {
    temp <- rep(0, p)
    pred.window <- data[(i - 1):(i - q), ]
    if (q == 1) {
      pred.window <- t(as.matrix(pred.window))
    }

    for (j in 1:q) {
      temp <- temp + filterList[[j]] %*% pred.window[j, ]
    }

    predMat[i, ] <- temp
  }

  residMat <- data - predMat
  residMat <- residMat[-c(1:q), ]

  loglikeTerm <- determinant(var(residMat), logarithm = TRUE)$modulus[1]
  penalty <- 2 * sum(L_q_vec + 1) / (n - q)
  AIC <- loglikeTerm + penalty

  list(
    psi = psi,
    filters = filterList,
    qvec = L_q_vec,
    aic = AIC,
    det = loglikeTerm,
    penalty = penalty
  )
}

#' @rdname GVAR_fit
#' @export
G.VAR.Fit <- GVAR_fit
