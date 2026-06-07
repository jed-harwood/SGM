## Wrapper functions for modular TAR-GAR simulation studies.
##
## This script consolidates the TAR-GAR order 1, 2, and 3 simulation
## workflows into reusable pieces while preserving the original defaults.
## The eBIC parameter count is centralized as:
##   ar_order * (q + 1) + 1 + d + net.size
## where d is the node dimension.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

meanNA <- function(x) {
  mean(x, na.rm = TRUE)
}

setup_targar_simulation_parallel <- function(num.thread) {
  if (num.thread > 1) {
    if (!requireNamespace("doParallel", quietly = TRUE) ||
        !requireNamespace("foreach", quietly = TRUE)) {
      stop("Packages `doParallel` and `foreach` are required for parallel simulation.")
    }
  }
  invisible(num.thread)
}

legacy_min_index <- function(x) {
  idx = which(x == min(x, na.rm = TRUE), arr.ind = TRUE)
  ## Preserve the original simulation scripts' indexing. When eBIC has tied
  ## minima, the old code used index.c[1] and index.c[2] on the arr.ind matrix.
  matrix(c(idx[1], idx[2]), nrow = 1, dimnames = list(NULL, colnames(idx)))
}

relative_sq_error <- function(est, truth) {
  sum((est - truth)^2) / sum(truth^2)
}

rand_graph <- function(d, edge.prob, self.prob = edge.prob, min = 0.5,
                       max = 1, selfloop = FALSE, isolate = FALSE) {
  A = matrix(0, d, d)
  A[upper.tri(A)] = runif(d * (d - 1) / 2, min = min, max = max) *
    sample(0:1, d * (d - 1) / 2, replace = TRUE,
           prob = c(1 - edge.prob, edge.prob))
  A = A + t(A)

  if (selfloop) {
    diag(A) = runif(d, min = min, max = max) *
      sample(0:1, d, replace = TRUE, prob = c(1 - self.prob, self.prob))
  }

  if (!isolate) {
    deg = apply(A > 0, 1, sum)
    for (i in 1:(d - 1)) {
      if (deg[i] == 0) {
        j = sample((i + 1):d, 1)
        A[i, j] = runif(1, min = min, max = max)
        A[j, i] = A[i, j]
        deg[j] = 1
      }
    }
  }

  A
}

laplacian <- function(A) {
  diag(apply(A, 1, sum)) - A
}

laplacian_norm <- function(A) {
  d = nrow(A)
  deg = apply(A, 1, sum)
  L.norm = matrix(0, d, d)

  for (i in seq_len(d)) {
    if (deg[i] > 0) {
      L.norm[i, i] = 1 - A[i, i] / deg[i]
    }
    for (j in seq_len(d)) {
      if ((j != i) && deg[j] > 0 && A[i, j] > 0) {
        L.norm[i, j] = -A[i, j] / sqrt(deg[i] * deg[j])
      }
    }
  }

  L.norm
}

targar_poly_spec <- function(eta, lambda, ar_order, q) {
  specs = matrix(0, nrow = length(lambda), ncol = ar_order)
  for (r in seq_len(ar_order)) {
    block = eta[((r - 1) * (q + 1) + 1):(r * (q + 1))]
    for (j in 0:q) {
      specs[, r] = specs[, r] + block[j + 1] * lambda^j
    }
  }
  specs
}

eta_comparison <- function(eta, theta1, ar_order, q) {
  eta.comp = rep(NA_real_, ar_order * (q + 1))
  for (r in seq_len(ar_order)) {
    offset = (r - 1) * (q + 1)
    for (j in 0:q) {
      eta.comp[offset + j + 1] = eta[offset + j + 1] / theta1^j
    }
  }
  eta.comp
}

targar_truth <- function(eta, theta0, theta1, L, ar_order) {
  d = nrow(L)
  q = length(eta) / ar_order - 1
  eig = eigen(L, symmetric = TRUE)
  lambda = eig$values
  Q = t(eig$vectors)
  specs = targar_poly_spec(eta, lambda, ar_order = ar_order, q = q)

  Omega = t(Q) %*% diag((theta0 + theta1 * lambda)^2, d) %*% Q
  Sigma = t(Q) %*% diag((theta0 + theta1 * lambda)^(-2), d) %*% Q
  R.list = vector("list", ar_order)
  for (r in seq_len(ar_order)) {
    R.list[[r]] = t(Q) %*% diag(specs[, r], d) %*% Q
  }
  names(R.list) = paste0("R", seq_len(ar_order))

  if (ar_order == 1) {
    if (!all(abs(specs[, 1]) < 1)) {
      stop("TAR-GAR(1) process is not stationary for the supplied eta.")
    }
    Gamma0 = t(Q) %*%
      diag((theta0 + theta1 * lambda)^(-2) *
             (1 - specs[, 1]^2)^(-1), d) %*% Q
    Gamma1 = t(Q) %*%
      diag((theta0 + theta1 * lambda)^(-2) *
             (1 - specs[, 1]^2)^(-1) * specs[, 1], d) %*% Q
    return(c(list(Gamma0 = Gamma0, Gamma1 = Gamma1, Omega = Omega,
                  Sigma = Sigma), R.list))
  }

  if (ar_order == 2) {
    R1.spec = specs[, 1]
    R2.spec = specs[, 2]
    if ((!all(abs(R2.spec) < 1) &&
         (!all((R1.spec + R2.spec) < 1)) &&
         (!all((R2.spec - R1.spec) < 1)))) {
      warning("TAR-GAR(2) process may not be stationary.", call. = FALSE)
    }
    Gamma0 = t(Q) %*%
      diag((1 - R1.spec^2 -
              2 * R2.spec * (1 - R2.spec)^(-1) * R1.spec^2 -
              R2.spec^2)^(-1) * (theta0 + theta1 * lambda)^(-2), d) %*%
      Q
    Gamma1 = t(Q) %*%
      diag((1 - R2.spec)^(-1) * R1.spec *
             (1 - R1.spec^2 -
                2 * R2.spec * (1 - R2.spec)^(-1) * R1.spec^2 -
                R2.spec^2)^(-1) *
             (theta0 + theta1 * lambda)^(-2), d) %*% Q
    Gamma2 = R.list[[1]] %*% Gamma1 + R.list[[2]] %*% Gamma0
    Gamma0.tilde = rbind(cbind(Gamma0, Gamma1),
                         cbind(Gamma1, Gamma0))
    return(c(list(Gamma0 = Gamma0, Gamma1 = Gamma1, Gamma2 = Gamma2,
                  Gamma0.tilde = Gamma0.tilde, Omega = Omega,
                  Sigma = Sigma), R.list))
  }

  if (ar_order == 3) {
    R1.spec = specs[, 1]
    R2.spec = specs[, 2]
    R3.spec = specs[, 3]
    if ((!all(abs(R3.spec) < 1) ||
         (!all((-R1.spec - R2.spec) > R3.spec - 1)) ||
         (!all((R1.spec - R2.spec) > -R3.spec - 1))) ||
        (!all((R1.spec * R3.spec + R2.spec) > R3.spec^2 - 1))) {
      warning("TAR-GAR(3) process may not be stationary.", call. = FALSE)
    }

    M.spec = 1 - R1.spec^2 - R2.spec^2 - R3.spec^2 -
      (2 * R1.spec * R2.spec + 2 * R2.spec * R3.spec) *
      ((1 - R2.spec - R3.spec * (R1.spec + R2.spec))^(-1) *
         (R1.spec + R3.spec * R2.spec)) -
      2 * (R1.spec * R3.spec) * (R1.spec + R3.spec) *
      (1 - R2.spec - R3.spec * (R1.spec + R2.spec))^(-1) *
      (R1.spec + R2.spec * R3.spec) -
      2 * R1.spec * R2.spec * R3.spec
    M.spec.neg0.5 = numeric(d)
    M.spec.neg0.5[M.spec > 1e-6] = M.spec[M.spec > 1e-6]^(-0.5)
    M.neg0.5 = t(Q) %*% diag(M.spec.neg0.5, d) %*% Q
    Gamma0 = M.neg0.5 %*% Sigma %*% M.neg0.5
    Gamma1 = t(Q) %*%
      diag((1 - R2.spec - R3.spec * (R1.spec + R2.spec))^(-1) *
             (R1.spec + R3.spec * R2.spec), d) %*% Q %*% Gamma0
    Gamma2 = (R.list[[1]] + R.list[[3]]) %*% Gamma1 +
      R.list[[2]] %*% Gamma0
    Gamma3 = R.list[[1]] %*% Gamma2 + R.list[[2]] %*% Gamma1 +
      R.list[[3]] %*% Gamma0
    Gamma0.tilde = rbind(cbind(Gamma0, Gamma1, Gamma2),
                         cbind(Gamma1, Gamma0, Gamma1),
                         cbind(Gamma2, Gamma1, Gamma0))
    Gamma0.tilde = (t(Gamma0.tilde) + Gamma0.tilde) / 2
    return(c(list(Gamma0 = Gamma0, Gamma1 = Gamma1, Gamma2 = Gamma2,
                  Gamma3 = Gamma3, Gamma0.tilde = Gamma0.tilde,
                  Omega = Omega, Sigma = Sigma), R.list))
  }

  stop("`ar_order` must be 1, 2, or 3.")
}

generate_targar_data <- function(n, truth, ar_order, n_rep) {
  if (!requireNamespace("mnormt", quietly = TRUE)) {
    stop("Package `mnormt` is required for TAR-GAR simulation.")
  }
  d = nrow(truth$Sigma)
  data.rep = vector("list", n_rep)

  for (i in seq_len(n_rep)) {
    U.c = mnormt::rmnorm(n, mean = rep(0, d), varcov = truth$Sigma)
    data.c = matrix(0, n, d)

    if (ar_order == 1) {
      data.c[1, ] = mnormt::rmnorm(1, mean = rep(0, d),
                                   varcov = truth$Gamma0)
      if (n > 1) {
        for (tt in 2:n) {
          data.c[tt, ] = data.c[tt - 1, ] %*% truth$R1 + U.c[tt, ]
        }
      }
    } else if (ar_order == 2) {
      temp.gen = mnormt::rmnorm(1, mean = rep(0, 2 * d),
                                varcov = truth$Gamma0.tilde)
      data.c[1, ] = temp.gen[(d + 1):(2 * d)]
      data.c[2, ] = temp.gen[1:d]
      if (n > 2) {
        for (tt in 3:n) {
          data.c[tt, ] = data.c[tt - 1, ] %*% truth$R1 +
            data.c[tt - 2, ] %*% truth$R2 + U.c[tt, ]
        }
      }
    } else if (ar_order == 3) {
      temp.gen = mnormt::rmnorm(1, mean = rep(0, 3 * d),
                                varcov = truth$Gamma0.tilde)
      data.c[1, ] = temp.gen[(2 * d + 1):(3 * d)]
      data.c[2, ] = temp.gen[(d + 1):(2 * d)]
      data.c[3, ] = temp.gen[1:d]
      if (n > 3) {
        for (tt in 4:n) {
          data.c[tt, ] = data.c[tt - 1, ] %*% truth$R1 +
            data.c[tt - 2, ] %*% truth$R2 +
            data.c[tt - 3, ] %*% truth$R3 + U.c[tt, ]
        }
      }
    }

    data.rep[[i]] = data.c
  }

  data.rep
}

targar_default_config <- function(ar_order = 1) {
  ar_order = as.integer(ar_order)
  if (!ar_order %in% 1:3) {
    stop("`ar_order` must be 1, 2, or 3.")
  }

  base = list(
    ar_order = ar_order,
    d = 100,
    model = "LN",
    n_rep = 100,
    generation_n = NULL,
    graph_seed = 1,
    edge.prob = NULL,
    graph_min = 0.5,
    graph_max = 1,
    selfloop = FALSE,
    isolate = FALSE,
    theta0 = 1,
    theta1 = 2,
    eta = NULL,
    lambda.C = c(1.5, 1, 0.5, 0.25, 0.1),
    C.thre = NULL,
    eps.thre = 1e-6,
    max_iter = 50000,
    deg_max_iter = 50000,
    lap_z_max_iter = 50000,
    eta_max_iter = 1000,
    output_dir = "TARGAR/results",
    save_results = TRUE,
    keep_fits = TRUE
  )

  if (ar_order == 1) {
    base$n = 500
    base$generation_n = base$n + 100
    base$q = 1
    base$num.thread = 25
    base$num.pass = 3
    base$stationary = TRUE
    base$case = "TARGAR_order1_q1"
  } else if (ar_order == 2) {
    base$n = 500
    base$generation_n = base$n
    base$q = 3
    base$num.thread = 20
    base$num.pass = 3
    base$stationary = TRUE
    base$case = "TARGAR_order2_q3"
  } else {
    base$n = 100
    base$generation_n = base$n
    base$q = 3
    base$num.thread = 25
    base$num.pass = 5
    base$stationary = TRUE
    base$case = "TARGAR_order3_q3"
  }

  base$edge.prob = 2 / base$d
  base$fit_ar_order = base$ar_order
  base$fit_q = base$q
  base
}

targar_on_targar_default_config <- function(dgp_ar_order = 1) {
  dgp_ar_order = as.integer(dgp_ar_order)
  config = targar_default_config(dgp_ar_order)

  config$fit_ar_order = 1
  config$fit_q = 1
  config$num.pass = 3
  config$stationary = TRUE

  if (dgp_ar_order == 1) {
    config$n = 100
    config$generation_n = config$n
    config$q = 3
    config$num.thread = 25
    config$case = "TARGAR_p1q1_on_TARGAR_order1_q3"
  } else if (dgp_ar_order == 2) {
    config$n = 500
    config$generation_n = config$n
    config$q = 3
    config$num.thread = 20
    config$case = "TARGAR_p1q1_on_TARGAR_order2_q3"
  } else if (dgp_ar_order == 3) {
    config$n = 500
    config$generation_n = config$n
    config$q = 3
    config$num.thread = 25
    config$case = "TARGAR_p1q1_on_TARGAR_order3_q3"
  } else {
    stop("`dgp_ar_order` must be 1, 2, or 3.")
  }

  config$edge.prob = 2 / config$d
  config
}

targar_threshold_constants <- function(d) {
  if (d == 100) {
    exp(seq(log(1), log(0.05), length.out = 10))
  } else if (d == 250) {
    exp(seq(log(1), log(0.075), length.out = 10))
  } else {
    exp(seq(log(1), log(0.1), length.out = 10))
  }
}

prepare_targar_graph <- function(config) {
  d = config$d
  set.seed(config$graph_seed)
  A = rand_graph(d, config$edge.prob, min = config$graph_min,
                 max = config$graph_max, selfloop = config$selfloop,
                 isolate = config$isolate)
  net.tr = A > 0
  deg = apply(A, 1, sum)

  if (config$model == "L") {
    L = laplacian(A)
    v0 = rep(1, d)
  } else {
    L = laplacian_norm(A)
    v0 = matrix(sqrt(deg) / sqrt(sum(deg)), d, 1)
  }
  L.lam = eigen(L, symmetric = TRUE)$values
  lambda.max = max(L.lam)

  list(A = A, net.tr = net.tr, deg = deg, L = L, L.lam = L.lam,
       lambda.max = lambda.max, v0 = v0)
}

prepare_targar_simulation <- function(config, graph_setup = NULL) {
  d = config$d
  n = config$n
  generation_n = config$generation_n %||% n
  ar_order = config$ar_order
  fit_ar_order = config$fit_ar_order %||% ar_order
  q = config$q
  n_eff = n - fit_ar_order
  if (n_eff <= 0) {
    stop("`n - fit_ar_order` must be positive.")
  }
  if (generation_n < n) {
    stop("`generation_n` must be at least `n`.")
  }

  graph_setup = graph_setup %||% prepare_targar_graph(config)
  A = graph_setup$A
  net.tr = graph_setup$net.tr
  deg = graph_setup$deg
  L = graph_setup$L
  L.lam = graph_setup$L.lam
  lambda.max = graph_setup$lambda.max
  v0 = graph_setup$v0

  if (is.null(config$eta)) {
    stop("`eta` must be supplied as a numeric vector after `lambda.max` is available.")
  }
  if (is.function(config$eta)) {
    stop("`eta` must be a numeric vector, not a function.")
  }
  eta = as.numeric(config$eta)
  if (length(eta) != ar_order * (q + 1)) {
    stop("`eta` must have length ar_order * (q + 1).")
  }
  if (any(!is.finite(eta))) {
    stop("`eta` must contain only finite numeric values.")
  }
  eta.comp = eta_comparison(eta, theta1 = config$theta1,
                            ar_order = ar_order, q = q)
  truth = targar_truth(eta, theta0 = config$theta0,
                       theta1 = config$theta1, L = L,
                       ar_order = ar_order)

  C.thre = config$C.thre %||% targar_threshold_constants(d)
  lambda.v = config$lambda.C * sqrt(log(d) / n_eff)
  net.thre = C.thre * sqrt(log(d) / n_eff)
  rho.v = pmax(lambda.v, 0.01)

  data_seed = config$data_seed %||% (5 * d + n)
  set.seed(data_seed)
  data.generated = generate_targar_data(generation_n, truth = truth,
                                        ar_order = ar_order,
                                        n_rep = config$n_rep)
  data = lapply(data.generated, function(x) {
    x[seq_len(n), , drop = FALSE]
  })

  list(config = config, A = A, net.tr = net.tr, deg = deg, L = L,
       L.lam = L.lam, v0 = v0, eta = eta, eta.comp = eta.comp,
       truth = truth, lambda.v = lambda.v, net.thre = net.thre,
       C.thre = C.thre, rho.v = rho.v, data = data, data_seed = data_seed,
       generation_n = generation_n, lambda.max = lambda.max, n_eff = n_eff)
}

fit_targar_replicates <- function(setup, use_parallel = TRUE) {
  if (!requireNamespace("SGM", quietly = TRUE)) {
    stop("Package `SGM` is required. Install the local package before running simulations.")
  }
  config = setup$config
  n_rep = config$n_rep
  fit_ar_order = config$fit_ar_order %||% config$ar_order
  fit_q = config$fit_q %||% config$q

  fit_one = function(i) {
    data.i = setup$data[[i]][seq_len(config$n), , drop = FALSE]
    SGM::TARGAR_fit(
      data.i,
      lambda.v = setup$lambda.v,
      net.thre = setup$net.thre,
      rho.v = setup$rho.v,
      num.pass = config$num.pass,
      model = config$model,
      p = fit_ar_order,
      q = fit_q,
      eps.thre = config$eps.thre,
      max_iter = config$max_iter,
      deg_max_iter = config$deg_max_iter,
      lap_z_max_iter = config$lap_z_max_iter,
      eta_max_iter = config$eta_max_iter,
      stationary = config$stationary,
      verbose = FALSE
    )
  }

  if (use_parallel && config$num.thread > 1) {
    if (!requireNamespace("doParallel", quietly = TRUE) ||
        !requireNamespace("foreach", quietly = TRUE)) {
      stop("Packages `doParallel` and `foreach` are required for parallel simulation.")
    }
    cl = parallel::makeCluster(config$num.thread)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    foreach::`%dopar%`(
      foreach::foreach(i = seq_len(n_rep), .packages = "SGM"),
      {
        fit_one(i)
      }
    )
  } else {
    lapply(seq_len(n_rep), fit_one)
  }
}

targar_loglike <- function(S, theta0, theta1 = 1, L, n) {
  d = nrow(S)
  temp = diag(theta0, d) + theta1 * L
  temp1 = sum(diag(S %*% temp %*% temp))
  temp2 = determinant(temp, logarithm = TRUE)$modulus[[1]]
  -n / 2 * temp1 + n * temp2 - n * d / 2 * log(2 * pi)
}

targar_bic <- function(loglike, n, k) {
  k * log(n) - 2 * loglike
}

targar_ebic_parameter_count <- function(ar_order, q, d, net.size) {
  ar_order * (q + 1) + 1 + d + net.size
}

targar_ebic_penalty <- function(d, net.size, gamma) {
  P.total = d * (d - 1) / 2
  2 * gamma * (lfactorial(P.total) - lfactorial(net.size) -
                 lfactorial(P.total - net.size))
}

extract_refit_L <- function(item, use_0S = TRUE) {
  if (use_0S && !is.null(item$result.0S)) {
    return(item$result.0S$L)
  }
  if (!is.null(item$result.0.post)) {
    return(item$result.0.post$L)
  }
  item$L
}

extract_refit_theta0 <- function(item, use_0S = TRUE) {
  if (use_0S && !is.null(item$result.0S)) {
    return(item$result.0S$theta0)
  }
  if (!is.null(item$theta0)) {
    return(item$theta0)
  }
  item$theta0.ini
}

extract_refit_theta1 <- function(item, use_0S = TRUE) {
  if (use_0S && !is.null(item$result.0S)) {
    return(item$result.0S$theta1)
  }
  if (!is.null(item$result.0.post)) {
    return(item$result.0.post$theta1)
  }
  item$theta1
}

extract_refit_R <- function(item, r) {
  item[[paste0("R", r, ".0S")]] %||%
    if (!is.null(item$R_list)) item$R_list[[r]] else NULL
}

extract_targar_metrics <- function(results, setup) {
  config = setup$config
  n_rep = config$n_rep
  n_lambda = length(setup$lambda.v)
  n_thre = length(setup$net.thre)
  ar_order = config$ar_order
  q = config$q
  fit_ar_order = config$fit_ar_order %||% ar_order
  fit_q = config$fit_q %||% q
  d = config$d
  n_eff = setup$n_eff

  dims = c(n_rep, n_lambda, n_thre)
  conv.refit = array(NA, dim = dims)
  conv.refit.post = array(NA, dim = dims)
  theta.0S.err = array(NA, dim = dims)
  L.refit.0S.err = array(NA, dim = dims)
  v0.0S.err = array(NA, dim = dims)
  net.size.refit = array(NA, dim = dims)
  power.refit = array(NA, dim = dims)
  fdr.refit = array(NA, dim = dims)
  time.refit = array(NA, dim = dims)
  eta.refit.0S.err = array(NA, dim = c(dims, fit_ar_order * (fit_q + 1)))
  R.refit.0S.err = vector("list", ar_order)
  for (r in seq_len(ar_order)) {
    R.refit.0S.err[[r]] = array(NA, dim = dims)
  }
  names(R.refit.0S.err) =
    paste0("R", seq_len(ar_order), ".refit.0S.err")

  log.like.0S = array(NA, dim = dims)
  bic.0S = array(Inf, dim = dims)
  ebic.0S = array(Inf, dim = dims)
  ebic.param.0S = array(NA, dim = dims)

  gamma = if (d / n_eff > 0.5) 1 else 0.5

  for (i in seq_len(n_rep)) {
    for (j in seq_len(n_lambda)) {
      for (k in seq_len(n_thre)) {
        item = results[[i]]$refit[[j]][[k]]
        A.c = item$A.net
        net.size.c = sum(A.c) / 2
        net.size.refit[i, j, k] = net.size.c
        time.refit[i, j, k] = item$time
        conv.refit.post[i, j, k] = isTRUE(item$conv.2)
        conv.0S = isTRUE(item$conv.3) ||
          (config$model == "L" && isTRUE(item$conv.2))
        conv.refit[i, j, k] = conv.0S

        if (isTRUE(item$conv.2)) {
          power.refit[i, j, k] = sum(A.c * setup$net.tr) / sum(setup$net.tr)
          fdr.refit[i, j, k] = sum(A.c * (1 - setup$net.tr)) / sum(A.c)
        }

        if (conv.0S) {
          L.est = extract_refit_L(item, use_0S = TRUE)
          theta0.est = extract_refit_theta0(item, use_0S = TRUE)
          theta1.est = extract_refit_theta1(item, use_0S = TRUE)
          theta.0S.err[i, j, k] = (theta0.est - config$theta0)^2
          L.refit.0S.err[i, j, k] =
            relative_sq_error(L.est * theta1.est, setup$L * config$theta1)
          v0.0S.err[i, j, k] = sum((item$v0.est - setup$v0)^2)

          for (r in seq_len(ar_order)) {
            R.est = if (r <= fit_ar_order) {
              extract_refit_R(item, r)
            } else {
              matrix(0, d, d)
            }
            if (!is.null(R.est) && !is.null(setup$truth[[paste0("R", r)]])) {
              R.refit.0S.err[[r]][i, j, k] =
                relative_sq_error(R.est, setup$truth[[paste0("R", r)]])
            }
          }
          if (!is.null(item$eta)) {
            eta.truth = setup$eta.comp[seq_len(fit_ar_order * (fit_q + 1))]
            eta.refit.0S.err[i, j, k, ] = (item$eta - eta.truth)^2
          }

          log.like.0S[i, j, k] =
            targar_loglike(S = item$S, theta0 = theta0.est,
                           theta1 = theta1.est, L = L.est, n = n_eff)
          ebic.param.0S[i, j, k] =
            targar_ebic_parameter_count(ar_order = fit_ar_order, q = fit_q,
                                        d = d, net.size = net.size.c)
          bic.0S[i, j, k] =
            targar_bic(log.like.0S[i, j, k], n_eff, ebic.param.0S[i, j, k])
          ebic.0S[i, j, k] = bic.0S[i, j, k] +
            targar_ebic_penalty(d = d, net.size = net.size.c, gamma = gamma)
        }
      }
    }
  }

  selected = select_targar_0S_ebic_metrics(
    ebic.0S = ebic.0S,
    net.size.refit = net.size.refit,
    L.refit.0S.err = L.refit.0S.err,
    theta.0S.err = theta.0S.err,
    R.refit.0S.err = R.refit.0S.err,
    eta.refit.0S.err = eta.refit.0S.err,
    power.refit = power.refit,
    fdr.refit = fdr.refit,
    v0.0S.err = v0.0S.err
  )

  out = c(list(
    conv.refit = conv.refit,
    conv.refit.post = conv.refit.post,
    theta.0S.err = theta.0S.err,
    L.refit.0S.err = L.refit.0S.err,
    eta.refit.0S.err = eta.refit.0S.err,
    v0.0S.err = v0.0S.err,
    net.size.refit = net.size.refit,
    power.refit = power.refit,
    fdr.refit = fdr.refit,
    time.refit = time.refit,
    log.like.0S = log.like.0S,
    bic.0S = bic.0S,
    ebic.0S = ebic.0S,
    ebic.param.0S = ebic.param.0S,
    ebic.gamma = gamma,
    dgp_ar_order = ar_order,
    dgp_q = q,
    fit_ar_order = fit_ar_order,
    fit_q = fit_q,
    ebic.parameter.formula = "fit_ar_order * (fit_q + 1) + 1 + d + net.size"
  ), R.refit.0S.err, selected)
  out
}

select_targar_0S_ebic_metrics <- function(ebic.0S, net.size.refit,
                                          L.refit.0S.err, theta.0S.err,
                                          R.refit.0S.err, eta.refit.0S.err,
                                          power.refit, fdr.refit,
                                          v0.0S.err) {
  n_rep = dim(ebic.0S)[1]
  ar_order = length(R.refit.0S.err)
  eta_len = dim(eta.refit.0S.err)[4]

  size.ebic.0S = rep(NA_real_, n_rep)
  L.0S.ebic.err = rep(NA_real_, n_rep)
  theta.0S.ebic.err = rep(NA_real_, n_rep)
  eta.0S.ebic.err = matrix(NA_real_, nrow = n_rep, ncol = eta_len)
  power.0S.ebic = rep(NA_real_, n_rep)
  fdr.0S.ebic = rep(NA_real_, n_rep)
  F1.0S.ebic = rep(NA_real_, n_rep)
  v0.0S.ebic = rep(NA_real_, n_rep)
  ebic.0S.selec = rep(NA_real_, n_rep)
  ebic.0S.index = matrix(NA_integer_, nrow = n_rep, ncol = 2)
  colnames(ebic.0S.index) = c("lambda", "net.thre")

  R.0S.ebic.err = vector("list", ar_order)
  for (r in seq_len(ar_order)) {
    R.0S.ebic.err[[r]] = rep(NA_real_, n_rep)
  }
  names(R.0S.ebic.err) = paste0("R", seq_len(ar_order), ".0S.ebic.err")

  for (i in seq_len(n_rep)) {
    ebic.i = matrix(ebic.0S[i, , ],
                    nrow = dim(ebic.0S)[2],
                    ncol = dim(ebic.0S)[3])
    if (!all(is.infinite(ebic.i))) {
      idx = legacy_min_index(ebic.i)
      j = idx[1, 1]
      k = idx[1, 2]
      ebic.0S.index[i, ] = c(j, k)
      size.ebic.0S[i] = net.size.refit[i, j, k]
      ebic.0S.selec[i] = ebic.0S[i, j, k]
      L.0S.ebic.err[i] = L.refit.0S.err[i, j, k]
      theta.0S.ebic.err[i] = theta.0S.err[i, j, k]
      eta.0S.ebic.err[i, ] = eta.refit.0S.err[i, j, k, ]
      power.0S.ebic[i] = power.refit[i, j, k]
      fdr.0S.ebic[i] = fdr.refit[i, j, k]
      F1.0S.ebic[i] = 2 * power.refit[i, j, k] *
        (1 - fdr.refit[i, j, k]) /
        (1 - fdr.refit[i, j, k] + power.refit[i, j, k])
      v0.0S.ebic[i] = v0.0S.err[i, j, k]
      for (r in seq_len(ar_order)) {
        R.0S.ebic.err[[r]][i] = R.refit.0S.err[[r]][i, j, k]
      }
    }
  }

  c(list(size.ebic.0S = size.ebic.0S,
         L.0S.ebic.err = L.0S.ebic.err,
         theta.0S.ebic.err = theta.0S.ebic.err,
         eta.0S.ebic.err = eta.0S.ebic.err,
         power.0S.ebic = power.0S.ebic,
         fdr.0S.ebic = fdr.0S.ebic,
         F1.0S.ebic = F1.0S.ebic,
         v0.0S.ebic = v0.0S.ebic,
         ebic.0S.selec = ebic.0S.selec,
         ebic.0S.index = ebic.0S.index),
    R.0S.ebic.err)
}

summarize_targar_metrics <- function(metrics) {
  selected_names = grep("\\.0S\\.ebic(\\.err)?$|^ebic\\.0S\\.selec$|^size\\.ebic\\.0S$|^power\\.0S\\.ebic$|^fdr\\.0S\\.ebic$|^F1\\.0S\\.ebic$|^v0\\.0S\\.ebic$",
                        names(metrics), value = TRUE)
  summaries = lapply(metrics[selected_names], function(x) {
    if (is.matrix(x) || is.array(x)) {
      apply(x, length(dim(x)), summary, na.rm = TRUE)
    } else {
      summary(x, na.rm = TRUE)
    }
  })
  summaries$mean.ebic.0S.selec = meanNA(metrics$ebic.0S.selec)
  summaries$mean.ebic.0S.path = apply(metrics$ebic.0S, c(2, 3), meanNA)
  R_names = grep("^R[0-9]+\\.0S\\.ebic\\.err$", names(metrics), value = TRUE)
  summaries$R.0S.ebic.err = lapply(metrics[R_names], summary, na.rm = TRUE)
  summaries
}

run_targar_simulation <- function(config = targar_default_config(1),
                                  use_parallel = TRUE) {
  setup = prepare_targar_simulation(config)
  results = fit_targar_replicates(setup, use_parallel = use_parallel)
  metrics = extract_targar_metrics(results, setup)
  summaries = summarize_targar_metrics(metrics)

  out = list(setup = setup,
             results = if (isTRUE(config$keep_fits)) results else NULL,
             metrics = metrics,
             summaries = summaries)

  if (isTRUE(config$save_results)) {
    if (!dir.exists(config$output_dir)) {
      dir.create(config$output_dir, recursive = TRUE)
    }
    filename = file.path(
      config$output_dir,
      paste0(config$case, "_d", config$d, "_n", config$n,
             "_model", config$model, "_rep", config$n_rep,
             "_modular.RData")
    )
    save(out, file = filename)
    out$output_file = filename
  }

  out
}
