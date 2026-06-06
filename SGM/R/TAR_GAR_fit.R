###############################
## TAR-GAR(p,q) fitting
##
## Here p is the AR order, not the node dimension.
## The implementation supports AR orders p = 1, 2, and 3.
###############################

AutoCov <- function(data, order = 0, symmetrize = c("all", "order0", "none"),
                    transpose.lag = FALSE) {
  data = as.matrix(data)
  symmetrize = match.arg(symmetrize)
  n = nrow(data)
  d = ncol(data)

  if (order < 0 || order > n - 1) {
    stop("order must be an integer between {0,..., n-1}")
  }

  data.m = matrix(apply(data, 2, mean), n, d, byrow = TRUE)
  data.c = data - data.m
  data.c.pre = data.c[1:(n - order), , drop = FALSE]
  data.c.lag = data.c[(1 + order):n, , drop = FALSE]

  result = (t(data.c.pre) %*% data.c.lag) / (n - order)
  if (transpose.lag) {
    result = t(result)
  }
  if (symmetrize == "all" || (symmetrize == "order0" && order == 0)) {
    result = (result + t(result)) / 2
  }

  result
}

Stack.targar <- function(data, p) {
  data = as.matrix(data)
  n = nrow(data)
  d = ncol(data)

  if (!p %in% 1:3) {
    stop("p must be the AR order 1, 2, or 3")
  }
  if (n <= p) {
    stop("data must have more rows than the AR order p")
  }
  if (p == 1) {
    return(data)
  }

  dataStack = matrix(NA, nrow = n - p + 1, ncol = p * d)
  for (i in p:n) {
    dataStack[i - p + 1, ] = unlist(lapply(0:(p - 1), function(lag) {
      data[i - lag, ]
    }), use.names = FALSE)
  }

  dataStack
}

PseudoInverseSym <- function(S, eps.thre = 1e-6, power = -1) {
  S.eigen = eigen(S, symmetric = TRUE)
  S.lam = pmax(S.eigen$values, 0)
  S.lam.inv = numeric(length(S.lam))
  S.lam.inv[S.lam > eps.thre] = S.lam[S.lam > eps.thre]^power
  S.eigen$vectors %*% diag(S.lam.inv, length(S.lam.inv)) %*%
    t(S.eigen$vectors)
}

R.initial <- function(data, p, q = NULL, eps.thre = 1e-6) {
  data = as.matrix(data)
  d = ncol(data)

  if (!p %in% 1:3) {
    stop("p must be the AR order 1, 2, or 3")
  }

  if (p == 1) {
    Gamma0 = AutoCov(data, order = 0, symmetrize = "all")
    Gamma1 = AutoCov(data, order = 1, symmetrize = "all")
    Gamma0.neg0.5 = PseudoInverseSym(Gamma0, eps.thre = eps.thre,
                                     power = -0.5)
    R1 = Gamma0.neg0.5 %*% Gamma1 %*% Gamma0.neg0.5
    R1 = (R1 + t(R1)) / 2
    return(list(R1 = R1))
  }

  data.stacked = Stack.targar(data, p)
  Gamma0 = AutoCov(data.stacked, order = 0, symmetrize = "order0",
                   transpose.lag = TRUE)
  Gamma1 = AutoCov(data.stacked, order = 1, symmetrize = "order0",
                   transpose.lag = TRUE)
  Gamma0.inv = PseudoInverseSym(Gamma0, eps.thre = eps.thre, power = -1)

  ## For AR order p > 1, use the unsymmetrized companion estimate
  ## Gamma1 Gamma0^+ rather than the symmetrized whitening transform.
  A = Gamma1 %*% Gamma0.inv

  R.list = vector("list", p)
  for (k in 1:p) {
    R.list[[k]] = A[1:d, ((k - 1) * d + 1):(k * d), drop = FALSE]
  }
  names(R.list) = paste0("R", 1:p)

  R.list
}

S.theta0 <- function(data, R_list, p = length(R_list)) {
  data = as.matrix(data)
  n = nrow(data)

  if (n <= p) {
    stop("data must have more rows than the AR order p")
  }
  if (length(R_list) != p) {
    stop("R_list length must match AR order p")
  }

  residual = data[(p + 1):n, , drop = FALSE]
  for (k in 1:p) {
    lag.data = data[(p + 1 - k):(n - k), , drop = FALSE]
    residual = residual - lag.data %*% R_list[[k]]
  }

  S = stats::var(residual) * (n - p - 1) / (n - p)
  theta0 = (max(eigen(S, symmetric = TRUE)$values))^(-0.5)

  list(S = S, theta0 = theta0)
}

Phi.L <- function(lambda, q) {
  outer(lambda, 0:q, "^")
}

Weighted.d <- function(Phi, weights, diag.gamma) {
  2 * colSums(Phi * (weights * diag.gamma))
}

Weighted.D <- function(Phi, weights, diag.gamma) {
  2 * crossprod(Phi, Phi * (weights * diag.gamma))
}

Filter.from.eta <- function(eta.block, eigen.L, q) {
  lambda = eigen.L$values
  Phi = Phi.L(lambda, q)
  spec = as.numeric(Phi %*% eta.block)
  R = eigen.L$vectors %*% diag(spec, length(spec)) %*% t(eigen.L$vectors)
  (R + t(R)) / 2
}

Eta.data <- function(data, L) {
  data = as.matrix(data)
  n = nrow(data)
  d = ncol(data)
  data.m = matrix(apply(data, 2, mean), n, d, byrow = TRUE)
  data.c = data - data.m
  eigen.L = eigen(L, symmetric = TRUE)
  Q = t(eigen.L$vectors)
  data.tilde = data.c %*% t(Q)

  list(data.tilde = data.tilde, eigen.L = eigen.L, lambda = eigen.L$values)
}

QP.TAR.eta.p1 <- function(data, L, theta0, q = 1, eps = 1e-6,
                          stationary = TRUE) {
  eta.data = Eta.data(data, L)
  data.tilde = eta.data$data.tilde
  eigen.L = eta.data$eigen.L
  lambda = eta.data$lambda
  n = nrow(data.tilde)
  d = ncol(data.tilde)

  Gamma0 = t(data.tilde[2:n, , drop = FALSE]) %*%
    data.tilde[2:n, , drop = FALSE] / (n - 1)
  Gamma1 = t(data.tilde[2:n, , drop = FALSE]) %*%
    data.tilde[1:(n - 1), , drop = FALSE] / (n - 1)
  Gamma1 = (Gamma1 + t(Gamma1)) / 2

  Phi = Phi.L(lambda, q)
  weights = (theta0 + lambda)^2
  dvec = Weighted.d(Phi, weights, diag(Gamma1))
  Dmat = Weighted.D(Phi, weights, diag(Gamma0))

  Amat = matrix(0, nrow = q + 1, ncol = 2 * d)
  Amat[, 1:d] = t(Phi)
  Amat[, (d + 1):(2 * d)] = -t(Phi)
  bvec = rep(-1 + eps, 2 * d)

  result = try(quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat,
                                  bvec = bvec), silent = TRUE)
  if (inherits(result, "try-error")) {
    return(NULL)
  }

  eta = if (stationary) result$solution else result$unconstrained.solution
  R1 = Filter.from.eta(eta, eigen.L, q)

  list(result = result, eta = eta, R_list = list(R1 = R1), R1 = R1)
}

QP.TAR.eta.p2 <- function(data, L, theta0, q = 1, eps = 1e-6,
                          stationary = TRUE) {
  eta.data = Eta.data(data, L)
  data.tilde = eta.data$data.tilde
  eigen.L = eta.data$eigen.L
  lambda = eta.data$lambda
  n = nrow(data.tilde)
  d = ncol(data.tilde)

  y0 = data.tilde[3:n, , drop = FALSE]
  y1 = data.tilde[2:(n - 1), , drop = FALSE]
  y2 = data.tilde[1:(n - 2), , drop = FALSE]
  denom = n - 2

  Gamma0.1 = t(y1) %*% y1 / denom
  Gamma0.2 = t(y2) %*% y2 / denom
  Gamma1.01 = t(y0) %*% y1 / denom
  Gamma1.01 = (Gamma1.01 + t(Gamma1.01)) / 2
  Gamma1.12 = t(y1) %*% y2 / denom
  Gamma1.12 = (Gamma1.12 + t(Gamma1.12)) / 2
  Gamma2 = t(y0) %*% y2 / denom
  Gamma2 = (Gamma2 + t(Gamma2)) / 2

  Phi = Phi.L(lambda, q)
  weights = (theta0 + lambda)^2
  D11 = Weighted.D(Phi, weights, diag(Gamma0.1))
  D12 = Weighted.D(Phi, weights, diag(Gamma1.12))
  D22 = Weighted.D(Phi, weights, diag(Gamma0.2))
  Dmat = rbind(cbind(D11, D12), cbind(t(D12), D22))
  dvec = c(Weighted.d(Phi, weights, diag(Gamma1.01)),
           Weighted.d(Phi, weights, diag(Gamma2)))

  bsz = q + 1
  cons.block = t(Phi)
  Amat = matrix(0, nrow = 2 * bsz, ncol = 4 * d)
  Amat[1:bsz, 1:d] = -cons.block
  Amat[(bsz + 1):(2 * bsz), 1:d] = -cons.block
  Amat[1:bsz, (d + 1):(2 * d)] = cons.block
  Amat[(bsz + 1):(2 * bsz), (d + 1):(2 * d)] = -cons.block
  Amat[(bsz + 1):(2 * bsz), (2 * d + 1):(3 * d)] = -cons.block
  Amat[(bsz + 1):(2 * bsz), (3 * d + 1):(4 * d)] = cons.block
  bvec = rep(-1 + eps, 4 * d)

  result = try(quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat,
                                  bvec = bvec), silent = TRUE)
  if (inherits(result, "try-error")) {
    return(NULL)
  }

  eta = if (stationary) result$solution else result$unconstrained.solution
  eta1 = eta[1:bsz]
  eta2 = eta[(bsz + 1):(2 * bsz)]
  R1 = Filter.from.eta(eta1, eigen.L, q)
  R2 = Filter.from.eta(eta2, eigen.L, q)

  companion = rbind(cbind(R1, R2), cbind(diag(1, d), matrix(0, d, d)))
  if (!all(abs(eigen(companion, only.values = TRUE)$values) < 1)) {
    message("Not stationary!")
  }

  list(result = result, eta = eta, R_list = list(R1 = R1, R2 = R2),
       R1 = R1, R2 = R2)
}

D.maker.loss <- function(data, lambda, theta0, q = 1) {
  n = nrow(data)

  Gamma.01 = t(data[3:(n - 1), , drop = FALSE]) %*%
    data[4:n, , drop = FALSE] / (n - 3)
  Gamma.01 = (Gamma.01 + t(Gamma.01)) / 2
  Gamma.02 = t(data[2:(n - 2), , drop = FALSE]) %*%
    data[4:n, , drop = FALSE] / (n - 3)
  Gamma.02 = (Gamma.02 + t(Gamma.02)) / 2
  Gamma.03 = t(data[1:(n - 3), , drop = FALSE]) %*%
    data[4:n, , drop = FALSE] / (n - 3)
  Gamma.03 = (Gamma.03 + t(Gamma.03)) / 2

  D.list = vector(mode = "list", length = 3)
  for (i in 1:3) {
    D.list[[i]] = vector(mode = "list", length = 3)
    for (j in i:3) {
      D.list[[i]][[j]] =
        t(data[(4 - i):(n - i), , drop = FALSE]) %*%
        data[(4 - j):(n - j), , drop = FALSE] / (n - 3)
      if (i != j) {
        D.list[[i]][[j]] = (D.list[[i]][[j]] + t(D.list[[i]][[j]])) / 2
      }
    }
  }

  d.1 = rep(0, q + 1)
  d.2 = rep(0, q + 1)
  d.3 = rep(0, q + 1)

  D11 = matrix(0, nrow = q + 1, ncol = q + 1)
  D12 = D11
  D13 = D11
  D22 = D11
  D23 = D11
  D33 = D11

  for (j in 1:(q + 1)) {
    d.1[j] = 2 * sum(lambda^((j - 1)) * (theta0 + lambda)^2 *
                       diag(Gamma.01))
    d.2[j] = 2 * sum(lambda^((j - 1)) * (theta0 + lambda)^2 *
                       diag(Gamma.02))
    d.3[j] = 2 * sum(lambda^((j - 1)) * (theta0 + lambda)^2 *
                       diag(Gamma.03))
    for (k in 1:(q + 1)) {
      D11[j, k] = 2 * sum(lambda^((j - 1) + (k - 1)) *
                            (theta0 + lambda)^2 *
                            diag(D.list[[1]][[1]]))
      D12[j, k] = 2 * sum(lambda^((j - 1) + (k - 1)) *
                            (theta0 + lambda)^2 *
                            diag(D.list[[1]][[2]]))
      D13[j, k] = 2 * sum(lambda^((j - 1) + (k - 1)) *
                            (theta0 + lambda)^2 *
                            diag(D.list[[1]][[3]]))
      D22[j, k] = 2 * sum(lambda^((j - 1) + (k - 1)) *
                            (theta0 + lambda)^2 *
                            diag(D.list[[2]][[2]]))
      D23[j, k] = 2 * sum(lambda^((j - 1) + (k - 1)) *
                            (theta0 + lambda)^2 *
                            diag(D.list[[2]][[3]]))
      D33[j, k] = 2 * sum(lambda^((j - 1) + (k - 1)) *
                            (theta0 + lambda)^2 *
                            diag(D.list[[3]][[3]]))
    }
  }

  dvec = c(d.1, d.2, d.3)
  Dmat = rbind(cbind(D11, D12, D13),
               cbind(t(D12), D22, D23),
               cbind(t(D13), t(D23), D33))

  list(D = Dmat, d = dvec)
}

Pj <- function(lambda.val, eta.j) {
  q = length(eta.j) - 1
  sum(eta.j * lambda.val^(0:q))
}

constraint_fn <- function(eta, lambda_vals, q) {
  eta1 = eta[1:(q + 1)]
  eta2 = eta[(q + 2):(2 * q + 2)]
  eta3 = eta[(2 * q + 3):(3 * q + 3)]
  constraints = numeric()

  for (lam in lambda_vals) {
    P1 = Pj(lam, eta1)
    P2 = Pj(lam, eta2)
    P3 = Pj(lam, eta3)

    c1 = P1 + P2 + P3 - 1
    c2 = -P1 + P2 - P3 - 1
    c3 = -P3 * (P1 - P3) - P2 - 1
    c4 = P3 - 1
    c5 = -P3 - 1

    constraints = c(constraints, c1, c2, c3, c4, c5)
  }

  constraints
}

jacobian_fn <- function(eta, lambda_vals, q) {
  n_eta = length(eta)
  eta1.idx = 1:(q + 1)
  eta2.idx = (q + 2):(2 * q + 2)
  eta3.idx = (2 * q + 3):(3 * q + 3)
  J = matrix(0, nrow = 5 * length(lambda_vals), ncol = n_eta)

  for (i in seq_along(lambda_vals)) {
    lam = lambda_vals[i]
    powers = lam^(0:q)
    row.start = 5 * (i - 1)

    dP1 = numeric(n_eta)
    dP2 = numeric(n_eta)
    dP3 = numeric(n_eta)
    dP1[eta1.idx] = powers
    dP2[eta2.idx] = powers
    dP3[eta3.idx] = powers

    P1.val = sum(eta[eta1.idx] * powers)
    P3.val = sum(eta[eta3.idx] * powers)

    J[row.start + 1, ] = dP1 + dP2 + dP3
    J[row.start + 2, ] = -dP1 + dP2 - dP3
    J[row.start + 3, ] = -P3.val * dP1 - dP2 - P1.val * dP3 +
      2 * P3.val * dP3
    J[row.start + 4, ] = dP3
    J[row.start + 5, ] = -dP3
  }

  J
}

loss.eta <- function(D, d, eta) {
  as.numeric((1 / 2) * t(eta) %*% D %*% eta - d %*% eta)
}

eta.nlr.mod <- function(data, L, theta0, p, q, eps_thre = 1e-6,
                        max.iter = 1000, stationary = TRUE) {
  if (p != 3) {
    stop("NLP is only implemented for AR order p=3")
  }
  if (stationary && !requireNamespace("nloptr", quietly = TRUE)) {
    stop("Package `nloptr` is required to fit stationary TAR-GAR models with p = 3.")
  }

  eta.data = Eta.data(data, L)
  data.tilde = eta.data$data.tilde
  eigen.L = eta.data$eigen.L
  lambda = eta.data$lambda

  obj.loss = D.maker.loss(data = data.tilde, lambda = lambda,
                          theta0 = theta0, q = q)
  D.loss = obj.loss$D
  d.loss = obj.loss$d

  if (stationary) {
    eval_f <- function(eta) {
      list(objective = loss.eta(D.loss, d.loss, eta),
           gradient = as.numeric(D.loss %*% eta - d.loss))
    }
    eval_g_ineq <- function(eta) {
      list(constraints = constraint_fn(eta, lambda, q),
           jacobian = jacobian_fn(eta, lambda, q))
    }

    result = nloptr::nloptr(
      x0 = rep(0, p * (q + 1)),
      eval_f = eval_f,
      eval_g_ineq = eval_g_ineq,
      opts = list(algorithm = "NLOPT_LD_MMA",
                  xtol_rel = 1e-6,
                  maxeval = max.iter)
    )
    eta = result$solution
    loss = result$objective
    conv = result$status
  } else {
    eta = try(solve(D.loss, d.loss), silent = TRUE)
    if (inherits(eta, "try-error")) {
      return(NULL)
    }
    eta = as.numeric(eta)
    loss = loss.eta(D.loss, d.loss, eta)
    conv = TRUE
    result = NULL
  }

  bsz = q + 1
  R1 = Filter.from.eta(eta[1:bsz], eigen.L, q)
  R2 = Filter.from.eta(eta[(bsz + 1):(2 * bsz)], eigen.L, q)
  R3 = Filter.from.eta(eta[(2 * bsz + 1):(3 * bsz)], eigen.L, q)

  list(result = result, eta = eta,
       R_list = list(R1 = R1, R2 = R2, R3 = R3),
       R1 = R1, R2 = R2, R3 = R3, loss = loss, conv = conv)
}

fit.eta <- function(data, L, theta0, p, q, eps_thre = 1e-6,
                    stationary = TRUE, max.iter = 1000) {
  if (p == 1) {
    QP.TAR.eta.p1(data = data, L = L, theta0 = theta0, q = q,
                  eps = eps_thre, stationary = stationary)
  } else if (p == 2) {
    QP.TAR.eta.p2(data = data, L = L, theta0 = theta0, q = q,
                  eps = eps_thre, stationary = stationary)
  } else if (p == 3) {
    eta.nlr.mod(data = data, L = L, theta0 = theta0, p = p, q = q,
                eps_thre = eps_thre, max.iter = max.iter,
                stationary = stationary)
  } else {
    stop("p must be the AR order 1, 2, or 3")
  }
}

step.00 <- function(data, p, eps_thre = 1e-6) {
  R.initial(data = data, p = p, eps.thre = eps_thre)
}

step.0 <- function(data, R_list, p) {
  S.theta0(data = data, R_list = R_list, p = p)
}

step.2 <- function(data, L, theta0, p, q, eps_thre, stationary = TRUE,
                   eta_max_iter = 1000) {
  fit.eta(data = data, L = L, theta0 = theta0, p = p, q = q,
          eps_thre = eps_thre, stationary = stationary,
          max.iter = eta_max_iter)
}

add.R.fields <- function(x, R_list, suffix = "") {
  for (k in seq_along(R_list)) {
    x[[paste0("R", k, suffix)]] = R_list[[k]]
  }
  x
}

step.thre <- function(resList, lambda.v, net.thre) {
  A.net.e = vector(mode = "list", length = length(lambda.v))

  for (j in seq_along(lambda.v)) {
    A.net.e[[j]] = vector(mode = "list", length = length(net.thre))
    L.lambda = resList[[j]]$L * resList[[j]]$theta1

    for (k in seq_along(net.thre)) {
      temp = abs(L.lambda) > net.thre[k]
      diag(temp) = 0
      A.net.e[[j]][[k]] = temp
    }
  }

  A.net.e
}

step.lap.est <- function(data, resList, p, q, A.net.e, lambda.v, net.thre,
                         time.vec, model = "LN", eps_thre = 1e-6,
                         eps_abs = 1e-5, eps_rel = 1e-3,
                         max_iter = 50000, deg_max_iter = 50000,
                         lap_z_max_iter = max_iter,
                         eta_max_iter = 1000, verbose = FALSE,
                         stationary = TRUE) {
  results = vector(mode = "list", length = length(lambda.v))
  d = ncol(data)
  n_resid = nrow(data) - p
  ini.mat = matrix(0, d, d)
  ini.vec = rep(0, d)

  for (j in seq_along(lambda.v)) {
    results[[j]] = vector(mode = "list", length = length(net.thre))
    resList.j = resList[[j]]
    S.c = resList.j$S
    theta0.c = resList.j$theta0

    for (k in seq_along(net.thre)) {
      time.start = proc.time()[3]
      A.c = A.net.e[[j]][[k]]
      rho.refit = sqrt(log(d) / n_resid)

      refit.2 = try(ADMM_L2_Zero(
        SS = S.c, theta0 = theta0.c, v = ini.vec, rho = rho.refit,
        AA = A.c, model = model, Z_ini = ini.mat, W_ini = ini.mat,
        eps_thre = eps_thre, eps_abs = eps_abs, eps_rel = eps_rel,
        max_iter = max_iter, verbose = verbose
      ), silent = TRUE)
      conv.2 = !inherits(refit.2, "try-error")
      if (!conv.2) {
        refit.2 = NULL
      }

      v0.res = NULL
      v0.est = ini.vec
      conv.v0 = FALSE
      refit.3 = NULL
      conv.3 = FALSE

      if (conv.2 && model != "L") {
        v0.res = try(ADMM.Deg.L(
          refit.2$L, rho = 0.1, epsilon = sqrt(1 / (2 * d)),
          eps.abs = 1e-5, eps.rel = 1e-3, max.iter = deg_max_iter,
          verbose = verbose
        ), silent = TRUE)

        if (!inherits(v0.res, "try-error")) {
          v0.est = v0.res$v
          conv.v0 = v0.res$conv
        } else {
          v0.res = NULL
        }
      }

      if (conv.v0 && model != "L") {
        refit.3 = try(ADMM_Lap_Zero(
          SS = S.c, V0 = v0.est, rho = rho.refit, AA = A.c,
          model = model, ZZ_ini = ini.mat, WW_ini = ini.mat,
          phi_ini = 1e-6, eps_thre = eps_thre,
          eps_abs = eps_abs, eps_rel = eps_rel,
          max_iter = max_iter, Z_max_iter = lap_z_max_iter,
          Z_conv_abs = 1e-5, Z_conv_rel = 1e-3,
          verbose = verbose
        ), silent = TRUE)

        if (!inherits(refit.3, "try-error")) {
          conv.3 = TRUE
        } else {
          refit.3 = NULL
        }
      }

      if (conv.3) {
        L.unscaled = refit.3$L
        theta1.est = refit.3$theta1
        L.est = L.unscaled * theta1.est
        theta0.est = refit.3$theta0
      } else if (conv.2) {
        if (model != "L") {
          warning("Step 3 of refitting failed. Using Step 2 estimate for L instead.",
                  call. = FALSE)
        }
        L.unscaled = refit.2$L
        theta1.est = refit.2$theta1
        L.est = L.unscaled * theta1.est
        theta0.est = theta0.c
      } else {
        L.unscaled = NULL
        theta1.est = NA_real_
        L.est = NULL
        theta0.est = theta0.c
      }

      eta.refit = NULL
      if (!is.null(L.est)) {
        eta.refit = step.2(
          data = data, L = L.est, theta0 = theta0.est, p = p, q = q,
          eps_thre = eps_thre, stationary = stationary,
          eta_max_iter = eta_max_iter
        )
      }

      time.end = proc.time()[3]
      time.refit = time.end - time.start + time.vec[j]

      item = list(
        S = S.c,
        theta0.ini = theta0.c,
        theta0 = theta0.est,
        theta1 = theta1.est,
        L = L.unscaled,
        L.scaled = L.est,
        result.0.post = refit.2,
        conv.2 = conv.2,
        theta0.0S = if (conv.3) refit.3$theta0 else NA,
        result.0S = refit.3,
        v0.est = v0.est,
        conv.3 = conv.3,
        A.net = A.c,
        R_list = if (is.null(eta.refit)) NULL else eta.refit$R_list,
        eta = if (is.null(eta.refit)) NULL else eta.refit$eta,
        eta.0S = if (is.null(eta.refit)) NULL else eta.refit$eta,
        time = time.refit
      )
      if (!is.null(eta.refit)) {
        item = add.R.fields(item, eta.refit$R_list, suffix = ".0S")
      }
      results[[j]][[k]] = item
    }
  }

  results
}

validate_targar_fit_inputs <- function(data, p, q, lambda.v, net.thre, model,
                                       rho.v, num.pass, eps.thre, eps_abs,
                                       eps_rel, max_iter, deg_max_iter,
                                       lap_z_max_iter, eta_max_iter,
                                       stationary, verbose) {
  if (!is.matrix(data) || !is.numeric(data) || any(!is.finite(data))) {
    stop("`data` must be a finite numeric matrix.")
  }
  if (!is.numeric(p) || length(p) != 1 || !is.finite(p) ||
      p != as.integer(p) || !p %in% 1:3) {
    stop("`p` must be the AR order 1, 2, or 3.")
  }
  if (nrow(data) <= p) {
    stop("`data` must have more rows than the AR order `p`.")
  }
  if (!is.numeric(q) || length(q) != 1 || !is.finite(q) ||
      q != as.integer(q) || q < 1) {
    stop("`q` must be a positive integer.")
  }
  if (!is.numeric(lambda.v) || length(lambda.v) == 0 ||
      any(!is.finite(lambda.v)) || any(lambda.v <= 0)) {
    stop("`lambda.v` must be a non-empty numeric vector of positive values.")
  }
  if (!is.numeric(rho.v) || length(rho.v) != length(lambda.v) ||
      any(!is.finite(rho.v)) || any(rho.v <= 0)) {
    stop("`rho.v` must be a positive numeric vector with the same length as `lambda.v`.")
  }
  if (!is.numeric(net.thre) || length(net.thre) == 0 ||
      any(!is.finite(net.thre)) || any(net.thre < 0)) {
    stop("`net.thre` must be a non-empty numeric vector of non-negative values.")
  }
  if (!is.numeric(num.pass) || length(num.pass) != 1 ||
      !is.finite(num.pass) || num.pass <= 0 ||
      num.pass != as.integer(num.pass)) {
    stop("`num.pass` must be a positive integer scalar.")
  }
  if (any(!is.finite(c(eps.thre, eps_abs, eps_rel))) ||
      eps.thre <= 0 || eps_abs <= 0 || eps_rel <= 0) {
    stop("`eps.thre`, `eps_abs`, and `eps_rel` must be positive finite scalars.")
  }
  iter.values = c(max_iter, deg_max_iter, lap_z_max_iter, eta_max_iter)
  if (any(!is.finite(iter.values)) || any(iter.values <= 0) ||
      any(iter.values != as.integer(iter.values))) {
    stop("Iteration limits must be positive integer scalars.")
  }
  if (!is.logical(stationary) || length(stationary) != 1 || is.na(stationary)) {
    stop("`stationary` must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  invisible(model)
}

TAR.GAR.fit <- function(data, p, q = 1, lambda.v, net.thre,
                        rho.v = lambda.v, num.pass = 2, model = "LN",
                        eps.thre = 1e-6, eps_abs = 1e-5,
                        eps_rel = 1e-3, max_iter = 50000,
                        deg_max_iter = 50000,
                        lap_z_max_iter = max_iter,
                        eta_max_iter = 1000,
                        stationary = TRUE, verbose = FALSE,
                        max.iter.eta = NULL) {
  if (!is.null(max.iter.eta)) {
    eta_max_iter = max.iter.eta
  }

  data = as.matrix(data)
  model.input = model
  model.internal = normalize_gar_model(model)
  validate_targar_fit_inputs(
    data = data, p = p, q = q, lambda.v = lambda.v, net.thre = net.thre,
    model = model.internal, rho.v = rho.v, num.pass = num.pass,
    eps.thre = eps.thre, eps_abs = eps_abs, eps_rel = eps_rel,
    max_iter = max_iter, deg_max_iter = deg_max_iter,
    lap_z_max_iter = lap_z_max_iter, eta_max_iter = eta_max_iter,
    stationary = stationary, verbose = verbose
  )

  p = as.integer(p)
  q = as.integer(q)
  num.pass = as.integer(num.pass)
  max_iter = as.integer(max_iter)
  deg_max_iter = as.integer(deg_max_iter)
  lap_z_max_iter = as.integer(lap_z_max_iter)
  eta_max_iter = as.integer(eta_max_iter)

  eps_thre = eps.thre
  time.tot.beg = proc.time()[3]
  resList = vector(length = length(lambda.v), mode = "list")
  iniRes = resList
  time.vec = numeric(length(lambda.v))
  d = ncol(data)

  R.list = step.00(data = data, p = p, eps_thre = eps_thre)
  R.ini = R.list
  if (verbose) {
    message("Step 00 complete")
  }

  for (j in seq_along(lambda.v)) {
    if (verbose) {
      message(paste("Starting fitting for lambda =", lambda.v[j]))
    }
    time.start = proc.time()[3]
    eta = NULL
    fit.L.theta1 = NULL
    S.i = NULL
    theta0.i = NULL
    conv.step.1.i = FALSE

    for (i in 1:num.pass) {
      fit.S.theta0 = step.0(data = data, R_list = R.list, p = p)
      S.i = fit.S.theta0$S
      theta0.i = fit.S.theta0$theta0
      if (verbose) {
        message("step 0 complete")
      }

      v0 = if (model.internal == "L") rep(1, d) else rep(0, d)
      fit.L.theta1 = try(ADMM_L2(
        s = S.i, theta0 = theta0.i, v = v0, rho = rho.v[j],
        lambda = lambda.v[j], model = model.internal,
        Z_ini = matrix(0, d, d), W_ini = matrix(0, d, d),
        eps_thre = eps_thre, eps_abs = eps_abs, eps_rel = eps_rel,
        max_iter = max_iter, verbose = verbose
      ), silent = TRUE)
      if (inherits(fit.L.theta1, "try-error")) {
        stop(sprintf("ADMM_L2 failed at lambda index %d, pass %d", j, i))
      }
      L.i = fit.L.theta1$L * fit.L.theta1$theta1
      conv.step.1.i = fit.L.theta1$conv
      if (verbose) {
        message("step 1 complete")
      }

      fit.eta = step.2(
        data = data, L = L.i, theta0 = theta0.i, p = p, q = q,
        eps_thre = eps_thre, stationary = stationary,
        eta_max_iter = eta_max_iter
      )
      if (is.null(fit.eta)) {
        stop(sprintf("Eta fit failed at lambda index %d, pass %d", j, i))
      }
      eta = fit.eta$eta
      R.list = fit.eta$R_list
      if (verbose) {
        message("step 2/3 complete")
      }

      if (i == 1) {
        ini.item = list(
          R_ini = R.ini,
          L.ini = fit.L.theta1$L,
          theta1 = fit.L.theta1$theta1,
          eta.ini = eta,
          S.ini = S.i,
          theta0.ini = theta0.i,
          conv.step1 = conv.step.1.i
        )
        ini.item = add.R.fields(ini.item, R.ini, suffix = ".ini")
        iniRes[[j]] = ini.item
      }
    }

    time.end = proc.time()[3]
    time.vec[j] = time.end - time.start

    res.item = list(
      S = S.i,
      theta0 = theta0.i,
      L = fit.L.theta1$L,
      theta1 = fit.L.theta1$theta1,
      eta = eta,
      R_list = R.list,
      conv.step1 = conv.step.1.i
    )
    res.item = add.R.fields(res.item, R.list)
    resList[[j]] = res.item
    if (verbose) {
      message(paste("Pass", i, "complete"))
    }
    R.list = R.ini
  }

  if (verbose) {
    message("Starting refitting")
  }
  refit.zero = step.thre(resList, lambda.v, net.thre)
  refit.L = step.lap.est(
    data = data, resList = resList, p = p, q = q,
    A.net.e = refit.zero, lambda.v = lambda.v, net.thre = net.thre,
    time.vec = time.vec, model = model.internal, eps_thre = eps_thre,
    eps_abs = eps_abs, eps_rel = eps_rel, max_iter = max_iter,
    deg_max_iter = deg_max_iter, lap_z_max_iter = lap_z_max_iter,
    eta_max_iter = eta_max_iter, verbose = verbose,
    stationary = stationary
  )

  time.tot.end = proc.time()[3]
  time.total = time.tot.end - time.tot.beg

  result = list(result.pass = resList, ini = iniRes, refit = refit.L,
                A.net = refit.zero, time.total = time.total,
                nobs = nrow(data), n.resid = nrow(data) - p,
                p = p, q = q, lambda.v = lambda.v, rho.v = rho.v,
                net.thre = net.thre, model = model.input,
                model.internal = model.internal, fit.type = "TARGAR",
                stationary = stationary)
  class(result) = c("TARGAR_fit", "list")
  result
}

#' TAR-GAR fitting procedure
#'
#' @description
#' `TARGAR_fit()` fits a TAR-GAR(p,q) model to temporally dependent
#' multivariate data. It alternates between estimating AR filter matrices
#' `R1`, ..., `Rp` and estimating the latent graph Laplacian for the
#' innovation process, then refits the graph over the supplied `net.thre`
#' sequence.
#'
#' @param data An n by d numeric data matrix, with rows ordered in time.
#' @param lambda.v Positive tuning parameter sequence controlling graph
#'   sparsity.
#' @param net.thre Non-negative threshold sequence used to define graph
#'   zero-patterns for the final refit.
#' @param rho.v Positive ADMM parameter sequence. Defaults to `lambda.v`.
#' @param num.pass Positive integer number of alternating TAR-GAR passes before
#'   the final refit. Defaults to 2.
#' @param model One of `"LN"`, `"L"`, or `"LN.noselfloop"`, matching the GAR
#'   structure used for the innovation covariance.
#' @param p AR order. Must be 1, 2, or 3. Defaults to 1 for compatibility with
#'   the original package API.
#' @param q Polynomial order in the graph Laplacian. Must be a positive
#'   integer. Defaults to 1.
#' @param eps.thre Small positive numerical threshold. Defaults to 1e-6.
#' @param eps_abs Absolute ADMM convergence tolerance. Defaults to 1e-5.
#' @param eps_rel Relative ADMM convergence tolerance. Defaults to 1e-3.
#' @param max_iter Maximum number of ADMM iterations. Defaults to 50000.
#' @param deg_max_iter Maximum number of iterations for the degree-vector ADMM
#'   refit. Defaults to 50000.
#' @param lap_z_max_iter Maximum number of inner Z updates in the Laplacian
#'   refit. Defaults to `max_iter`.
#' @param eta_max_iter Maximum number of eta optimization iterations for the
#'   p = 3 nonlinear program. Defaults to 1000.
#' @param stationary Logical flag indicating whether to impose stationarity
#'   constraints while fitting the TAR filters. Defaults to TRUE.
#' @param verbose Logical flag indicating whether to print progress.
#' @param max.iter.eta Backward-compatible alias for `eta_max_iter`.
#'
#' @returns A list with class `"TARGAR_fit"` containing:
#' * `result.pass`: pre-refit results indexed by `lambda.v`.
#' * `ini`: first-pass initial estimates indexed by `lambda.v`.
#' * `refit`: final refit results indexed by `lambda.v` and `net.thre`.
#' * `A.net`: thresholded graph zero-patterns indexed by `lambda.v` and
#'   `net.thre`.
#' * `p`, `q`, `lambda.v`, `rho.v`, `net.thre`, `model`, and timing metadata.
#'
#' `model_selec()` can be used on the returned object to select a final model
#' with the same workflow used for `GAR1_fit()`.
#'
#' @examples
#' \dontrun{
#' library(SGM)
#' data("targar")
#'
#' y = targar$data
#' n = nrow(y)
#' d = ncol(y)
#'
#' lambda.v = c(1, 0.5) * sqrt(log(d) / (n - 1))
#' rho.v = pmax(lambda.v, 0.01)
#' net.thre = exp(seq(log(1), log(0.1), length.out = 5)) *
#'   sqrt(log(d) / (n - 1))
#'
#' fit = TARGAR_fit(
#'   y, lambda.v = lambda.v, net.thre = net.thre, rho.v = rho.v,
#'   p = 1, q = 1, num.pass = 2, model = "LN"
#' )
#' sel = model_selec(fit)
#'
#' sel$theta0
#' sel$theta1
#' sel$L
#' sel$R_list
#' }
#'
#' @export
TARGAR_fit <- function(data, lambda.v, net.thre, rho.v = lambda.v,
                       num.pass = 2, model = "LN", p = 1, q = 1,
                       eps.thre = 1e-6, eps_abs = 1e-5,
                       eps_rel = 1e-3, max_iter = 50000,
                       deg_max_iter = 50000,
                       lap_z_max_iter = max_iter,
                       eta_max_iter = 1000, stationary = TRUE,
                       verbose = FALSE, max.iter.eta = NULL) {
  TAR.GAR.fit(
    data = data, p = p, q = q, lambda.v = lambda.v, net.thre = net.thre,
    rho.v = rho.v, num.pass = num.pass, model = model,
    eps.thre = eps.thre, eps_abs = eps_abs, eps_rel = eps_rel,
    max_iter = max_iter, deg_max_iter = deg_max_iter,
    lap_z_max_iter = lap_z_max_iter, eta_max_iter = eta_max_iter,
    stationary = stationary, verbose = verbose,
    max.iter.eta = max.iter.eta
  )
}

#' @rdname TARGAR_fit
#' @export
fit_TAR_GAR <- function(data, p, q, lambda.v, net.thre, rho.v = lambda.v,
                        num.pass = 2, model = "LN", eps.thre = 1e-6,
                        eps_abs = 1e-5, eps_rel = 1e-3,
                        max_iter = 50000, deg_max_iter = 50000,
                        lap_z_max_iter = max_iter,
                        eta_max_iter = 1000, stationary = TRUE,
                        verbose = FALSE, max.iter.eta = NULL) {
  TARGAR_fit(
    data = data, lambda.v = lambda.v, net.thre = net.thre, rho.v = rho.v,
    num.pass = num.pass, model = model, p = p, q = q,
    eps.thre = eps.thre, eps_abs = eps_abs, eps_rel = eps_rel,
    max_iter = max_iter, deg_max_iter = deg_max_iter,
    lap_z_max_iter = lap_z_max_iter, eta_max_iter = eta_max_iter,
    stationary = stationary, verbose = verbose,
    max.iter.eta = max.iter.eta
  )
}
