## Goal: Create functions to apply GAR(1) fitting procedure, and model selection
## 2024/09/12

####################
## Required Packages
####################

#require(SGM) # R Package that has our C++ functions
#require(gmp) # R package used for handling numerical issues
#require(foreach) # For running in Parallel
#require(doParallel) # For running in Parallel

count_edges <- function(A) {
  if (is.null(A)) {
    return(NA_real_)
  }
  sum(A) / 2
}

normalize_gar_model <- function(model) {
  if (!is.character(model) || length(model) != 1 || is.na(model)) {
    stop("`model` must be a single character string.")
  }
  if (identical(model, "LN.noselfloop")) {
    return("LN.noloop")
  }
  if (model %in% c("LN", "L", "LN.noloop")) {
    return(model)
  }
  stop("`model` must be one of \"LN\", \"L\", or \"LN.noselfloop\".")
}

validate_gar1_fit_inputs <- function(S, nobs, lambda.v, net.thre, model, step, rho.v, eps_thre, eps_abs, eps_rel, max_iter_s1, max_iter_s2, max_iter_s3) {
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) != ncol(S) || any(!is.finite(S))) {
    stop("`S` must be a finite numeric square matrix.")
  }
  if (!isTRUE(all.equal(S, t(S), check.attributes = FALSE, tolerance = 1e-8))) {
    stop("`S` must be symmetric.")
  }
  if (!is.numeric(nobs) || length(nobs) != 1 || !is.finite(nobs) || nobs <= 0 || nobs != as.integer(nobs)) {
    stop("`nobs` must be a positive integer scalar.")
  }
  if (!is.numeric(lambda.v) || length(lambda.v) == 0 || any(!is.finite(lambda.v)) || any(lambda.v <= 0)) {
    stop("`lambda.v` must be a non-empty numeric vector of positive values.")
  }
  if (!is.numeric(rho.v) || length(rho.v) != length(lambda.v) || any(!is.finite(rho.v)) || any(rho.v <= 0)) {
    stop("`rho.v` must be a positive numeric vector with the same length as `lambda.v`.")
  }
  if (!is.numeric(net.thre) || length(net.thre) == 0 || any(!is.finite(net.thre)) || any(net.thre < 0)) {
    stop("`net.thre` must be a non-empty numeric vector of non-negative values.")
  }
  if (!is.numeric(step) || length(step) != 1 || !step %in% 1:3) {
    stop("`step` must be one of 1, 2, or 3.")
  }
  if (any(!is.finite(c(eps_thre, eps_abs, eps_rel))) || eps_thre <= 0 || eps_abs <= 0 || eps_rel <= 0) {
    stop("`eps_thre`, `eps_abs`, and `eps_rel` must be positive finite scalars.")
  }
  iter.values = c(max_iter_s1, max_iter_s2, max_iter_s3)
  if (any(!is.finite(iter.values)) || any(iter.values <= 0) || any(iter.values != as.integer(iter.values))) {
    stop("`max_iter_s1`, `max_iter_s2`, and `max_iter_s3` must be positive integer scalars.")
  }
  invisible(model)
}

best_finite_index <- function(x, label) {
  finite.idx = which(is.finite(x), arr.ind = TRUE)
  if (length(finite.idx) == 0) {
    stop(sprintf("No finite eBIC values were available for %s. Check convergence diagnostics before calling `model_selec()`.", label))
  }
  values = x[finite.idx]
  finite.idx[which.min(values), , drop = FALSE]
}


##############
## Module 1: 3-Step Estimation Procedure for GAR(1) Model
########

################
### Step 0a: Get sample covariance matrix and theta0 estimate
##########

fit_step_0a = function(S, nobs){
  n = nobs

  sigma = eigen(S, symmetric = TRUE, only.values = TRUE)
  theta0.e = sqrt(1 / max(sigma$values))

  temp = list("S" = S, "theta0" = theta0.e, "n" = n)
  return(temp)
}

###############
### Step 1: fit L given theta_0.e from Step 0 and obtain thresholded zero-patterns
#############

fit_step_1 = function(step0a, lambda.v, rho.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_s1, verbose){
  S = step0a$S
  theta0.e = step0a$theta0

  p = ncol(S)
  Z = matrix(0, p, p)
  W = Z

  step1.fit = vector("list", length(lambda.v))
  A.net = vector("list", length(lambda.v))

  for (j in seq_along(lambda.v)){
    step1.fit[[j]] = ADMM_L2(S, theta0.e, rep(0, p), rho.v[j], lambda.v[j], model, Z, W, eps_thre, eps_abs, eps_rel, max_iter_s1, verbose)

    if (!is.null(step1.fit[[j]]) && isTRUE(step1.fit[[j]]$conv)) {
      A.net[[j]] = vector("list", length(net.thre))

      for (k in seq_along(net.thre)){
        net.e = abs(step1.fit[[j]]$L) > net.thre[k]
        diag(net.e) = 0
        A.net[[j]][[k]] = net.e
      }
    }
  }

  conv.step1 = vapply(step1.fit, function(result.c) {
    !is.null(result.c) && isTRUE(result.c$conv)
  }, logical(1))

  if (!any(conv.step1)) {
    warning("Step 1 did not converge for any value of `lambda.v`.")
  }

  return(list("fit" = step1.fit, "A.net" = A.net, "conv" = conv.step1))
}

#################
### Step 2: refit L given the pattern from Step 1
############

fit_step_2 = function(step0a, step1, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_s2, verbose){
  S = step0a$S
  theta0.e = step0a$theta0
  p = ncol(S)
  n = step0a$n

  step2.fit = vector("list", length(lambda.v))
  Z = matrix(0, p, p)
  W = Z

  A.net = step1$A.net
  conv.step1 = step1$conv
  conv.step2 = matrix(FALSE, length(lambda.v), length(net.thre))

  for (j in seq_along(lambda.v)){
    if (isTRUE(conv.step1[j])) {
      step2.fit[[j]] = vector("list", length(net.thre))
      for (k in seq_along(net.thre)){
        step2.fit[[j]][[k]] = ADMM_L2_Zero(S, theta0.e, v = rep(0, p), rho = sqrt(log(p) / n), A = A.net[[j]][[k]], model, Z_ini = Z, W_ini = W,
                                           eps_thre, eps_abs, eps_rel, max_iter_s2, verbose)
        conv.step2[j, k] = !is.null(step2.fit[[j]][[k]]) && isTRUE(step2.fit[[j]][[k]]$conv)
      }
    }
  }

  if (any(conv.step1) && !any(conv.step2)) {
    warning("Step 2 did not converge for any combination of `lambda.v` and `net.thre`.")
  }

  return(list("fit" = step2.fit, "conv" = conv.step2))
}

######################
## Step 3 helpers
#############

fit_step_3a = function(step0a, step2, lambda.v, net.thre, eps_abs, eps_rel, verbose){
  v0.s3 = vector("list", length(lambda.v))
  conv.step3a = matrix(FALSE, length(lambda.v), length(net.thre))

  step2.fit = step2$fit
  p = ncol(step0a$S)

  for (j in seq_along(lambda.v)){
    v0.s3[[j]] = vector("list", length(net.thre))

    for (k in seq_along(net.thre)){
      result.c = step2.fit[[j]][[k]]

      if (!is.null(result.c) && isTRUE(result.c$conv)) {
        L.est = result.c$L
        temp = try(ADMM.Deg.L(L.est, rho = 0.1, epsilon = sqrt(1 / (2 * ncol(L.est))), eps.abs = 1e-5, eps.rel = 1e-3, max.iter = 100000, verbose = FALSE))

        if (inherits(temp, "try-error")){
          temp = list("v" = rep(0, p), "conv" = FALSE)
        }

        if (!isTRUE(temp$conv)){
          warning(sprintf("v0 did not converge at (lambda.index=%d, net.thre.index=%d). Replacing with zero vector.", j, k))
        }

        v0.s3[[j]][[k]] = temp$v
        conv.step3a[j, k] = isTRUE(temp$conv)
      }
    }
  }

  return(list("v0.s3" = v0.s3, "conv" = conv.step3a))
}

## Estimate theta0 and L simultaneously with 0-pattern from Step 1 and estimated v0 from Step 3a
fit_step_3b = function(step0a, step1, step2, step3a, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_s3, verbose){
  S = step0a$S
  n = step0a$n
  p = ncol(S)

  A.net = step1$A.net
  v0.s3 = step3a$v0.s3
  conv.step3a = step3a$conv

  step3b.fit = vector("list", length(lambda.v))
  theta0.s3 = matrix(NA, length(lambda.v), length(net.thre))
  conv.step3b = matrix(FALSE, length(lambda.v), length(net.thre))
  Z = matrix(0, p, p)
  W = Z
  phi = 0

  for (j in seq_along(lambda.v)){
    step3b.fit[[j]] = vector("list", length(net.thre))

    for (k in seq_along(net.thre)){
      if (isTRUE(conv.step3a[j, k])) {
        v0.e = v0.s3[[j]][[k]]
        step3b.fit[[j]][[k]] = ADMM_Lap_Zero(S, v0.e, rho = sqrt(log(p) / n), AA = A.net[[j]][[k]], model = model, ZZ_ini = Z, WW_ini = W, phi_ini = phi,
                                             eps_thre = eps_thre, eps_abs = eps_abs, eps_rel = eps_rel, max_iter = max_iter_s3, Z_max_iter = 100000,
                                             Z_conv_abs = 1e-5, Z_conv_rel = 1e-3, verbose = verbose)
        theta0.s3[j, k] = step3b.fit[[j]][[k]]$theta0
        conv.step3b[j, k] = !is.null(step3b.fit[[j]][[k]]) && isTRUE(step3b.fit[[j]][[k]]$conv)
      }
    }
  }

  if (any(conv.step3a) && !any(conv.step3b)) {
    warning("Step 3b did not converge for any combination of `lambda.v` and `net.thre`.")
  }

  return(list("fit" = step3b.fit, "theta0.s3" = theta0.s3, "conv" = conv.step3b))
}

fit_step_3 = function(step0a, step1, step2, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_s3, verbose){
  step_3a = fit_step_3a(step0a, step2, lambda.v, net.thre, eps_abs, eps_rel, verbose)
  step_3b = fit_step_3b(step0a, step1, step2, step_3a, lambda.v, net.thre, model, eps_thre, eps_abs, eps_rel, max_iter_s3, verbose)

  return(list(
    "step3a" = step_3a$v0.s3,
    "step3b" = step_3b$fit,
    "theta0.s3" = step_3b$theta0.s3,
    "v0.s3" = step_3a$v0.s3,
    "conv" = list("step3a" = step_3a$conv, "step3b" = step_3b$conv)
  ))
}


#################################
## Model Fitting function
################################
#' GAR(1) fitting procedure
#'
#' @description
#' `GAR1_fit` performs a three-step estimation procedure, using a penalized MLE approach, to estimate graph filter parameters `theta0` and `theta`, and the normalized graph Laplacian `L`.
#'
#' @param S An estimate of the covariance matrix, such as the MLE.
#' @param nobs The number of samples used to calculate `S`
#' @param lambda.v Tuning parameter to control sparsity of the estimated graph
#' @param net.thre Tuning parameter to control noisy entries in estimated graph
#' @param model
#' * `"LN"`: Fits a normalized graph Laplacian
#' * `"L"`: Fits a graph Laplacian
#' * `"LN.noselfloop"`: Fits a normalized graph laplacian assuming no self-loops.
#' @param step How many steps of the estimation procedure you want to run. Either 1, 2, or 3.
#' @param rho.v ADMM parameter (typically equal to `lambda.v`)
#' @param eps_thre Small positive number
#' @param eps_abs ADMM convergence criterion
#' @param eps_rel ADMM convergence criterion
#' @param max_iter_s1 Maximum number of iterations for Step 1.
#' @param max_iter_s2 Maximum number of iterations for Step 2.
#' @param max_iter_s3 Maximum number of iterations for Step 3b.
#'
#' @returns
#' A list object
#' * `S` The supplied covariance estimate.
#' * `nobs` The number of observations used to form `S`.
#' * `model` The fitted GAR model family.
#' * `step` The last step requested in the fitting procedure.
#' * `lambda.v` The sparsity tuning parameter sequence.
#' * `rho.v` The ADMM tuning parameter sequence.
#' * `net.thre` The graph-threshold sequence used in Steps 2 and 3.
#' * `theta0.init` The initial estimate of `theta0` from Step 0, used in Steps 1 and 2.
#' * `theta0.s3` A matrix of Step 3b estimates for `theta0` (NULL if `step<3`).
#' * `A.net` A list of adjacency matrices defining the zero-patterns created in Step 1 and used in Steps 2 and 3.
#' * `step1` A list containing the Step 1 fits, indexed by `lambda.v`.
#' * `step2` A list containing the Step 2 fits, indexed by `lambda.v` and `net.thre` (NULL if `step<2`).
#' * `step3a` A list containing the Step 3a `v0` estimates, indexed by `lambda.v` and `net.thre` (NULL if `step<3`).
#' * `step3b` A list containing the Step 3b joint fits, indexed by `lambda.v` and `net.thre` (NULL if `step<3`).
#' * `v0.s3` A list containing the Step 3a `v0` estimates (NULL if `step<3`).
#' * `conv` A list of convergence diagnostics with components `step1`, `step2`, `step3a`, and `step3b`.
#'
#' @example man-roxygen/GAR1_fit_example.R
#' @export
GAR1_fit = function(S, nobs, lambda.v, net.thre, model = "LN", step = 3, rho.v = lambda.v, eps_thre = 1e-6, eps_abs = 1e-5, eps_rel = 1e-3, max_iter_s1 = 10000, max_iter_s2 = 10000, max_iter_s3 = 10000, verbose = FALSE){
  model.input = model
  model.internal = normalize_gar_model(model)
  validate_gar1_fit_inputs(S, nobs, lambda.v, net.thre, model.internal, step, rho.v, eps_thre, eps_abs, eps_rel, max_iter_s1, max_iter_s2, max_iter_s3)

  step0a = fit_step_0a(S, nobs)
  if (verbose) {
    print("Step 0a complete")
  }

  step1 = fit_step_1(step0a, lambda.v, rho.v, net.thre, model.internal, eps_thre, eps_abs, eps_rel, max_iter_s1, verbose)
  if (verbose) {
    print("Step 1 complete")
  }

  if (step >= 2){
    step2 = fit_step_2(step0a, step1, lambda.v, net.thre, model.internal, eps_thre, eps_abs, eps_rel, max_iter_s2, verbose)
    if (verbose) {
      print("Step 2 complete")
    }
  } else {
    step2 = NULL
  }

  if (step == 3){
    step3 = fit_step_3(step0a, step1, step2, lambda.v, net.thre, model.internal, eps_thre, eps_abs, eps_rel, max_iter_s3, verbose)
    if (verbose) {
      print("Step 3 complete")
    }
  } else {
    step3 = NULL
  }

  resultList = list(
    "S" = step0a$S,
    "nobs" = nobs,
    "model" = model.input,
    "step" = step,
    "lambda.v" = lambda.v,
    "rho.v" = rho.v,
    "net.thre" = net.thre,
    "theta0.init" = step0a$theta0,
    "theta0.s3" = if (is.null(step3)) NULL else step3$theta0.s3,
    "A.net" = step1$A.net,
    "step1" = step1$fit,
    "step2" = if (is.null(step2)) NULL else step2$fit,
    "step3a" = if (is.null(step3)) NULL else step3$step3a,
    "step3b" = if (is.null(step3)) NULL else step3$step3b,
    "v0.s3" = if (is.null(step3)) NULL else step3$v0.s3,
    "conv" = list(
      "step1" = step1$conv,
      "step2" = if (is.null(step2)) NULL else step2$conv,
      "step3a" = if (is.null(step3)) NULL else step3$conv$step3a,
      "step3b" = if (is.null(step3)) NULL else step3$conv$step3b
    )
  )

  return(resultList)
}


##############
## Module 2: Fitting Process for GAR(1) Model
## Goal: Select model via eBIC and provide bootstrap goodness of fit
########

####################################
## eBIC and log-likelihood selection
######
#' Select tuning parameters for GAR(1) model
#'
#' @description
#' Given a fitted GAR(1) model path from `GAR1_fit`, uses the eBIC criterion to select the appropriate tuning parameters.
#'
#' @param resultList A list output from `GAR1_fit`.
#'
#' @returns A list object
#' * `selected.model`: A list containing the ADMM output for the selected model.
#' * `theta0`: The selected `theta0` estimate.
#' * `theta1`: The selected `theta1` estimate.
#' * `L`: The selected graph Laplacian estimate.
#' * `v0`: The selected `v0` estimate (NULL when `step=2`).
#' * `A.net.e`: A matrix encoding the (unweighted) graph chosen by the eBIC criterion.
#' * `index`: Index for the optimal tuning parameters (lambda, net.thre) for the eBIC-selected model.
#' * `lambda.v`: The selected `lambda.v` value.
#' * `net.thre`: The selected `net.thre` value.
#' * `ebic`: ebic score for the selected model.
#'
#' @export
model_selec = function(resultList){
  if (!is.list(resultList)) {
    stop("`resultList` must be a fitted model object.")
  }

  n = if (!is.null(resultList$nobs)) resultList$nobs else resultList$n
  step = if (!is.null(resultList$step)) resultList$step else 3
  model = if (!is.null(resultList$model)) resultList$model else "TARGAR"

  if (model %in% c("LN", "L", "LN.noselfloop", "LN.noloop")){
    if (is.null(resultList$S) || is.null(resultList$nobs)) {
      stop("`resultList` must be an object returned by `GAR1_fit()`.")
    }
    if (step < 2) {
      stop("`model_selec()` requires a `GAR1_fit()` object fitted with `step >= 2`.")
    }
    if (!is.numeric(resultList$nobs) || length(resultList$nobs) != 1 || !is.finite(resultList$nobs) || resultList$nobs <= 0) {
      stop("`resultList$nobs` must be a positive finite scalar.")
    }
    if (is.null(resultList$lambda.v) || is.null(resultList$net.thre) || is.null(resultList$A.net)) {
      stop("`resultList` is missing required tuning metadata or zero-patterns.")
    }

    S = resultList$S
    p = nrow(S)

    theta0.ini = resultList$theta0.init
    A.net = resultList$A.net
    step2.fit = resultList$step2
    step3b.fit = resultList$step3b
    theta0.s3 = resultList$theta0.s3
    v0.s3 = resultList$v0.s3

    n.lambda = length(resultList$lambda.v)
    n.net.thr = length(resultList$net.thre)
    if (n.lambda == 0 || n.net.thr == 0) {
      stop("`resultList` does not contain any tuning-parameter combinations to select from.")
    }

    if (p / n > 0.5){
      gamma = 1
    } else {
      gamma = 0.5
    }

    P.total = p * (p - 1) / 2

    conv.post = matrix(NA, nrow = n.lambda, ncol = n.net.thr)
    conv.0S = conv.post

    log.post.like = conv.post
    log.0S.like = conv.post

    bic.post = conv.post
    bic.0S = conv.post

    ebic.post = conv.post
    ebic.0S = conv.post

    for (j in seq_len(n.lambda)){
      for (k in seq_len(n.net.thr)){
        if (is.null(A.net[[j]]) || is.null(A.net[[j]][[k]])) {
          next
        }
        A.net.c = A.net[[j]][[k]]
        net.size.c = count_edges(A.net.c)
        if (!is.finite(net.size.c)) {
          next
        }
        ebic.term = 2 * gamma * (lfactorial(P.total) - lfactorial(net.size.c) - lfactorial(P.total - net.size.c))

        if (step == 2){
          result.c = step2.fit[[j]][[k]]
          if (!is.null(result.c) && (conv.post[j, k] = result.c$conv) == TRUE){
            L.est = result.c$L
            theta1.e = result.c$theta1
            log.post.like[j, k] = LogLike(S, theta0.ini, theta1.e, L.est, n)

            if (model == "LN" || model == "LN.noloop"){
              bic.post[j, k] = BIC(log.post.like[j, k], n, net.size.c + 1 + p)
            } else {
              bic.post[j, k] = BIC(log.post.like[j, k], n, net.size.c + 1)
            }
            ebic.post[j, k] = bic.post[j, k] + ebic.term
          }
        }

        if (step == 3){
          result.c = step3b.fit[[j]][[k]]
          if (!is.null(result.c) && (conv.0S[j, k] = result.c$conv) == TRUE){
            theta0.e = theta0.s3[j, k]
            L.est = result.c$L
            theta1.e = result.c$theta1
            log.0S.like[j, k] = LogLike(S, theta0.e, theta1.e, L.est, n)
            bic.0S[j, k] = BIC(log.0S.like[j, k], n, net.size.c + 1)
            ebic.0S[j, k] = bic.0S[j, k] + ebic.term
          }
        }
      }
    }

    v0.opt = NULL
    if (step == 2){
      index.c = best_finite_index(ebic.post, "Step 2")

      resultOptimal = step2.fit[[index.c[1]]][[index.c[2]]]
      A.net.opt = A.net[[index.c[1]]][[index.c[2]]]
      ebic.opt = ebic.post[index.c]
      theta0.opt = theta0.ini
    }

    if (step == 3){
      index.c = best_finite_index(ebic.0S, "Step 3")
      resultOptimal = step3b.fit[[index.c[1]]][[index.c[2]]]
      A.net.opt = A.net[[index.c[1]]][[index.c[2]]]
      v0.opt = v0.s3[[index.c[1]]][[index.c[2]]]
      ebic.opt = ebic.0S[index.c]
      theta0.opt = theta0.s3[index.c[1], index.c[2]]
    }

    if (is.null(resultOptimal) || !is.list(resultOptimal) || is.null(resultOptimal$L) || is.null(resultOptimal$theta1)) {
      stop("The selected GAR model is incomplete. Check convergence diagnostics before calling `model_selec()`.")
    }
    resultOptimal$theta0 = theta0.opt
    resultOptimal$v0 = v0.opt
  } else {
    if (is.null(n) || is.null(resultList$refit)) {
      stop("`resultList` does not contain the metadata needed for model selection.")
    }
    p = nrow(resultList$refit[[1]][[1]]$A.net)
    n.lambda.v = length(resultList$refit)
    n.net.thre = length(resultList$refit[[1]])
    loglike.0S = matrix(NA, nrow = n.lambda.v, ncol = n.net.thre)
    bic.0S = loglike.0S
    ebic.0S = loglike.0S

    if (p / n > 0.5){
      gamma = 1
    } else {
      gamma = 0.5
    }

    P.total = p * (p - 1) / 2

    for (j in seq_len(n.lambda.v)){
      for (k in seq_len(n.net.thre)){
        result.c = resultList$refit[[j]][[k]]

        S.c = result.c$S
        A.c = result.c$A.net
        net.size.c = sum(A.c) / 2
        L.est = result.c$result.0S$L
        theta0.est = result.c$result.0S$theta0
        theta1.est = result.c$result.0S$theta1

        ebic.term = 2 * gamma * (lfactorial(P.total) - lfactorial(net.size.c) - lfactorial(P.total - net.size.c))

        loglike.0S[j, k] = LogLike(S = S.c, theta0 = theta0.est, theta1 = theta1.est, L.est, n - 1)
        bic.0S[j, k] = BIC(loglike.0S[j, k], n - 1, net.size.c + 3)
        ebic.0S[j, k] = bic.0S[j, k] + ebic.term
      }
    }

    index.c = which(ebic.0S == min(ebic.0S, na.rm = TRUE), arr.ind = TRUE)
    ebic.opt = ebic.0S[index.c[1], index.c[2]]
    resultOptimal = resultList$refit[[index.c[1]]][[index.c[2]]]
    A.net.opt = resultOptimal$A.net
    theta0.opt = resultOptimal$result.0S$theta0
    v0.opt = resultOptimal$v0.est
    resultOptimal$theta0 = theta0.opt
  }

  colnames(index.c) = c("lambda", "net.thre")
  rownames(index.c) = "index"

  retList = list(
    "selected.model" = resultOptimal,
    "theta0" = theta0.opt,
    "theta1" = resultOptimal$theta1,
    "L" = resultOptimal$L,
    "v0" = v0.opt,
    "A.net.e" = A.net.opt,
    "ebic" = ebic.opt,
    "index" = index.c,
    "lambda.v" = if (!is.null(resultList$lambda.v)) resultList$lambda.v[index.c[1]] else NULL,
    "net.thre" = if (!is.null(resultList$net.thre)) resultList$net.thre[index.c[2]] else NULL
  )
  return(retList)
}
