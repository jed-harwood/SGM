#############
##### goodness of fit by parametric bootstrap
##### Added 2024-09-27: Use Step 1 Estimates

########################################
#### Functions used for GF.fit #########
########################################

###### generate data according to GAR(1) model: Sigma^{-1}=Omega=(theta0+theta1*L)^2
GenData.L2 <- function(nobs, theta0, theta1, L, rep = 1) {
  ## nobs: sample size
  ## theta0>0, theta1>0
  ## L: (normalized) Laplacian, p by p matrix
  ## rep: number of replicates
  ## return: list of nobs by p data matrices; p by p concentration matrix: Omega; p by p covariance matrix: Sigma

  library(mnormt)
  data.rep = NULL

  p = nrow(L)
  Omega = (theta0 * diag(1, p) + theta1 * L) %*% (theta0 * diag(1, p) + theta1 * L)
  Sigma = solve(Omega)

  for (i in 1:rep) {
    data = rmnorm(nobs, mean = rep(0, p), varcov = Sigma)
    data.rep[[i]] = data
  }

  list(data = data.rep, Omega = Omega, Sigma = Sigma)
}

validate_gar1_gf_inputs <- function(S, nobs, lambda.v, rho.v, eps_thre, eps_abs, eps_rel, max_iter, num.thread, rep.boot, seed) {
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) != ncol(S) || any(!is.finite(S))) {
    stop("`S` must be a finite numeric square matrix.")
  }
  if (!isTRUE(all.equal(S, t(S), check.attributes = FALSE, tolerance = 1e-8))) {
    stop("`S` must be symmetric.")
  }
  if (!is.numeric(nobs) || length(nobs) != 1 || !is.finite(nobs) || nobs <= 0 || nobs != as.integer(nobs)) {
    stop("`nobs` must be a positive integer scalar.")
  }
  if (!is.numeric(lambda.v) || length(lambda.v) != 1 || !is.finite(lambda.v) || lambda.v <= 0) {
    stop("`lambda.v` must be a positive numeric scalar.")
  }
  if (!is.numeric(rho.v) || length(rho.v) != 1 || !is.finite(rho.v) || rho.v <= 0) {
    stop("`rho.v` must be a positive numeric scalar.")
  }
  if (any(!is.finite(c(eps_thre, eps_abs, eps_rel))) || eps_thre <= 0 || eps_abs <= 0 || eps_rel <= 0) {
    stop("`eps_thre`, `eps_abs`, and `eps_rel` must be positive finite scalars.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || !is.finite(max_iter) || max_iter <= 0 || max_iter != as.integer(max_iter)) {
    stop("`max_iter` must be a positive integer scalar.")
  }
  if (!is.numeric(num.thread) || length(num.thread) != 1 || !is.finite(num.thread) || num.thread <= 0 || num.thread != as.integer(num.thread)) {
    stop("`num.thread` must be a positive integer scalar.")
  }
  if (!is.numeric(rep.boot) || length(rep.boot) != 1 || !is.finite(rep.boot) || rep.boot <= 0 || rep.boot != as.integer(rep.boot)) {
    stop("`rep.boot` must be a positive integer scalar.")
  }
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
    stop("`seed` must be a finite numeric scalar.")
  }
  invisible(TRUE)
}

resolve_loglike_function <- function() {
  if (exists("LogLike", mode = "function", inherits = TRUE)) {
    return(get("LogLike", mode = "function", inherits = TRUE))
  }
  if ("SGM" %in% loadedNamespaces()) {
    return(get("LogLike", envir = asNamespace("SGM"), mode = "function"))
  }
  stop("`LogLike` could not be found. Load the `SGM` package or source the file that defines `LogLike` before calling `GAR1_gf()`.")
}

bootstrap.like <- function(L, theta0, theta1, nobs, lambda.v, rho.v = lambda.v,
                           eps_thre = 1e-6, eps_abs = 1e-5, eps_rel = 1e-3,
                           max_iter = 10000, rep.boot = 100, num.thread = 1,
                           seed = 1, verbose = FALSE) {
  set.seed(seed)
  if (verbose) {
    print("Generating bootstrap samples...")
  }
  temp = GenData.L2(nobs, theta0, theta1, L, rep.boot)
  data.boot = temp$data

  step0_fun = fit_step_0a
  step1_fun = fit_step_1
  log_like_fun = resolve_loglike_function()

  worker_fun = function(i) {
    S.c = var(data.boot[[i]]) * ((nobs - 1) / nobs)
    step0.rep = step0_fun(S.c, nobs)
    step1.rep = step1_fun(
      step0a = step0.rep,
      lambda.v = lambda.v,
      rho.v = rho.v,
      net.thre = 0,
      model = "LN",
      eps_thre = eps_thre,
      eps_abs = eps_abs,
      eps_rel = eps_rel,
      max_iter_s1 = max_iter,
      verbose = FALSE
    )

    fit.rep.i = step1.rep$fit[[1]]
    if (is.null(fit.rep.i) || !isTRUE(fit.rep.i$conv) || length(fit.rep.i$theta1) == 0) {
      return(NA_real_)
    }

    log_like_fun(S.c, step0.rep$theta0, fit.rep.i$theta1, fit.rep.i$L, nobs)
  }

  if (num.thread <= 1) {
    return(vapply(seq_len(rep.boot), worker_fun, numeric(1)))
  }

  foreach::`%dopar%`(
    foreach::foreach(
      i = seq_len(rep.boot),
      .combine = "c",
      .packages = "SGM",
      .export = c("worker_fun", "data.boot")
    ),
    worker_fun(i)
  )
}




#' Goodness of Fit Test
#' 
#' @description
#' This function provides a goodness of fit test to see whether a GAR(1) model with the normalized Laplacian is applicable. It is valid when the signal dimension is at most the number of observations available.
#' 
#' @param S Estimate for covariance matrix, such as the MLE
#' @param nobs The number of observations used to calculate `S`
#' @param lambda.v Tuning parameter for GAR(1). Positive number.
#' @param rho.v ADMM parameter. Positive number.
#' @param eps_thre Small positive number.
#' @param eps_abs ADMM convergence criterion.
#' @param eps_rel ADMM convergence criterion. 
#' @param max_iter Number of iterations to run the Step 1 fit.
#' @param num.thread Number of threads to use for computing.
#' @param rep.boot Number of bootstrap samples to generate for the test.
#' @param seed Random seed used for reproducibility.
#' 
#' @returns p-value for the goodness of fit test 
#' 
#' @export
GAR1_gf = function(S, nobs, lambda.v, rho.v = lambda.v, eps_thre = 1e-6, eps_abs = 1e-5,
                   eps_rel = 1e-3, max_iter = 10000, num.thread = 1, rep.boot = 100, seed = 1) {
  validate_gar1_gf_inputs(S, nobs, lambda.v, rho.v, eps_thre, eps_abs, eps_rel, max_iter, num.thread, rep.boot, seed)

  doParallel::registerDoParallel(num.thread)
  on.exit(doParallel::stopImplicitCluster(), add = TRUE)

  step0 = fit_step_0a(S, nobs)
  step1 = fit_step_1(
    step0a = step0,
    lambda.v = lambda.v,
    rho.v = rho.v,
    net.thre = 0,
    model = "LN",
    eps_thre = eps_thre,
    eps_abs = eps_abs,
    eps_rel = eps_rel,
    max_iter_s1 = max_iter,
    verbose = FALSE
  )

  init.fit = step1$fit[[1]]
  if (is.null(init.fit) || !isTRUE(init.fit$conv) || length(init.fit$theta1) == 0) {
    stop("The observed-data Step 1 fit in `GAR1_gf()` did not converge. Try increasing `max_iter` or adjusting `lambda.v`/`rho.v`.")
  }

  theta0.est = step0$theta0
  theta1.est = init.fit$theta1
  L.est = init.fit$L

  log.like.boot = bootstrap.like(
    L = L.est,
    theta0 = theta0.est,
    theta1 = theta1.est,
    nobs = nobs,
    lambda.v = lambda.v,
    rho.v = rho.v,
    eps_thre = eps_thre,
    eps_abs = eps_abs,
    eps_rel = eps_rel,
    max_iter = max_iter,
    rep.boot = rep.boot,
    num.thread = num.thread,
    seed = seed,
    verbose = FALSE
  )

  if (!any(is.finite(log.like.boot))) {
    stop("All bootstrap refits in `GAR1_gf()` failed to converge. Try increasing `max_iter` or reducing `lambda.v`.")
  }

  log.like.obs = resolve_loglike_function()(S, theta0.est, theta1.est, L.est, nobs)
  mean(log.like.obs > log.like.boot, na.rm = TRUE)
}
