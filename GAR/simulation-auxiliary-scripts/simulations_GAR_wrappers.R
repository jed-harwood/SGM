setup_simulation_parallel <- function(num.thread) {
  registerDoFuture()
  plan(multisession, workers = num.thread)
  options(future.globals.maxSize = 2 * 1024^3)
}


simulation_data_setup <- function(p, n, rep, model, theta0, theta1) {
  filename = paste("sep_GAR1_", model, "_p", p, "_n", n, "_seqlamnet.Rdata", sep = "")

  set.seed(1)
  edge.prob = 2 / p
  A.tr = Rand.Graph(p = p, edge.prob = edge.prob, self.prob = 2 * edge.prob, min = 0.5, max = 1, selfloop = FALSE, isolate = FALSE)
  net.tr = (A.tr > 0)
  diag(net.tr) = 0
  deg = apply(A.tr, 1, sum)
  L = Laplacian(A.tr)
  LN = Laplacian.Norm(A.tr)

  if (model == "LN" || model == "LN.noloop") {
    L.tr = LN
  } else {
    L.tr = L
  }

  theta0.tr = theta0
  theta1.tr = theta1

  if (model == "LN" || model == "LN.noloop") {
    v0 = deg^{0.5}
    v0 = matrix(v0 / sqrt(sum(v0^2)), p, 1)
  } else {
    v0 = rep(1, p)
  }
  v0.tr = v0

  set.seed(5)
  if (model == "LN" || model == "LN.noloop") {
    temp = GenData.L2(n, theta0, theta1, LN, rep)
  } else {
    temp = GenData.L2(n, theta0, theta1, L, rep)
  }

  data = temp$data
  Omega.tr = temp$Omega
  Sigma.tr = temp$Sigma
  A.ggm.tr = abs(Omega.tr) > 1e-6
  diag(A.ggm.tr) = 0

  list(
    p = p,
    n = n,
    rep = rep,
    model = model,
    theta0 = theta0,
    theta1 = theta1,
    filename = filename,
    edge.prob = edge.prob,
    A.tr = A.tr,
    net.tr = net.tr,
    deg = deg,
    L = L,
    LN = LN,
    L.tr = L.tr,
    theta0.tr = theta0.tr,
    theta1.tr = theta1.tr,
    v0.tr = v0.tr,
    data = data,
    Omega.tr = Omega.tr,
    Sigma.tr = Sigma.tr,
    A.ggm.tr = A.ggm.tr
  )
}


simulation_gar_setup <- function(p, n) {
  C.v = c(1, 0.5)
  lambda.v = C.v * sqrt(log(p) / n)
  rho.v = pmax(lambda.v, 0.01)

  if (p == 100) {
    C.thre = exp(seq(log(1), log(0.05), length.out = 10))
  } else if (p == 250) {
    C.thre = exp(seq(log(1), log(0.075), length.out = 10))
  } else {
    C.thre = exp(seq(log(1), log(0.1), length.out = 10))
  }
  net.thre = C.thre * sqrt(log(p) / n)

  list(C.v = C.v, lambda.v = lambda.v, rho.v = rho.v, C.thre = C.thre, net.thre = net.thre)
}


simulation_fit_gar <- function(data, n, rep, model, lambda.v, net.thre, rho.v) {
  results.GAR = foreach(i = 1:rep, .maxcombine = max(rep, 2)) %dopar% {
    print(paste("Fitting GAR for replicate: ", i))

    S.i = var(data[[i]]) * ((n - 1) / n)
    GAR.i.res = SGM::GAR1_fit(S = S.i, nobs = n, lambda.v = lambda.v, net.thre = net.thre,
                              model = model, rho.v = rho.v, eps_thre = 1e-6, eps_abs = 1e-5,
                              eps_rel = 1e-3, max_iter_s3 = 10000, verbose = FALSE)

    res.i = list("S" = S.i, "conv" = GAR.i.res$conv, "step3b" = GAR.i.res$step3b,
                 "modelList" = GAR.i.res, "A.net" = GAR.i.res$A.net,
                 "step1" = GAR.i.res$step1)
    res.i
  }

  results.GAR
}


simulation_process_gar_results <- function(results.GAR, rep, net.tr, L.tr, theta0.tr, theta1.tr, v0.tr, A.ggm.tr, Sigma.tr, Omega.tr) {
  GAR.ebic = vector(mode = "list", length = rep)
  L.ebic.err = rep(NA, rep)
  theta0.ebic.err = rep(NA, rep)
  v0.ebic.err = rep(NA, rep)
  power.ebic.vec = rep(NA, rep)
  fdr.ebic.vec = rep(NA, rep)
  F1.ebic.vec = rep(NA, rep)
  L.s1.err = rep(NA, rep)
  theta0.s1.err = rep(NA, rep)

  Sigma.gar.err = rep(NA, rep)
  Omega.gar.err = rep(NA, rep)
  power.ggm.gar = rep(NA, rep)
  FDR.ggm.gar = rep(NA, rep)
  F1.ggm.gar = rep(NA, rep)

  for (i in 1:rep) {
    GAR.models.i = results.GAR[[i]]$modelList
    GAR.ebic.i = SGM::model_selec(resultList = GAR.models.i)
    GAR.ebic[[i]] = GAR.ebic.i$selected.model

    A.ebic.i = GAR.ebic.i$A.net.e
    L.ebic.i = GAR.ebic.i$L * GAR.ebic.i$theta1
    theta0.ebic.i = GAR.ebic.i$theta0
    v0.ebic.i = GAR.ebic.i$v0
    v0.ebic.i = v0.ebic.i / sqrt(sum(v0.ebic.i^2))

    net.size.ebic = sum(A.ebic.i > 0) / 2
    L.ebic.err[i] = sum((L.ebic.i - L.tr * theta1.tr)^2) / sum((theta1.tr * L.tr)^2)
    theta0.ebic.err[i] = abs(theta0.ebic.i - theta0.tr)^2
    v0.ebic.err[i] = sum(abs(v0.ebic.i - v0.tr)^2)
    power.ebic.vec[i] = sum(A.ebic.i * net.tr) / sum(net.tr)
    fdr.ebic.vec[i] = sum(A.ebic.i * (1 - net.tr)) / sum(A.ebic.i)
    F1.ebic.vec[i] = (2 * (1 - fdr.ebic.vec[i]) * power.ebic.vec[i]) / (1 - fdr.ebic.vec[i] + power.ebic.vec[i])

    L.gar.s1 = results.GAR[[i]]$step1[[1]]$L * results.GAR[[i]]$step1[[1]]$theta1
    theta0.gar.s1 = results.GAR[[i]]$step1[[1]]$theta0

    theta0.s1.err[i] = abs(theta0.gar.s1 - theta0.tr)^2
    L.s1.err[i] = sum((L.gar.s1 - L.tr * theta1.tr)^2) / sum((theta1.tr * L.tr)^2)

    temp.ggm = GGM.GAR(GAR.ebic.i$L, theta0.ebic.i, GAR.ebic.i$theta1)
    Sigma.gar.err[i] = sum((temp.ggm$Sigma - Sigma.tr)^2) / sum(Sigma.tr^2)
    Omega.gar.err[i] = sum((temp.ggm$Omega - Omega.tr)^2) / sum(Omega.tr^2)

    temp.L = (-1) * GAR.ebic.i$selected.model$W
    temp.ggm = GGM.GAR(temp.L, theta0.ebic.i, GAR.ebic.i$theta1)

    A.ggm = abs(temp.ggm$Omega) > 1e-6
    diag(A.ggm) = 0

    FDR.ggm.gar[i] = sum(A.ggm * (1 - A.ggm.tr)) / sum(A.ggm)
    power.ggm.gar[i] = sum(A.ggm * A.ggm.tr) / sum(A.ggm.tr)
    F1.ggm.gar[i] = (2 * (1 - FDR.ggm.gar[i]) * power.ggm.gar[i]) / (1 - FDR.ggm.gar[i] + power.ggm.gar[i])
  }

  list(
    GAR.ebic = GAR.ebic,
    L.ebic.err = L.ebic.err,
    theta0.ebic.err = theta0.ebic.err,
    v0.ebic.err = v0.ebic.err,
    power.ebic.vec = power.ebic.vec,
    fdr.ebic.vec = fdr.ebic.vec,
    F1.ebic.vec = F1.ebic.vec,
    L.s1.err = L.s1.err,
    theta0.s1.err = theta0.s1.err,
    Sigma.gar.err = Sigma.gar.err,
    Omega.gar.err = Omega.gar.err,
    power.ggm.gar = power.ggm.gar,
    FDR.ggm.gar = FDR.ggm.gar,
    F1.ggm.gar = F1.ggm.gar
  )
}


simulation_fit_oracle <- function(data, rep, p, n, theta0.tr, v0.tr, net.tr, model, theta1.tr, L.tr) {
  Z <- matrix(0, p, p)
  W <- Z

  result.orc.tr = NULL

  print(paste("Fitting GAR-Oracle with true v0, theta0 and Graph-Zero-pattern "))
  result.orc.tr = foreach(i = 1:rep, .maxcombine = max(rep, 2)) %dopar% {
    print(paste("fit replicate ", i))
    S = var(data[[i]]) * ((n - 1) / n)

    SGM::ADMM_L2_Zero(S, theta0.tr, v0.tr, rho = sqrt(log(p) / n), A = net.tr, model = model, Z_ini = Z, W_ini = W,
                      eps_thre = 1e-6, eps_abs = 1e-5, eps_rel = 1e-3, max_iter = 10000, verbose = FALSE)
  }

  L.oracle.err = rep(NA, rep)
  for (i in 1:rep) {
    L.orc.i = result.orc.tr[[i]]$L * result.orc.tr[[i]]$theta1
    L.oracle.err[i] = sum((L.orc.i - theta1.tr * L.tr)^2) / sum((theta1.tr * L.tr)^2)
  }

  list(result.orc.tr = result.orc.tr, L.oracle.err = L.oracle.err)
}


simulation_fit_glasso <- function(data, rep, p, n, A.ggm.tr) {
  C.glasso = exp(seq(log(1.5), log(.15), length.out = 100))
  lambda.glasso = C.glasso * sqrt(log(p) / n)

  result.glasso = vector("list", rep)
  result.glasso = foreach(i = 1:rep, .maxcombine = max(length(rep), 2)) %dopar% {
    result.glasso.i = vector("list", length(lambda.glasso))

    S = var(data[[i]]) * ((n - 1) / n)
    for (j in 1:length(lambda.glasso)) {
      print(paste("glasso replicate ", i))
      result.glasso.i[[j]] = glasso(s = S, nobs = n, rho = lambda.glasso[j])
    }
    result.glasso.i
  }

  A.glasso.net = vector("list", rep)
  edge.zero = vector("list", rep)
  net.glasso.size = vector("list", rep)
  fdr.glasso = matrix(NA, nrow = rep, ncol = length(lambda.glasso))
  power.glasso = matrix(NA, nrow = rep, ncol = length(lambda.glasso))
  F1.glasso = matrix(NA, nrow = rep, ncol = length(lambda.glasso))

  for (i in 1:rep) {
    A.glasso.net.i = vector("list", length(lambda.glasso))
    net.glasso.size.i = rep(NA, length(lambda.glasso))
    edge.zero.i = vector("list", length(lambda.glasso))

    for (j in 1:length(lambda.glasso)) {
      temp = (result.glasso[[i]][[j]]$wi + t(result.glasso[[i]][[j]]$wi)) / 2
      A.glasso.net.i[[j]] = abs(temp) > 1e-6
      net.glasso.size.i[j] = (sum(A.glasso.net.i[[j]]) - p) / 2
      edge.zero.i[[j]] = which(A.glasso.net.i[[j]] == FALSE, arr.ind = TRUE)

      A.glasso.c = A.glasso.net.i[[j]]
      diag(A.glasso.c) = 0

      fdr.glasso[i, j] = sum((A.glasso.c * (1 - A.ggm.tr))) / sum(A.glasso.c)
      power.glasso[i, j] = sum(A.glasso.c * A.ggm.tr) / sum(A.ggm.tr)
      F1.glasso[i, j] = (2 * (1 - fdr.glasso[i, j]) * power.glasso[i, j]) / (1 - fdr.glasso[i, j] + power.glasso[i, j])
    }

    A.glasso.net[[i]] = A.glasso.net.i
    net.glasso.size[[i]] = net.glasso.size.i
    edge.zero[[i]] = edge.zero.i
  }

  result.glasso.refit = vector("list", rep)
  result.glasso.refit = foreach(i = 1:rep, .maxcombine = max(length(rep), 2)) %dopar% {
    result.glasso.refit.i = vector("list", length(lambda.glasso))
    S = var(data[[i]]) * (n - 1) / n

    for (j in 1:length(lambda.glasso)) {
      result.glasso.refit.i[[j]] = glasso(s = S, rho = 1e-10, nobs = n, zero = edge.zero[[i]][[j]], penalize.diagonal = FALSE)
    }

    result.glasso.refit.i
  }

  list(
    C.glasso = C.glasso,
    lambda.glasso = lambda.glasso,
    result.glasso = result.glasso,
    A.glasso.net = A.glasso.net,
    edge.zero = edge.zero,
    net.glasso.size = net.glasso.size,
    fdr.glasso = fdr.glasso,
    power.glasso = power.glasso,
    F1.glasso = F1.glasso,
    result.glasso.refit = result.glasso.refit
  )
}


simulation_process_glasso_results <- function(glasso_fit, data, rep, p, n, Sigma.tr, Omega.tr) {
  lambda.glasso = glasso_fit$lambda.glasso
  net.glasso.size = glasso_fit$net.glasso.size
  result.glasso.refit = glasso_fit$result.glasso.refit
  fdr.glasso = glasso_fit$fdr.glasso
  power.glasso = glasso_fit$power.glasso
  F1.glasso = glasso_fit$F1.glasso

  if (p / n > 0.5) {
    gamma = 1
  } else {
    gamma = 0.5
  }

  Sigma.glasso.err = matrix(NA, rep, length(lambda.glasso))
  Omega.glasso.err = Sigma.glasso.err

  log.like.glasso = matrix(NA, nrow = rep, ncol = length(lambda.glasso))
  bic.glasso = matrix(Inf, nrow = rep, ncol = length(lambda.glasso))
  ebic.glasso = matrix(Inf, nrow = rep, ncol = length(lambda.glasso))
  P.total = p * (p - 1) / 2

  Sigma.sele.glasso = rep(NA, rep)
  Omega.sele.glasso = rep(NA, rep)
  fdr.sele.glasso = rep(NA, rep)
  power.sele.glasso = rep(NA, rep)
  F1.sele.glasso = rep(NA, rep)
  ebic.sele.glasso = rep(NA, rep)
  size.sele.glasso = rep(NA, rep)

  for (i in 1:rep) {
    S = var(data[[i]]) * ((n - 1) / n)
    for (j in 1:length(lambda.glasso)) {
      temp = (result.glasso.refit[[i]][[j]]$wi + t(result.glasso.refit[[i]][[j]]$wi)) / 2

      Omega.glasso.err[i, j] = sum((temp - Omega.tr)^2) / sum(Omega.tr^2)
      Sigma.e = solve(temp)
      Sigma.glasso.err[i, j] = sum((Sigma.e - Sigma.tr)^2) / sum(Sigma.tr^2)

      log.like.glasso[i, j] = LogLike.glasso(S = S, Omega = temp, n = n)
      bic.glasso[i, j] = -2 * log.like.glasso[i, j] + (net.glasso.size[[i]][j] + p) * log(n)
      ebic.glasso[i, j] = bic.glasso[i, j] + gamma * (lfactorial(P.total) - lfactorial(net.glasso.size[[i]][j]) - lfactorial(P.total - net.glasso.size[[i]][j]))
    }
  }

  for (i in 1:rep) {
    if (!all(is.infinite(bic.glasso[i, ]))) {
      j.sele = which.min(ebic.glasso[i, ])

      Sigma.sele.glasso[i] = Sigma.glasso.err[i, j.sele]
      Omega.sele.glasso[i] = Omega.glasso.err[i, j.sele]
      fdr.sele.glasso[i] = fdr.glasso[i, j.sele]
      power.sele.glasso[i] = power.glasso[i, j.sele]
      F1.sele.glasso[i] = F1.glasso[i, j.sele]
      size.sele.glasso[i] = net.glasso.size[[i]][j.sele]
      ebic.sele.glasso[i] = ebic.glasso[i, j.sele]
    }
  }

  list(
    gamma = gamma,
    Sigma.glasso.err = Sigma.glasso.err,
    Omega.glasso.err = Omega.glasso.err,
    log.like.glasso = log.like.glasso,
    bic.glasso = bic.glasso,
    ebic.glasso = ebic.glasso,
    P.total = P.total,
    Sigma.sele.glasso = Sigma.sele.glasso,
    Omega.sele.glasso = Omega.sele.glasso,
    fdr.sele.glasso = fdr.sele.glasso,
    power.sele.glasso = power.sele.glasso,
    F1.sele.glasso = F1.sele.glasso,
    ebic.sele.glasso = ebic.sele.glasso,
    size.sele.glasso = size.sele.glasso
  )
}


simulation_results_tables <- function(gar_results, oracle_results, glasso_results) {
  table1.metric.names = c("theta0.err", "v0.err", "L.err", "Power", "FDR", "F1")
  table1.metric = c(mean(gar_results$theta0.ebic.err), mean(gar_results$v0.ebic.err), mean(gar_results$L.ebic.err), mean(gar_results$power.ebic.vec),
                    mean(gar_results$fdr.ebic.vec), mean(gar_results$F1.ebic.vec))
  names(table1.metric) = table1.metric.names

  table2.metric.names = c("theta0.ini", "L.Step1.err", "L.Oracle.err")
  table2.metric = matrix(NA, nrow = 1, ncol = 3)
  colnames(table2.metric) = table2.metric.names
  table2.metric[1, ] = c(mean(gar_results$theta0.s1.err), mean(gar_results$L.s1.err), mean(oracle_results$L.oracle.err))

  table3.metric = matrix(NA, nrow = 2, ncol = 5)
  colnames(table3.metric) = c("Sigma.err", "Omega.err", "FDR", "Power", "F1")
  rownames(table3.metric) = c("GAR", "GLASSO")
  table3.metric[1, ] = c(mean(gar_results$Sigma.gar.err), mean(gar_results$Omega.gar.err), mean(gar_results$FDR.ggm.gar),
                        mean(gar_results$power.ggm.gar), mean(gar_results$F1.ggm.gar))
  table3.metric[2, ] = c(mean(glasso_results$Sigma.sele.glasso), mean(glasso_results$Omega.sele.glasso), mean(glasso_results$fdr.sele.glasso),
                        mean(glasso_results$power.sele.glasso), mean(glasso_results$F1.sele.glasso))

  list(table1.metric = table1.metric, table2.metric = table2.metric, table3.metric = table3.metric)
}
