#########################################
## TEMP APPLICATION: ANALYSIS PIPELINE ##
## AUTHOR: JEDIDIAH HARWOOD
## DATE: 7/21/25
#########################################

## G-VAR baseline reference: Isufi, E., Loukas, A., Perraudin, N., and
## Leus, G. (2019). Forecasting Time Series With VARMA Recursions on Graphs.
## IEEE Transactions on Signal Processing, 67(18), 4870-4885.
## https://doi.org/10.1109/TSP.2019.2929930

## Packages 
library(dplyr)
library(readr)
library(ggplot2)
library(zoo)
library(foreach)
library(doParallel)
library(tidyverse)
library(igraph)
library(raster)
library(RSpectra)
library(Matrix)
library(sf)
library(sfnetworks)

if (!requireNamespace("SGM", quietly = TRUE)) {
  stop("The SGM package must be installed before running the temperature application.")
}




temp_detrend = function(path_train, path_test, plot_temp = FALSE, run_prep = T, standardize=F, forecast=FALSE){
  
  
  if (run_prep){
  ## Read in the dataset
  temp = read_csv(path_train,
                  col_types = list("f", "f", "d", "d", "d", "D", "d","d"))
  
  
  ## Explore temperature by station; randomly sample 10 stations
  if (!plot_temp){
  set.seed(123)
  selec.stat = sample(unique(temp$STATION), 10)
  pltDat = temp %>% filter(STATION %in% selec.stat)
  plot_og = ggplot(data = pltDat, aes(x = DATE, y = TEMP, color = STATION)) + 
    geom_line() + 
    ggtitle("Average Daily Temperature for 10 Randomly Selected Stations")
  }
  
  ## Make each station its own column with each rows as a given date
  tempWide = temp %>% dplyr::select(TEMP, STATION, DATE) %>% 
    pivot_wider(names_from = STATION, values_from = TEMP) %>% 
    arrange(DATE) 
  
  ## Remove any stations with more than 50 NA values
  num.missing = colSums(is.na(tempWide[,-1])) %>% 
    sort(decreasing = T)
  ind.rem = which(num.missing < 50)[1]
  rem.stat = num.missing[1:(ind.rem)]
  ## Store names of excessively odd STATIONS
  rem.stat.names = names(rem.stat)
  
  ## Remove STATIONS that have around more than 50 missing values
  tempWide_cln = tempWide %>% dplyr::select(-all_of(rem.stat.names))
  
  ## Imputation for NA values using cubic spline smoother
  tempWide_cln = zoo(tempWide_cln[,-1], order.by = tempWide_cln$DATE)
  tempMat = na.spline(tempWide_cln) %>% as.matrix() 
  tempMat = tempMat[-c(1, nrow(tempMat)), ]
  
  ## Detrending using LOESS smoother
  ## Remove seasonality component 
  ## and trend by loess smoothing applied to each column
  tempMat.tilde = matrix(NA, nrow = nrow(tempMat), 
                         ncol = ncol(tempMat),
                         dimnames = list(rownames(tempMat),
                                         colnames(tempMat)))
  time = rownames(tempMat) %>% 
    as.Date() %>% 
    as.numeric()
  for (j in 1:ncol(tempMat)){
  
    ## Fit loess model
    fit.loess = loess(tempMat[,j] ~ time)
    
    ## Extract Residuals as new time series
    #tmp.res = ts(residuals(fit.loess), frequency = 365)
    #tmp.dseason = stl(tmp.res, s.window = "periodic", l.window = 64)
    if (!standardize){
      tempMat.tilde[, j] = residuals(fit.loess) #tmp.dseason$time.series[,3] # 
    } else{
      resid.spline = residuals(fit.loess)
      tempMat.tilde[, j] = resid.spline/sd(resid.spline)
    }
    
    
  }
  
  if (plot_temp){
    ## Analyze trend
    set.seed(123)
    plot_detrend = matplot(tempMat.tilde[, sample(1:ncol(tempMat), 10)],
            type = "l",
            main = "Loess-detrended Temperature")
  }
  
  ## Now, repeat for the testing year; 
  temp.test = read_csv(path_test)
  
  ## Make each station its own column with each rows as a given date
  tempWide.test = temp.test %>% dplyr::select(TEMP, STATION, DATE) %>% 
    pivot_wider(names_from = STATION, values_from = TEMP) %>% 
    arrange(DATE) 
  
  ## Select columns that appear in year 2020
  tempWide.cln.test = try(tempWide.test %>% dplyr::select(DATE, all_of(colnames(tempMat.tilde))))
  if (inherits(tempWide.cln.test, "try-error")){
    message("Error: The testing set doesn't have the same stations as the training set!")
    return(NULL)
  }
  
  ## Implement column-wise imputation, using cubic spline smoothing
  tempWide.cln.test = zoo(tempWide.cln.test[,-1], order.by = tempWide.cln.test$DATE)
  tempMat.test= na.spline(tempWide.cln.test) %>% as.matrix() 
  
  ## Remove seasonality component 
  ## -and trend by loess smoothing applied to each column
  tempMat.tilde.test = matrix(NA, nrow = nrow(tempMat.test), 
                            ncol = ncol(tempMat.test),
                            dimnames = list(rownames(tempMat.test),
                                            colnames(tempMat.test)))
  time = rownames(tempMat.test) %>% 
    as.Date() %>% 
    as.numeric()
  for (j in 1:ncol(tempMat)){
    
    ## Fit third-order polynomial
    #fit.pm3 = lm(tempMat[,j] ~ poly(time, degree = 3, raw = T))
    
    ## Fit loess model
    fit.loess = loess(tempMat.test[,j] ~ time)
    ## Extract Residuals as new time series
    #tmp.res = ts(residuals(fit.loess), frequency = 365)
    #tmp.dseason = stl(tmp.res, s.window = "periodic", l.window = 64)
    if (!standardize){
      tempMat.tilde.test[, j] = residuals(fit.loess) #tmp.dseason$time.series[,3] # 
    } else{
      resid.spline = residuals(fit.loess)
      tempMat.tilde.test[, j] = resid.spline/sd(resid.spline)
    }
  }
  } else{
    load(path_train)
    if (str_detect(path_train, "2017")){
      tempMat.tilde = temp17
      tempMat.tilde.test = temp17
    } else if (str_detect(path_train, "2015")){ 
      tempMat.tilde = temp15
      tempMat.tilde.test = temp15
    } else if (str_detect(path_train, "2014")){ 
      tempMat.tilde = temp14
      tempMat.tilde.test = temp14
    } else if (str_detect(path_train, "2013")){ 
      tempMat.tilde = temp13
      tempMat.tilde.test = temp13
    } else if (str_detect(path_train, "2012")){ 
      tempMat.tilde = temp12
      tempMat.tilde.test = temp12
    } else if (str_detect(path_train, "2011")){ 
      tempMat.tilde = temp11
      tempMat.tilde.test = temp11
    } else if (str_detect(path_train, "2018")){
      tempMat.tilde = temp18
      tempMat.tilde.test = temp18
    } else if (str_detect(path_train, "2019")){
      tempMat.tilde = temp19
      tempMat.tilde.test = temp19
    } else{
    tempMat.tilde = temp.even
    load(path_test)
    tempMat.tilde.test = temp.odd
    }
    
    ## Standardize
    if (standardize){
      for (j in 1:ncol(tempMat.tilde)){
        tempMat.tilde[,j] = tempMat.tilde[,j]/sd(tempMat.tilde[,j])
        tempMat.tilde.test[,j] = tempMat.tilde.test[,j]/sd(tempMat.tilde.test[,j])
      }
    }
  }
  
  ## If conducting forecasting, split the training year into
  ## two parts so we can conduct forecasting
  if (forecast){
    ## First 2/3-months of the year for training, last for testing.
    tempMat.tilde.test = tempMat.tilde[-c(1:219),]
    tempMat.tilde = tempMat.tilde[c(1:219),]
    
  }
  
  
  resList = list("train" = tempMat.tilde, "test" = tempMat.tilde.test)
  return(resList)
}




fit_targar = function(train, C.lambda = c(8, 4, 1.5, 1, 0.5, 0.25, 0.1), C.thre = exp(seq(log(0.75),log(0.075), length.out=12)), n.cores = 1, stationary=TRUE){
  
  ## Fit the sequence of TAR-GAR(p,q) models to the data, p=1,2,3 and q = 1,2,3
  n.train = nrow(train)
  p = ncol(train)
  num.pass = 3
  
  ## Register cores
  registerDoParallel(n.cores)
  
  ## Initialize TAR-GAR
  TARGAR_list = vector(mode = "list", length = 3)
  
  TARGAR_list = foreach(i=1:3, .maxcombine = 3) %:%
    foreach(j=1:3, .maxcombine = 3, .packages = "SGM") %dopar% {
      lambda.v=C.lambda*sqrt(log(p)/(n.train-i)) 
      rho.v=pmax(lambda.v, 0.01)
      net.thre=C.thre*sqrt(log(p)/(n.train-i)) 
      
      before.time <- proc.time()[3]
      temp <- try(SGM::fit_TAR_GAR(
        data = train,
        p = i,
        q = j,
        lambda.v = lambda.v,
        net.thre = net.thre,
        rho.v = rho.v,
        num.pass = num.pass,
        stationary = stationary
      ))
      end.time = proc.time()[3]
      
      sprintf(paste0("TARGAR(p=", i , " q=", j, ")", " Completed."))
      
      if (inherits(temp, "try-error")){
        print(temp)
      } else{
        temp
      }
    }
  return(TARGAR_list)
}


TARGAR_selec = function(TARGAR_list, p, n.train, C.v, C.thre, heatmap = FALSE){
  
  ## Calculate eBIC Scores
  ebic.tar = array(NA, dim = c(3, 3, length(C.v), length(C.thre)))
  
  ## Matrix to store proportion of models converged for given i,j
  conv.prop.mat = array(NA, dim = c(3,3))
  conv.prop.list = list() 
  
  ## List to store selected models
  tar.models.selec = vector("list", length = 3)
  
  if(p/(n.train-1)>0.5){## e.g., for p=100,n=100
    gamma=1  ## eBIC parameter: p/n~1 set as 1; when p/n<0.5: set as 0.5 
  }else{
    gamma=0.5
  }
  
  P.total=p*(p-1)/2
  
  sgm_loglike = getFromNamespace("LogLike", "SGM")
  sgm_bic = getFromNamespace("BIC", "SGM")
  
  index_selec_list = vector(mode = "list", length = 3)
  
  for (i in 1:3){
    tar.models.selec[[i]] = vector(mode = "list", length = 3)
    conv.prop.list[[i]] = vector(mode = "list", length = 3)
    index_selec_list[[i]] = vector(mode = "list", length = 3)
    for (j in 1:3){
      ## Select model with p=i, q = j
      targar.ij = TARGAR_list[[i]][[j]]$refit
      
      ## Proportion of converged vec
      conv.mat.ij = matrix(NA, nrow = length(C.v), ncol = length(C.thre))
      
      for(l in 1:length(C.v)){ # for lth lambda
        
        for(k in 1:length(C.thre)){ # for kth net.thre
          
          result.c = targar.ij[[l]][[k]]
          result.0S = result.c$result.0S # Step 3 of GAR
          S.c = result.c$S # Final pass estimate of Innovation Covariance matrix
          
          conv.mat.ij[l,k] = if (!is.null(result.0S) && !is.null(result.0S$conv)) {
            result.0S$conv
          } else {
            isTRUE(result.c$conv.3) || isTRUE(result.c$conv.2)
          }
          L.0S = result.c$L
          theta1.0S = result.c$theta1
          theta0.0S = result.c$theta0
          if (is.null(L.0S) && !is.null(result.0S)) {
            L.0S = result.0S$L
            theta1.0S = result.0S$theta1
            theta0.0S = result.0S$theta0
          }
          if (is.null(L.0S) && !is.null(result.c$result.0.post)) {
            L.0S = result.c$result.0.post$L
            theta1.0S = result.c$result.0.post$theta1
            theta0.0S = result.c$theta0.ini
          }
          net.size.c = sum(result.c$A.net)/2 # Network size
          if (is.null(L.0S) || is.null(theta0.0S) || is.null(theta1.0S) ||
              !is.finite(theta0.0S) || !is.finite(theta1.0S)) {
            ebic.tar[i,j,l,k] = NA
            next
          }
          ##############
          ebic.term=2*gamma*(lfactorial(P.total)-lfactorial(net.size.c)-lfactorial(P.total-net.size.c))
          ###############
          
          log.like.0S.c = sgm_loglike(S=S.c, theta0=theta0.0S , theta1 = theta1.0S, L = L.0S, n = n.train-i)
          bic.0S.c = sgm_bic(log.like.0S.c, n.train-i, net.size.c+1+i*(j+1)) ##here sample size is n-1!! (Do not account for the p components of v0); Use net.size + 1 (theta0) + q + 1 (etas)
          ebic.tar[i,j,l,k] = bic.0S.c + ebic.term
          
        }
        
      }
      
      conv.prop.list[[i]][[j]]=  conv.mat.ij
      conv.prop.mat[i,j] = mean(conv.mat.ij)
      
      index.selec = which(ebic.tar[i,j,,] == min(ebic.tar[i,j,,]), arr.ind = T)
      index_selec_list[[i]][[j]] = index.selec
      tar.models.selec[[i]][[j]] = TARGAR_list[[i]][[j]]$refit[[index.selec[1]]][[index.selec[2]]]
      
    }
  }
  
  ## Create heatmaps for convergence of parameters, with selected model starred
  ## Also, extract the adjacency matrix for the selected models
  conv.heatmap.mat = vector(mode = "list", length = 2)
  adj.list = vector(mode = "list", length = 3)
  net.size.mat = matrix(NA, 3, 3)
  
  for (i in 1:3){
    temp.obj = matrix(0, nrow = length(C.v), ncol = length(C.thre))
    conv.heatmap.mat[[i]] = vector(mode = "list", length = 3)
    adj.list[[i]] = vector(mode = "list", length = 3)
    for (j in 1:3){
      temp.obj = conv.prop.list[[i]][[j]]
      conv.heatmap.mat[[i]][[j]] = as.matrix(temp.obj) + 0
      index.c = index_selec_list[[i]][[j]]
      adj.list[[i]][[j]] = TARGAR_list[[i]][[j]]$refit[[index.c[1]]][[index.c[2]]]$A.net
      net.size.mat[i,j] = sum(adj.list[[i]][[j]])/2
    }
  }
  ## Object to store the Plots
  pltList = vector(mode = "list", length=3)
  
  ## Object to store the overlap index
  overlap_targar = matrix(NA, nrow = 9, ncol = 9)
  colnames(overlap_targar) = c(paste0("TARGAR(", 1, ",", 1:3, ")"),
                               paste0("TARGAR(", 2, ",", 1:3, ")"),
                               paste0("TARGAR(", 3, ",", 1:3, ")"))
  rownames(overlap_targar) = c(paste0("TARGAR(", 1, ",", 1:3, ")"),
                               paste0("TARGAR(", 2, ",", 1:3, ")"),
                               paste0("TARGAR(", 3, ",", 1:3, ")"))
  
  
  ## Create the plots and calculate overlap index
  for (i in 1:3){
    pltList[[i]] = vector(mode="list", length = 3)
    for (j in 1:3){
      
      ## Overlap Indices
      A.ij = adj.list[[i]][[j]]
      for (i_other in 1:3){
        for (j_other in 1:3){
          if ((i_other == i) && (j_other == j)){
            next
          }
          A_c.ij = adj.list[[i_other]][[j_other]] # Model to compare to
          
          ## Size of intersection
          intersect.c = A.ij[upper.tri(A.ij, diag = F)] * A_c.ij[upper.tri(A_c.ij, diag = F)]
          
          ## Size of union
          union.c = A.ij[upper.tri(A.ij, diag = F)] | A_c.ij[upper.tri(A_c.ij, diag = F)]
          
          ## overlap index
          overlap_targar[3*(i-1)+j, 3*(i_other-1)+j_other] = sum(intersect.c)/sum(union.c)
        }
      }
      
      ## Heatmaps
      prop.mat = conv.heatmap.mat[[i]][[j]]
      index.ij = index_selec_list[[i]][[j]]
      prop_df = data.frame("Proportion" = c(prop.mat),
                           "C.v" = rep(C.v, times = length(C.thre)),
                           "C.thre" = rep(C.thre, each = length(C.v))
                           )
      prop_df$C.v <- factor(round(prop_df$C.v,2), levels = sort(unique(round(prop_df$C.v,2))))
      prop_df$C.thre <- factor(round(prop_df$C.thre, 2), levels = sort(unique(round(prop_df$C.thre, 2))))
      
      sel_df = data.frame(
        "C.v" = round(C.v[index.ij[1]],2),
        "C.thre" = round(C.thre[index.ij[2]],2)
      )
      sel_df$C.v = factor(sel_df$C.v, levels = levels(prop_df$C.v))
      sel_df$C.thre = factor(sel_df$C.thre, levels = levels(prop_df$C.thre))
      
      # Plot using ggplot2
      plt_ij = ggplot(prop_df, aes(x = C.v, y = C.thre, fill = Proportion)) +
        geom_tile(color = "white") +  # Optional gridlines
        scale_fill_gradient(low = "red", high = "green") +
        labs(title = paste0("Heat Map of Convergence, ", "TARGAR(", i, ",", j, ")"), fill = "Converged") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_point(data = sel_df, aes(x=C.v,  y=C.thre, shape = "Selected Model"), inherit.aes = F) + 
        #geom_point(aes(x = factor(C.v[index.ij[1]]), y = factor(C.thre[index.ij[2]]), shape = "Selected Model"),
        #           size = 4, color = "blue", stroke = 1.2) +  # plotted point
        scale_shape_manual(name = "", values = c("Selected Model" = 8)) + 
        coord_fixed()
      pltList[[i]][[j]] = plt_ij
    }
  }
  
  ## List to return
  resList = list(
    "targar.selec" = tar.models.selec,
    "index.selec" = index_selec_list,
    "conv.prop.mat" = conv.prop.mat,
    "ebic" = ebic.tar,
    "heatmaps" = pltList,
    "overlap" = overlap_targar,
    "netsize" = net.size.mat
  )
  
  return(resList)
}

  


## Function to make predictions
make_forecast_targar = function(modelList, train, test, gvar=FALSE, pq.ebic = NULL, year = 2020){
  ######
  ## Forecasting errors
  
  tar_gar_1s = matrix(NA, nrow = 3, ncol = 3)
  gvar_7nn_1s = tar_gar_1s # New <-- 3/4/2026
  tar_gar_2s = tar_gar_1s
  tar_gar_3s = tar_gar_1s
  
  ## Training errors (1-step)
  tar_gar_1s_train = matrix(NA, nrow = 3, ncol = 3)
  gvar_7nn_1s_train = matrix(NA, nrow = 3, ncol = 3)
  
  for (i in 1:3){
    for (j in 1:3){
      ## Extract coefficient list
      model.ij = modelList[[i]][[j]] # Current model (order i, q=j)
      
      R.list = list(model.ij$R1.0S)
      
      if (i>=2){
        R.list[[2]] = model.ij$R2.0S
      }
      if (i ==3){
        R.list[[3]] = model.ij$R3.0S
      }
      
      ## Forecasting errors
      tar_gar_1s[i,j] = SGM::TARGAR_pred(model = R.list, testdata = test, n.ahead = 1, order = i)$rmse.node[4]
      tar_gar_2s[i,j] = SGM::TARGAR_pred(model = R.list, testdata = test, n.ahead = 2, order = i)$rmse.node[4]
      tar_gar_3s[i,j] = SGM::TARGAR_pred(model = R.list, testdata = test, n.ahead = 3, order = i)$rmse.node[4]
      
      ## Training Error
      tar_gar_1s_train[i,j] =  SGM::TARGAR_pred(model = R.list, testdata = train, n.ahead = 1, order = i)$rmse.node[4]
      
      ### New as of 3/4/2026: Add the G-VAR for corresponding orders; 
      if (exists("temp_knn_matrix_path", mode = "function")) {
        path_7nn = temp_knn_matrix_path(year, 7)
      } else if (year == 2020){
        path_7nn =paste0("data/knn_adj_matrices/adjacency_matrix_k_7.mtx")
      } else{
        path_7nn = paste0("data/knn_adj_matrices/", "tr", year, "/adjacency_matrix_k_7.mtx")
      }
      
      ## Fit the G_VAR(i,j, 7nn)
      A.7nn =  Matrix::readMM(path_7nn)
      gvar.ij = SGM::GVAR_fit(data = train, A = A.7nn, q = i, L_q_vec = rep(j,i))
      
      ## Extract the filter matrices and make predictions
      R.list.ij = gvar.ij$filters
      
      ## Forecasting Error
      gvar_7nn_1s[i,j] = SGM::TARGAR_pred(R.list.ij, test, n.ahead = 1, order = i)$rmse.node[4]
      
      ## Training Error
      gvar_7nn_1s_train[i,j] = SGM::TARGAR_pred(R.list.ij, train, n.ahead = 1, order = i)$rmse.node[4]
    }
  }
  
  ## If gvar and pq.ebic are not null, fit gvar
  if(gvar && !is.null(pq.ebic)){
    ## Read-in the appropriate adjacency matrix
    if (exists("temp_load_knn_adjacencies", mode = "function")) {
      knn = temp_load_knn_adjacencies(year = year, k = c(3, 7, 10), as_matrix = FALSE)
      A.3nn = knn[["3NN"]]
      A.7nn = knn[["7NN"]]
      A.10nn = knn[["10NN"]]
      
    } else if ((year == 2020)){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_10.mtx")
      
    } else if(year == 2019){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_10.mtx")
      
    } else if(year==2018){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_10.mtx")
      
    } else if (year == 2017){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_10.mtx")
      
    } else if (year == 2016){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_10.mtx")
    } else if (year == 2015){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_10.mtx")
    } else if (year == 2014){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_10.mtx")
    } else if (year == 2013){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_10.mtx")
    } else if (year == 2012){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_10.mtx")
    } else if (year == 2011){
      A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_3.mtx")
      A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_7.mtx")
      A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_10.mtx")
    }
    
    ## Now, fit the appropriate g-var model
    order_vec = rep(pq.ebic[2], pq.ebic[1])
    gvar_year_3nn = SGM::GVAR_fit(data = train, A = A.3nn, q = pq.ebic[1], L_q_vec = order_vec)
    R.list.3nn = gvar_year_3nn$filters
    
    gvar_year_7nn = SGM::GVAR_fit(data = train, A = A.7nn, q = pq.ebic[1], L_q_vec = order_vec)
    R.list.7nn = gvar_year_7nn$filters
    
    gvar_year_10nn = SGM::GVAR_fit(data = train, A = A.10nn, q = pq.ebic[1], L_q_vec = order_vec)
    R.list.10nn = gvar_year_10nn$filters
    
    ## Next, fit the 1-step-ahead predictions for the G.VAR models on test set
    gvar_pred_3nn = SGM::TARGAR_pred(R.list.3nn, test, n.ahead = 1, order = pq.ebic[1])$rmse.node[4]
    gvar_pred_7nn = SGM::TARGAR_pred(R.list.7nn, test, n.ahead = 1, order = pq.ebic[1])$rmse.node[4]
    gvar_pred_10nn = SGM::TARGAR_pred(R.list.10nn, test, n.ahead = 1, order = pq.ebic[1])$rmse.node[4]
    
    
  }
  
  resList = list("1-step" = tar_gar_1s,
                 "2-step" = tar_gar_2s,
                 "3-step" = tar_gar_3s)
  if(gvar && !is.null(pq.ebic)){
    resList[["gvar"]] = list("3nn.ebic" = gvar_pred_3nn,
                             "7nn.ebic" = gvar_pred_7nn,
                             "10nn.ebic" = gvar_pred_10nn,
                             "7nn.pq" = gvar_7nn_1s,
                             "7nn.pq.train" = gvar_7nn_1s_train)
  }
  
  return(resList)
}



## Function for Spectral Clustering
## Function to run the spectral clustering
spectral_cluster <- function(A, K, normalized = TRUE, nstart = 10, seed = 1) {
  ## A: Adjacency matrix
  ## K: Number of Clusters
  ## normalized: use normalized Laplacian matrix
  ## nstart: 
  stopifnot(nrow(A) == ncol(A))
  A <- as.matrix(A)
  diag(A) <- 0
  A <- (A + t(A)) / 2  # ensure symmetry
  
  # 1) Build igraph
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # 2) Graph Laplacian (normalized = TRUE → L_sym = I - D^{-1/2} A D^{-1/2})
  L <- laplacian_matrix(g, normalized = normalized, sparse = TRUE)
  
  # 3) K smallest eigenvectors of L
  U <- NULL; evals <- NULL
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    # Fast for large graphs
    es <- RSpectra::eigs_sym(L, k = K, which = "SM")  # "SM" = smallest magnitude
    U <- es$vectors
    evals <- es$values
  } else {
    # Fallback (dense)
    E <- eigen(as.matrix(L), symmetric = TRUE)         # returns eigenvalues in decreasing order
    idx <- order(E$values, decreasing = FALSE)[1:K]    # pick K smallest
    U <- E$vectors[, idx, drop = FALSE]
    evals <- E$values[idx]
  }
  
  # 4) Row-normalize (Ng–Jordan–Weiss)
  if (normalized) {
    rn <- sqrt(rowSums(U^2))
    rn[rn == 0] <- 1
    Tm <- U / rn
  } else {
    Tm <- U
  }
  
  # 5) k-means in the embedding
  set.seed(seed)
  km <- stats::kmeans(Tm, centers = K, nstart = nstart)
  mem <- km$cluster
  
  # Optionally wrap as an igraph "communities" object for plotting
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


spec_cluster_do = function(TARGAR_models, path_train, train, K=7, year=2020, use_base_model=F, output_dir = "results"){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  ## Generate the Lat and Lon Coordinates
  
  ########### <---------- New as of 11/20/2025
  ### Read-in the temp data
  if (str_detect(path_train, ".csv")){
    temp = read_csv(path_train,
                    col_types = list("f", "f", "d", "d", "d", "D", "d","d"))
    
  } else{
    if (year >= 2011 && year <= 2015){
      path_train_fix = if (exists("temp_temperature_path", mode = "function")) {
        temp_temperature_path("temp_CA_2010_2015.csv")
      } else {
        "data/temperature/temp_CA_2010_2015.csv"
      }
      temp = read_csv(path_train_fix,
                      col_types = list("f", "f", "d", "d", "d", "D", "d","d"))
    } else{
    path_train_fix = str_replace(path_train, "Rda", "csv")
    temp = read_csv(path_train_fix,
                    col_types = list("f", "f", "d", "d", "d", "D", "d","d"))
    }
  }
  ########### <------------- End New as of 11/20/2025
  
  ### Load Koppen information
  koppen_path <- if (exists("temp_koppen_raster_path", mode = "function")) {
    temp_koppen_raster_path()
  } else {
    "data/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"
  }
  koppen_raster <- raster(koppen_path)
  
  ### Create dataset for lat and lon
  latlong = temp %>%
    dplyr::filter(STATION %in% colnames(train)) %>%
    dplyr::select(STATION, LATITUDE, LONGITUDE) %>%
    distinct() %>%
    dplyr::rename(name= STATION,
                  x= LONGITUDE,
                  y = LATITUDE)  
  stations_sf <- st_as_sf(latlong, coords = c("x", "y"), crs = 4326)
  
  # Extract climate codes from the raster
  latlong$koppen_code <- raster::extract(koppen_raster, stations_sf)
  latlong$koppen_code = fct_recode(as.factor(latlong$koppen_code), BWh="4", BWk="5", BSh="6", BSk="7", Csa="8", Csb="9", Dsb="18")
  latlong$koppen_code = fct_na_value_to_level(latlong$koppen_code, "uncategorized")
  
  ## NEW AS OF 11/18/2025 <-------------
  # Manually fill-in koppen zones for uncategorized stations
  stations_manual = list("99401899999"="Csa",
                         "72493723289"="Csb",
                         "99402899999"="Csb",
                         "99404199999"="Csb",
                         "72494023234"="Csb",
                         "72594024213"="Csb",
                         "72292800369"="BSk")
  for (i in 1:nrow(latlong)){
    if (latlong$name[i] %in% names(stations_manual)){
      print(i)
      print(latlong$name[i])
      latlong$koppen_code[i] = stations_manual[[as.character(latlong$name[i])]]
    }
  }
  # drop levels for unused factors
  latlong$koppen_code = droplevels(latlong$koppen_code)
  
  ## END NEW AS of 11/28/2025 <---------------
  
  ### Now create the graph objects
  ##### TARGAR
  ebic.pq = apply(TARGAR_models$ebic, MARGIN = c(1,2), function(x){min(x, na.rm = T)})
  index.selec = which(ebic.pq == min(ebic.pq, na.rm = T), arr.ind = T)
  i.selec = index.selec[1]
  j.selec = index.selec[2]
  TARGAR.ebic = TARGAR_models$targar.selec[[i.selec]][[j.selec]]
  
  if (use_base_model){
    TARGAR.ebic = TARGAR_models$targar.selec[[1]][[1]]
  }
  
  A.tar = TARGAR.ebic$A.net
  graph.tar = graph_from_adjacency_matrix(A.tar,
                                          add.colnames = NULL,
                                          diag = FALSE,
                                          mode = "undirected")
  V(graph.tar)$koppen <- latlong$koppen_code
  
  
  
  ###### kNN Graphs
  k_vec = c(3, 7, 10)
  
  if (exists("temp_load_knn_adjacencies", mode = "function")) {
    knn = temp_load_knn_adjacencies(year = year, k = c(3, 7, 10), as_matrix = FALSE)
    A.3nn = knn[["3NN"]]
    A.7nn = knn[["7NN"]]
    A.10nn = knn[["10NN"]]
    
  } else if ((year == 2020)){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/adjacency_matrix_k_10.mtx")
    
  } else if(year == 2019){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2019/adjacency_matrix_k_10.mtx")
    
  } else if(year==2018){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2018/adjacency_matrix_k_10.mtx")
    
  } else if (year == 2017){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2017/adjacency_matrix_k_10.mtx")
    
  } else if (year == 2016){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2016/adjacency_matrix_k_10.mtx")
  } else if (year == 2015){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2015/adjacency_matrix_k_10.mtx")
  } else if (year == 2014){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2014/adjacency_matrix_k_10.mtx")
  } else if (year == 2013){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2013/adjacency_matrix_k_10.mtx")
  } else if (year == 2012){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2012/adjacency_matrix_k_10.mtx")
  } else if (year == 2011){
    A.3nn = Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_3.mtx")
    A.7nn = Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_7.mtx")
    A.10nn = Matrix::readMM("data/knn_adj_matrices/tr2011/adjacency_matrix_k_10.mtx")
  }
  
  A.3nn[A.3nn > 1e-6] = 1
  A.7nn[A.7nn > 1e-6] = 1
  A.10nn[A.10nn > 1e-6] = 1
  
  graph_3nn = graph_from_adjacency_matrix(A.3nn, mode = "undirected", diag = F)
  graph_7nn = graph_from_adjacency_matrix(A.7nn, mode = "undirected", diag = F)
  graph_10nn = graph_from_adjacency_matrix(A.10nn, mode = "undirected", diag = F)
  
  
  V(graph_3nn)$koppen = latlong$koppen_code ## Add the koppen zone
  V(graph_7nn)$koppen = latlong$koppen_code ## Add the koppen zone
  V(graph_10nn)$koppen = latlong$koppen_code ## Add the koppen zone
  
  ## Plot settings
  ##########################################
  ##### Layouts
  lay_3nn = layout_with_fr(graph_3nn)
  lay_7nn = layout_with_fr(graph_7nn)
  lay_10nn = layout_with_fr(graph_10nn)
  lay_tar = layout_with_fr(graph.tar)
  
  #### Colors and frequencies
  col_vec = rainbow(8)
  vertex_col_list = c()
  freq_vec = c() ## <-- New
  j = 1
  for (i in unique(V(graph.tar)$koppen)){
    vertex_col_list[[i]] = col_vec[j]
    freq_vec[j] = table(V(graph.tar)$koppen)[i] ## <-- New
    j = j + 1
  }
  
  col_vec = c()
  for (i in V(graph.tar)$koppen){
    col_vec = c(col_vec, vertex_col_list[[i]])
  }
  V(graph.tar)$col = col_vec
  ##########################################
  
  spec_3nn = spectral_cluster(A = graph_3nn, K = K, normalized = T, nstart = 25)
  spec_7nn = spectral_cluster(A = graph_7nn, K = K, normalized = T, nstart = 25)
  spec_10nn = spectral_cluster(A = graph_10nn, K = K, normalized = T, nstart = 25)
  spec_tar = spectral_cluster(A = graph.tar, K = K, normalized = T, nstart = 25)
  
  ## Calculate Purity and summary statistics from clusters
  majority_3nn = c()
  majority_7nn = c()
  majority_10nn = c()
  majority_tar = c()
  
  majority2_3nn = c()
  majority2_7nn = c()
  majority2_10nn = c()
  majority2_tar = c()
  
  cluster_info = vector(mode = "list", length = K)
  for (i in 1:K){
    ## Obtain indices of nodes in ith cluster
    idx_cluster_3nn = which(spec_3nn$membership == i)
    idx_cluster_7nn = which(spec_7nn$membership == i)
    idx_cluster_10nn = which(spec_10nn$membership == i)
    idx_cluster_tar = which(spec_tar$membership == i)
    
    ## Calculate length of cluster i for each graph
    cluster_3nn_i_len = length(idx_cluster_3nn)
    cluster_7nn_i_len = length(idx_cluster_7nn)
    cluster_10nn_i_len = length(idx_cluster_10nn)
    cluster_tar_i_len = length(idx_cluster_tar)
    
    ## Calculate most prevalent class in ith cluster
    majority_3nn[i] = sort(table(V(graph_3nn)$koppen[idx_cluster_3nn]), decreasing = T)[1]
    majority_7nn[i] = sort(table(V(graph_7nn)$koppen[idx_cluster_7nn]), decreasing = T)[1]
    majority_10nn[i] = sort(table(V(graph_10nn)$koppen[idx_cluster_10nn]), decreasing = T)[1]
    majority_tar[i] = sort(table(V(graph.tar)$koppen[idx_cluster_tar]), decreasing = T)[1]
    
    ## Calculate the second most prevalent class in ith cluster
    majority2_3nn[i] = sort(table(V(graph_3nn)$koppen[idx_cluster_3nn]), decreasing = T)[2]
    majority2_7nn[i] = sort(table(V(graph_7nn)$koppen[idx_cluster_7nn]), decreasing = T)[2]
    majority2_10nn[i] = sort(table(V(graph_10nn)$koppen[idx_cluster_10nn]), decreasing = T)[2]
    majority2_tar[i] = sort(table(V(graph.tar)$koppen[idx_cluster_tar]), decreasing = T)[2]
    
    
    
    ## Store Majority Size and Corresponding Label for Each cluster
    ### FIX: Returns NULL; No actual names; <-- Fixed as of 11/1025
    majority_3nn_zone = names(sort(table(V(graph_3nn)$koppen[idx_cluster_3nn]), decreasing = T)[1])
    majority_7nn_zone = names(sort(table(V(graph_7nn)$koppen[idx_cluster_7nn]), decreasing = T)[1])
    majority_10nn_zone = names(sort(table(V(graph_10nn)$koppen[idx_cluster_10nn]), decreasing = T)[1])
    majority_tar_zone = names(sort(table(V(graph.tar)$koppen[idx_cluster_tar]), decreasing = T)[1])
    
    ### Second Majority Class
    majority_3nn_zone_2 = names(sort(table(V(graph_3nn)$koppen[idx_cluster_3nn]), decreasing = T)[2])
    majority_7nn_zone_2 = names(sort(table(V(graph_7nn)$koppen[idx_cluster_7nn]), decreasing = T)[2])
    majority_10nn_zone_2 = names(sort(table(V(graph_10nn)$koppen[idx_cluster_10nn]), decreasing = T)[2])
    majority_tar_zone_2 = names(sort(table(V(graph.tar)$koppen[idx_cluster_tar]), decreasing = T)[2])
    
    
    ## Store summary statistics in nested list: cluster_info >> graph
    info_3nn = list("size" = cluster_3nn_i_len, "majority size" = majority_3nn[i], "majority zone" = majority_3nn_zone, "majority zone 2" = majority_3nn_zone_2, "second count" = majority2_3nn[i], "membership"=idx_cluster_3nn, "graph"=graph_3nn)
    info_7nn = list("size" = cluster_7nn_i_len, "majority size" = majority_7nn[i], "majority zone" = majority_7nn_zone, "majority zone 2" = majority_7nn_zone_2, "second count" = majority2_7nn[i], "membership"=idx_cluster_7nn, "graph"=graph_7nn)
    info_10nn = list("size" = cluster_10nn_i_len, "majority size" = majority_10nn[i], "majority zone" = majority_10nn_zone, "majority zone 2" = majority_10nn_zone_2, "second count" = majority2_10nn[i], "membership"=idx_cluster_10nn, "graph"=graph_10nn) ## Glitch fixed as of 11/21/2025
    info_tar = list("size" = cluster_tar_i_len, "majority size" = majority_tar[i], "majority zone" = majority_tar_zone, "majority zone 2" = majority_tar_zone_2, "second count"=majority2_tar[i], "membership"=idx_cluster_tar, "graph"=graph.tar)
    
    cluster_info[[i]] = list("3nn" = info_3nn, 
                             "7nn" = info_7nn,
                             "10nn" = info_10nn,
                             "tar" = info_tar)
  }
  
  purity_3nn = sum(majority_3nn)/length(spec_3nn$membership)
  #Purity: 0.4242424
  
  purity_7nn = sum(majority_7nn)/length(spec_7nn$membership)
  #Purity: [1] 0.4318182
  
  purity_10nn = sum(majority_10nn)/length(spec_10nn$membership)
  #Purity:  [1] 0.4469697
  purity_tar = sum(majority_tar)/length(spec_tar$membership)
  #Purity:  0.5227273
  
  save(cluster_info, file = file.path(output_dir, paste0("cluster_info_temp_", year, ".Rda")))
  print("Cluster information saved.")
  
  
  ## 3nn Plot
  pdf(file = file.path(output_dir, paste0("3nn_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
  plot(spec_3nn$communities, spec_3nn$graph,
       col = col_vec, cex = 0.5, vertex.label = "", layout = lay_3nn, 
       vertex.size = 7, main = paste0("Clustering of 3NN Graph (Purity = ",  round(purity_3nn, 3), ")"))
  legend("topright", legend = paste(unique(V(graph.tar)$koppen), paste0("(", freq_vec, ")")), ## <---- new
         pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)")
  dev.off()
  
  ## 7nn Plot
  pdf(file = file.path(output_dir, paste0("7nn_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
  plot(spec_7nn$communities, spec_7nn$graph,
       col = col_vec, cex = 0.5, vertex.label = "", layout = lay_7nn, 
       vertex.size = 7, main = paste0("Clustering of 7NN Graph (Purity = ",  round(purity_7nn, 3), ")"))
  legend("topright", legend = paste(unique(V(graph.tar)$koppen), paste0("(", freq_vec, ")")),  ## <---- New
         pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)")
  dev.off()
  
  ## 10nn Plot
  pdf(file = file.path(output_dir, paste0("10nn_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
  plot(spec_10nn$communities, spec_10nn$graph,
       col = col_vec, cex = 0.5, vertex.label = "", layout = lay_10nn, 
       vertex.size = 7, main = paste0("Clustering of 10NN Graph (Purity = ",  round(purity_10nn, 3), ")"))
  legend("topright", legend = paste(unique(V(graph.tar)$koppen), paste0("(", freq_vec, ")")), ## <--- New
         pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)")
  dev.off()
  
  ## TAR plot
  if (use_base_model){
    model_name = "TAR-GAR(1,1)"
  } else{
    model_name = paste0("TAR-GAR(", i.selec, ",", j.selec, ")")
  }
  pdf(file = file.path(output_dir, paste0("tar_cluster_temp", year, ".pdf")), width = 8.92, height = 4.86)
  plot(spec_tar$communities, spec_tar$graph,
       col = col_vec, cex = 0.5, vertex.label = "", layout = lay_tar, 
       vertex.size = 8, main = paste0("Clustering of ", model_name, " Graph (Purity = ",  round(purity_tar, 3), ")"))
  legend("topright", legend = paste(unique(V(graph.tar)$koppen), paste0("(", freq_vec, ")")), ## <--- new
         pch = 21, pt.bg = unlist(vertex_col_list), title = "True Koppen Zone (Freq)", col = 'black')
  dev.off()
  
  purity_vec = c(purity_tar, purity_3nn, purity_7nn, purity_10nn)
  names(purity_vec) = c("TARGAR", "3NN", "7NN", "10NN")
  return(purity_vec)
}


### Main function

TARGAR_on_temp = function(path_train, path_test, 
                          C.v=c(8, 4, 1.5, 1, 0.5, 0.25, 0.1), 
                          C.thre=exp(seq(log(0.75),log(0.075), length.out=12)), 
                          plot_temp = FALSE, 
                          heatmap = FALSE, 
                          n.cores = 1, 
                          verbose = FALSE,
                          run_prep = TRUE,
                          num.cluster=7,
                          year = 2020,
                          standardize=FALSE,
                          use_base_model = FALSE,
                          forecast=FALSE, 
                          gvar = FALSE,
                          stationary = TRUE,
                          output_dir = "results"
                          ){
  
  ## Preprocess the data
  temp_list = temp_detrend(path_train, path_test, plot_temp, run_prep = run_prep, standardize = standardize, forecast = forecast)
  temp.train = temp_list[["train"]]
  temp.test = temp_list[["test"]]
  if (verbose){
    print("Preprocessing completed...")
  }
  
  ## Fit the TAR-GAR Models 
  TARGAR_Fits = fit_targar(train = temp.train, C.lambda = C.v, C.thre = C.thre, n.cores = n.cores, stationary = stationary)
  if (verbose){
    print("Model Training Completed...")
  }
  
  ## Model Selection and Heatmaps, other TARGAR metrics 
  p = ncol(temp.train)
  n.train = nrow(temp.train)
  TARGAR_models = TARGAR_selec(TARGAR_Fits, p, n.train, C.v, C.thre, heatmap = FALSE)
  
  if (verbose){
    print("Model Selection Completed")
  }
  ## Model Forecasting
  model_selected = TARGAR_models$targar.selec
  ebic.model = apply(TARGAR_models$ebic, c(1,2), min) # model selected by ebic (p,q)
  pq.ebic = which(ebic.model == min(ebic.model, na.rm = T), arr.ind = T) # New
  forecasts = make_forecast_targar(modelList = model_selected, train = temp.train, test = temp.test, gvar = gvar, pq.ebic = pq.ebic, year = year) # New
  
  if (verbose){
    print("Forecasts made")
  }
  
  ## Purity via Koppen Zones
  purity = spec_cluster_do(TARGAR_models = TARGAR_models, path_train = path_train, train = temp.train, K = num.cluster, year = year, use_base_model=use_base_model, output_dir = output_dir)
  
  ## Results List
  resList = TARGAR_models
  resList[["forecasts"]] = forecasts
  resList[["train"]] = temp.train
  resList[["test"]] = temp.test
  resList[["purity"]] = purity
  
  return(resList)
}
  
  
