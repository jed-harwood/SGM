####################################################
## GAR Post-processing for Stock Data Application ##
## Outputs Plots and Results from the Paper       ##
## -data Obtained via Loggle R Package            ##
####################################################

## Load required Script
source("GenData.R")
library(SGM)
library(igraph)



## Given GAR model fit, obtain information for plot results
get_info_gar = function(modelList, n){
  
  ## Extract length of tuning parameters 
  n.thre = length(modelList$result.0S[[1]])
  n.lambda = length(modelList$result.0S)
  
  ## Initialize Storage
  net.size = matrix(NA, n.lambda, n.thre)
  log.likelihood.0S = matrix(NA, n.lambda, n.thre)
  
  ## Get index of selected model
  GAR.index.selec = model_selec(modelList, n)$index
  
  
  
  ## Calculate necessary quantities
  S = modelList$S # sample covariance matrix
  for (i in 1:n.lambda){ # ith lambda
    for (j in 1:n.thre){ # jth net.thre
      
      ## Extract estimated parameters
      L.ij = modelList$result.0S[[i]][[j]]$L
      theta0.ij = modelList$result.0S[[i]][[j]]$theta0
      theta1.ij = modelList$result.0S[[i]][[j]]$theta1
      
      ## Store the network size
      net.size[i,j] = sum(modelList$A.0.net[[i]][[j]])/2
      log.likelihood.0S[i,j] = LogLike(S = S, theta0 = theta0.ij, theta1 = theta1.ij, L = L.ij, n = n)
      
    }
  }
  
  GAR.size.selec = net.size[GAR.index.selec[1], GAR.index.selec[2]]
  
  resList = list("size" = net.size, "loglike" = log.likelihood.0S, 
                 "selec" = GAR.index.selec, "size.selec" = GAR.size.selec)
  return(resList)
}

## Given GLASSO model fit, obtain information for plot results

get_info_glasso = function(modelList, S, p, n){
  
  ## Number of lambda results
  n.lambda = length(modelList)
  
  ## Initialize Storage Results
  net.glasso.size = rep(NA, n.lambda) # network sizes
  loglike.glasso = rep(NA, n.lambda) # log-likelihood
  ebic.glasso = rep(NA, n.lambda) # eBIC scores
  
  
  ## quantities for eBIC selection
  P.total=p*(p-1)/2
  gam=0.5 ## eBIC parameter: now p<<n, choose small gamma 
  
  
  ## Calculate log-likelihood, eBIC, and network sizes
  for (j in 1:n.lambda){
    A.glasso.net.j=abs(modelList[[j]]$wi)>1e-6
    net.glasso.size[j]=(sum(A.glasso.net.j)-p)/2
    
    ## log-likelihood
    loglike.glasso[j]=LogLike.glasso(S=S, Omega=modelList[[j]]$wi, n=n)
    
    ## bic
    bic.j = -2*loglike.glasso[j]+(net.glasso.size[j]+p)*log(n)
    
    ## ebic
    ebic.term=2*gam*(lfactorial(P.total)-lfactorial(net.glasso.size[j])-lfactorial(P.total-net.glasso.size[j]))
    ebic.glasso[j] = bic.j + ebic.term
    
  }
  
  ## Determine which model minimized the eBIC
  index.lambda.min = which.min(ebic.glasso)
  glasso.ebic.size = net.glasso.size[index.lambda.min]
 
  ## return results
  resList = list("size" = net.glasso.size, "loglike" = loglike.glasso,
                 "index.selec" = index.lambda.min, "size.selec" = glasso.ebic.size)
  return(resList)
}


## Make plot for GAR estimated model
GAR_plot_stocks = function(modelList, n, sp.num){
  
  ## Conduct model selection
  GAR.selec = model_selec(modelList, n = n)
  
  ## Pull the eBIC selected parameters
  L.est = (-1) * GAR.selec$model.selec$W # (thresholded version)
  theta0.est = GAR.selec$model.selec$theta0
  theta1.est = GAR.selec$model.selec$theta1
  v0.est = GAR.selec$model.selec$v0
  adj.est = GAR.selec$A.net.e
  
  
  ## Create stock-wise plots and weights
  A.est=-diag(v0.est)%*%L.est%*%diag(v0.est) # off diagonal proportional to the underlying adjacency matrix
  diag(A.est)=0
  A.est=A.est/max(A.est)  ## rescale to [0,1] so the elements can be interpreted as weights

  ## stock-wise plot: only draw "large" edges 
  adj.est.plot= abs(L.est)>=0.1 ##only plot edges with L-entries > a threshold
  diag(adj.est.plot)=0

  ##sector-wise proportions of intra-sector edges and proportion of inter-sector edges (see loggle paper for details)
  sector.num <- length(sp.num)
  sp.ind <- c(0, cumsum(sp.num))
  edge.sector <- matrix(NA, sector.num, sector.num)
  edge.weight.sector <- matrix(NA, sector.num, sector.num)
  colnames(edge.sector)= names(sp.num)
  rownames(edge.sector)= names(sp.num)
  
  for(i in 1:sector.num) {##edge number/weights with/between sectors 
    for (j in 1:sector.num){
      if(i==j){
        edge.sector[i,i]<-(sum(adj.est[(sp.ind[i]+1):sp.ind[i+1],(sp.ind[i]+1):sp.ind[i+1]]))/2  ## only depend on the 0-1 adjacency matrix 
        edge.weight.sector[i,i]<-(sum(A.est[(sp.ind[i]+1):sp.ind[i+1],(sp.ind[i]+1):sp.ind[i+1]]))/2  ## depends on the weighted adjacency matrix
      }else{
        edge.sector[i,j] <- sum(adj.est[(sp.ind[i]+1):sp.ind[i+1],(sp.ind[j]+1):sp.ind[j+1]])
        edge.weight.sector[i,j] <- sum(A.est[(sp.ind[i]+1):sp.ind[i+1],(sp.ind[j]+1):sp.ind[j+1]])
      }
    }
  }
  
  inter.num=outer(sp.num, sp.num, FUN="*")
  sec.conn=edge.sector/inter.num  ## sector connectivity: only depend on the 0-1 adjacency matrix  
  diag(sec.conn)=diag(edge.sector)/(sp.num*(sp.num-1)/2) ##intra-section connectivity
  
  ## Result for table 4
  table.connect = round(sec.conn*100,1) ## more intra-sector connectivity than inter-section connectivity 

  
  ## Plot sector-wise network
  ##weighted edge connectivity :  depends on the weighted adjacency matrix
  sec.weight.conn=edge.weight.sector/inter.num  ## sector connectivity 
  diag(sec.weight.conn)=diag(edge.weight.sector)/(sp.num*(sp.num-1)/2) ##intra-section connectivity
  #round(sec.weight.conn*1000 ,3) ## only relative size meaningful as A.est is only up to a multiplier 
  
  ## section-wise plot: edge width proportional to between sector  connectivity  
  #sec.conn.u=sec.conn ## use unweighted sector connectivity
  sec.conn.u=sec.weight.conn ## use weighted sector connectivity
  
  adj.sector=(sec.conn.u>0)
  diag(adj.sector)=0
  net.sector <- graph.adjacency(adj.sector, mode = "undirected", diag = FALSE)
  
  V(net.sector)$color <- rainbow(sector.num)
  V(net.sector)$size<-diag(sec.conn.u)/max(diag(sec.conn.u))*50 ## vertices size proportional to with-sector connectivity
  E(net.sector)$color <- "darkgray"
  temp=sec.conn.u[lower.tri(sec.conn.u)]/max(sec.conn.u[lower.tri(sec.conn.u)])*5
  temp=temp[adj.sector[lower.tri(adj.sector)]==1]
  E(net.sector)$width <-temp ## edge size proportional to between size connectivity
  
  pdf("stock-GAR1LN-sector-wise-weighted-network.pdf")
  set.seed(301)
  plot(net.sector, vertex.label = names(sp.num), layout = layout.fruchterman.reingold, 
       main = "GAR1 estimated network: sector wise ")
  legend("bottomright",  names(sp.num), pch = 19, 
         col = rainbow(sector.num), cex = 0.5) 
  dev.off()
  
  
  
  return(table.connect)
}




### Main wrapper function
stock_results = function(modelList, model = "GAR", p, n, sp.num, S){
  
  if (model == "GAR"){
    resList = get_info_gar(modelList = modelList, n = n)
    table = GAR_plot_stocks(modelList = modelList, n = n, sp.num = sp.num)
    resList$table = table
  } else if (model == "glasso"){
    resList = get_info_glasso(modelList = modelList, S = S, p = p, n=n)
  } else{
    message("Model not supported.")
  }
  
  
  return(resList)
}










