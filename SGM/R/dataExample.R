#' Example dataset: stocks
#' 
#' @format A 283 by 1007 matrix
#' 
#' @description
#' A matrix with information on 283 stocks across five GICS sectors: 58 from Information Technology, 72 from Consumer Discretionary, 32 from Consumer Staples, 59 from Financials, and 62 from Industrials. The dataset spans January 1, 2007, to January 1, 2011, covering the global financial crisis, with 1007 closing prices per stock.  The numbers reflect the (standardized) log-returns of the S&P 500 stocks, and are based on raw data obtained via the R package `quantmod`.
#' 
#' @docType data
#' 
#' @usage data("stocks") 
#' 
#' @references Peng et al. 2024
#' 
#' @examples
#' data("stocks")
#' n=nrow(stocks)
#' p=ncol(stocks)
#' ### Set the model to fit: GAR(1)-normalized Laplacian model
#' model="LN" 
#' 
#' ### Set tuning parameters: lambda and net.thre sequence
#' C.v=c(8,4,1)
#' lambda.v=C.v*sqrt(log(p)/n)
#' C.thre=exp(seq(log(1),log(0.1), length.out=12))
#' net.thre=C.thre*sqrt(log(p)/n)
#' ### Set ADMM parameter
#' rho.v=pmax(lambda.v, 0.01)
#' 
#' ### Calculate sample covariance
#' S = var(stocks)*(n-1)/n
#' 
#' ### Determine if GAR(1) model is appropriate
#' GAR1_gf(S, n, lambda.v[1])
#' 
#' ### Run GAR1_fit (up to step 3)
#' resList = GAR1_fit(S, n, lambda.v, net.thre, model, 3, rho.v)
#' 
#' ### Conduct model selection via eBIC
#' optModel = model_sele(resList, n, step = 3, model = "LN")
#' 
#' @source Ryan JA, Ulrich JM (2024). quantmod: Quantitative Financial Modelling Framework. R package version 0.4.26, <https://CRAN.R-project.org/package=quantmod>.

"stocks"