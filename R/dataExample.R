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
#' data(stocks)
#' @source Ryan JA, Ulrich JM (2024). quantmod: Quantitative Financial Modelling Framework. R package version 0.4.26, <https://CRAN.R-project.org/package=quantmod>.

"stocks"