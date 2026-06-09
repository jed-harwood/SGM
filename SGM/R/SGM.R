#' SGM
#' 
#' A package for learning latent graphs from stationary signals.  This package provides 3 functions
#' 
#' @section SGM functions: 
#' put them here
#' 
#' @references
#' Isufi, E., Loukas, A., Perraudin, N., and Leus, G. (2019).
#' Forecasting Time Series With VARMA Recursions on Graphs.
#' \emph{IEEE Transactions on Signal Processing}, 67(18), 4870-4885.
#' \doi{10.1109/TSP.2019.2929930}
#' 
#' @docType _PACKAGE
#' @name SGM
#' @useDynLib SGM, .registration=TRUE
#' @import Rcpp
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach "%dopar%"
#' @importFrom mnormt rmnorm
NULL
