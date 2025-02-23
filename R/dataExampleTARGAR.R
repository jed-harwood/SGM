#' Example dataset: targar
#'
#' A dataset generated from a TARGAR model.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{data}{A 100 by 100 matrix with data simulated from an underlying TARGAR model, with parameters: \eqn{\theta_0 = 1}, \eqn{\theta_1 = 2}, \eqn{\eta_0 = 0.5}, \eqn{\eta_1 = 0.4 / \lambda_{\max}(L)}.}
#'   \item{A.tr}{The true weighted adjacency matrix (100 by 100).}
#'   \item{LN}{The true normalized Laplacian matrix (100 by 100).}
#'   \item{theta0}{The true value of \eqn{\theta_0}, a positive number.}
#'   \item{theta1}{The true value of \eqn{\theta_1}, a positive number.}
#'   \item{eta0}{The true value of \eqn{\eta_0}.}
#'   \item{eta1}{The true value of \eqn{\eta_1}.}
#' }
#'
#' @usage data("targar")
#' @docType data
#' @references Peng et al. (2024)
#' @source Simulated data from a TARGAR model
"targar"
