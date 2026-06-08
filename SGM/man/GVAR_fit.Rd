\name{GVAR_fit}
\alias{GVAR_fit}
\alias{G.VAR.Fit}
\title{Fit an Isufi-style graph VAR model}
\usage{
GVAR_fit(data, A, q, L_q_vec)

G.VAR.Fit(data, A, q, L_q_vec)
}
\arguments{
\item{data}{An n by d numeric data matrix, where d is the signal dimension.}

\item{A}{A d by d adjacency matrix.}

\item{q}{The autoregressive order.}

\item{L_q_vec}{Polynomial orders in the graph Laplacian for each lag.}
}
\value{
A list with graph filter coefficients, filter matrices, and AIC.
}
\description{
Fits a graph VAR model with filters parameterized by powers of a supplied graph Laplacian.
}
