% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ADMM_GAR1_fit_proc_package.R
\name{model_selec}
\alias{model_selec}
\title{Select tuning parameters for GAR(1) model}
\usage{
model_selec(resultList, n, step = 3, model = "LN")
}
\arguments{
\item{resultList}{A list output from \code{GAR1_fit}}

\item{n}{An integer referring to the number of observations}

\item{step}{2 or 3; How many steps were used to fit the model.  Requires that at least step 2 was used for \code{GAR1_fit}.}

\item{model}{Which model was fitted
\itemize{
\item \code{"LN"} Normalized graph Laplacian
\item \code{"L"} Graph Laplacian
\item \code{"LN.noselfloop"} Normalized graph Laplacian without self-loops.
\item \verb{"TARGAR} Time-Varying GAR
}}
}
\value{
A list object
\itemize{
\item \code{model.selec}: A list containing the model with the model chosen by the eBIC criterion.
\item \code{A.net.e}: A matrix encoding the (unweighted) graph chosen by the eBIC criterion.
\item \code{index}: Index for the optimal tuning parameters (lambda, net.thre) for the eBIC-selected model.
\item \code{ebic}: ebic score for the selected model.
}
}
\description{
Given a list of models, uses eBIC criterion to select the appropriate tuning parameters for GAR(1) and conducts a goodness of fit test.
}
