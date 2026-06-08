\name{ARp_pred}
\alias{ARp_pred}
\title{Forecast using univariate AR(p) fits}
\usage{
ARp_pred(arList, testdata, n.ahead = 1)
}
\arguments{
\item{arList}{A list of fitted AR(p) objects, such as outputs from \code{ar.yw()}.}

\item{testdata}{A numeric matrix of test observations.}

\item{n.ahead}{Number of steps ahead to forecast.}
}
\value{
A list containing predictions and forecasting error summaries.
}
\description{
Forecasts each component of a multivariate time series using component-wise AR(p) fits.
}
