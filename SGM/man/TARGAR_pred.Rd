\name{TARGAR_pred}
\alias{TARGAR_pred}
\title{Forecast using TAR-GAR filter matrices}
\usage{
TARGAR_pred(model, testdata, order = 1, n.ahead = 1)
}
\arguments{
\item{model}{A list containing TAR-GAR filter matrices.}

\item{testdata}{A numeric matrix of test observations.}

\item{order}{The autoregressive order of the model.}

\item{n.ahead}{Number of steps ahead to forecast.}
}
\value{
A list containing predictions and forecasting error summaries.
}
\description{
Forecasts a multivariate time series using estimated TAR-GAR filter matrices.
}
