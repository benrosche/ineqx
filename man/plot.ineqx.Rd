% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.ineqx.R
\name{plot.ineqx}
\alias{plot.ineqx}
\title{plot function}
\usage{
\method{plot}{ineqx}(ineqx.out, type, yscale = 1)
}
\arguments{
\item{ineqx.out}{ineqx.out object from ineqx()}

\item{type}{Character string. Plot type. Choose from: dMuP, dMuT, dSigmaP, dSigmaT, dWP, dWT, dBP, dBT, dPP, dPT, dT, dPA}

\item{yscale}{Either 1 or 2. Choose value on y-axis. 1: Effect in units of ystat. 2: Effect as \% of value at reference time (ref)}
}
\value{
Returns a ggplot2 object
}
\description{
[...]
}
\details{
The y-axis can be scaled in two ways:

\code{yscale=1} (default): absolute points of Var or CV2.\cr
\code{yscale=2}: \% of the value of Var or CV2 at reference time. \cr
(Total Var or CV2 with \code{type="dT"} and within-group & between-group Var or CV2 with \code{type="dW"} & \code{type="dB"}, respectively.)
}
\examples{
data(incdat)
i1 <- ineqx(y=inc, group=SES, time=year, ref=1, dat=incdat)
plot(i1, type="dMuSigma2")

}
\author{
Benjamin Rosche <benjamin.rosche@gmail.com>
}
