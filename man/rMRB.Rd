% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rMRB.R
\name{rMRB}
\alias{rMRB}
\title{Sampling from a Moving-Resting bridge}
\usage{
rMRB(time, start_pt, end_pt, lamM, lamR, sigma, s0)
}
\arguments{
\item{time}{time points at which observations are to be simulated}

\item{start_pt}{the start point location of Brownian bridge
(numeric for one dimension, vector for multiple dimension)}

\item{end_pt}{the end point location of Brownian brige
(numeric for one dimension, vector for multiple dimension)}

\item{lamM}{rate parameter of the exponential duration while moving}

\item{lamR}{rate parameter of the exponential duration while resting}

\item{sigma}{volatility parameter of the Brownian motion while moving}

\item{s0}{the state at time 0, must be one of "m" or "r", for moving and
resting, respectively}
}
\value{
A \code{data.frame} whose first column is the time points and whose
other columns are coordinates of the locations.
}
\description{
A moving-resting process consists of two states: moving and resting.
The transition between the two states is modeled by an alternating
renewal process, with exponentially distributed duration. An animal
stays at the same location while resting, and moves according to a
Brownian motion while moving. A moving-resting bridge is an extension
of Brownian bridge that uses moving-resting process to bridge given
starting point and ending point.
}
\examples{
tgrid <- seq(0, 10, length=1001)
## make it irregularly spaced
tgrid <- sort(sample(tgrid, 800))

## a 2-dim moving-resting bridge starting from c(0, 0)
## and ending at c(10, -10)
dat <- rMRB(tgrid, start_pt = c(0, 0), end_pt = c(10, -10),
            lamM = 1, lamR = 1, sigma = 1, "m")
par(mfrow = c(2, 1))
plot(dat[,1], dat[,2], xlab="t", ylab="X(t)", type='l')
plot(dat[,1], dat[,3], xlab="t", ylab="X(t)", type='l')
par(mfrow = c(1, 1))

}
