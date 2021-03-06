% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/movres.R
\name{approxNormalOrder}
\alias{approxNormalOrder}
\alias{approxNormalOrder2}
\title{Auxiliary for Preparing Discrete Distribution
used to approximating Standard Normal Distribution}
\usage{
approxNormalOrder(m)

approxNormalOrder2(m, width)
}
\arguments{
\item{m}{int, the number of order statistics used}

\item{width}{the width between two consecutive grid points.}
}
\value{
A numeric matrix with first column is support of discrete
distribution and second column is corresponding p.m.f..
}
\description{
Auxiliary for preparing discrete distribution used to
approximate standard normal. This function generates
order statistics of standard normal with same probability
assigned. Then, the discrete distribution is standardized
to variance one and mean zero.
}
\details{
This function use \code{EnvStats::evNormOrdStats} to get
the order statisics of standard normal distribution. The
same probability is assigned for each order statistics.
}
\seealso{
\code{EnvStats::evNormOrdStats} for order
statisics of standard normal. \code{\link{fitMRMEapprox}}
for fit MRME with approximated measurement error.
}
\author{
Chaoran Hu
}
