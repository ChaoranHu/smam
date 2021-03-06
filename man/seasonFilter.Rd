% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seasonalTools.R
\name{seasonFilter}
\alias{seasonFilter}
\title{Subsetting data during given season for each year (seasonal analysis toolbox)}
\usage{
seasonFilter(data, startDate, endDate)
}
\arguments{
\item{data}{The data be filtered, which has the same format
as the output from \code{\link{transfData}}.}

\item{startDate, endDate}{Start point and end point of
time interval during a year, which has the format "MM-DD".}
}
\value{
A \code{data.frame} with inputted data and additional
column 'BATCH' indicates which subset of inputted data is located
within given time interval. In column 'BATCH', different integers
stands for different segments and 0 stands for outside given time
interval.
}
\description{
Return subsets of data from each year, which is in given
time interval between \code{startDate} and \code{endDate}.
}
\author{
Chaoran Hu
}
