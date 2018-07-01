## transfData: transfData from raw table to standard table. 
## dateFilter: subset data during given date interval.
## seasonFilter: subset data during given season for each year.
## prepareSeasonalFit: prepare 'seasonFilter'-like data to the
##                     data 'nllk_seasonal_parallel' needed.
## nllk_seasonal_parallel: nllk for seasonal analysis dataset
##                         with parallel function.
## fitMovResHan.seasonal.parallel: fit three states model for
##                                 seasonal analysis.


##' Transfer raw dataset to the standard dataset (seasonal analysis toolbox)
##'
##' Transfer the raw location dataset of animal to the standard dataset,
##' which is acceptable in this packages. The raw dataset contains at least
##' four components: 1. \code{t1}: data information. 2. \code{dt..hr.}: the
##' difference of time between two sample points. 3. \code{e1}: the GPS
##' coordinate of east-west. 4. \code{n1}: the GPS coordinate of north-south.
##' (These weird variable names are from the original GPS data. We will
##' change them in later version.)
##'
##' @param data The raw dataset.
##' @param dateFormat Charater string indicates the format of date variable.
##' @param roundValue Round GPS coordinate to \code{roundValue} with unit meter.
##' @param lengthUnit Charater string indicates the length unit of GPS coordinate,
##' which can be "m" or "km"(default). Usually, we recommend not change the
##' default setup of this parameter. Otherwise, numerical computation problem
##' will happen.
##'
##' @return A \code{data.frame} containing the following components, which is
##' standard format of dataset in this package:
##'
##' \itemize{
##'   \item date: tells us the date of collecting this sample point.
##'   \item cumTime: cumulative time line. The collection of this data starts
##' from time 0 in this time line. (Time unit is hours.)
##'   \item centerE: the centered east-west GPS coordinate with the center is
##' the starting point (when \code{cumTime[1]}).
##'   \item centerN: the centered north-south GPS coordinate with the center is
##' the starting point (when \code{cumTime[1]}).
##' }
##'
##' @seealso \code{\link{as.Date}} has parameter \code{format}, which is the
##' same as the parameter \code{dateFormat} in this function.
##' @author Chaoran Hu
##' @export
transfData <- function(data, dateFormat, roundValue, lengthUnit = "km") {
    stopifnot(lengthUnit %in% c("m", "km"))
    date <- as.Date(data$t1, format = dateFormat)
    cumTime <- cumsum(data$dt..hr.)
    centerE <- data$e1 - data$e1[1]
    centerN <- data$n1 - data$n1[1]
    result <- data.frame(date, cumTime, centerE, centerN)
    result[, 3:4] <- round(result[, 3:4] / roundValue) * roundValue
    if (lengthUnit == "km") {
        result[, 3:4] <- result[, 3:4] / 1000
        return(result)
    }
    result
}


## subsetting data during given date interval
## input: data:      data be filtered
##        startDate: the begin date of date interval
##        endDate:   the end date of date interval
## note:  both startDate and endDate follow format "YYYY-MM-DD"
## output: a data.frame contain only sample points from
##         given time interval.
dateFilter <- function(data, startDate, endDate) {
    startDate <- as.Date(startDate)
    endDate <- as.Date(endDate)
    time.label <- which((data[, 1] >= startDate) & (data[, 1] <= endDate))
    data[time.label, ]
}


##' Subsetting data during given season for each year (seasonal analysis toolbox)
##'
##' Return subsets of data from each year, which is in given
##' time interval. The time interval here is for one year and
##' defined by \code{startDate} and \code{endDate}.
##'
##' @param data The data be filtered, which has the same format
##' as the output from \code{\link{transfData}}.
##' @param startDate,endDate Start point and end point of
##' time interval during a year, which has the format "MM-DD".
##'
##' @return A \code{list} with each element is the one-year subset
##' of data, which is in the given time interval. For example,
##' suppose that we have data cross three years, this function
##' will return a list with three elements and each element
##' contains the subset of data of certain year.
##'
##' @author Chaoran Hu
##' @export
seasonFilter <- function(data, startDate, endDate) {
    year <- unique(format(data[, 1], "%Y"))
    n.year <- length(year)
    startMth <- as.numeric(substr(startDate, 1, 2))
    endMth <- as.numeric(substr(endDate, 1, 2))

    if (startMth > endMth) {
        result <- vector('list', n.year - 1)
        
        for(i in 1:(n.year-1)) {
            startday <- paste(year[i], "-", startDate, sep = "")
            endday <- paste(year[i+1], "-", endDate, sep = "")
            result[[i]] <- dateFilter(data, startday, endday)
        }
    } else {
        result <- vector('list', n.year)

        for(i in 1:n.year) {
            startday <- paste(year[i], "-", startDate, sep = "")
            endday <- paste(year[i], "-", endDate, sep = "")
            result[[i]] <- dateFilter(data, startday, endday)
        }
    }
    
    ## delete list elements with length 0
    zero.label <- sapply(result, nrow) != 0
    result[zero.label]
}


## delete the date column in the output of 'seasonFilter'.
## convert each element of data list to matrix
## generate diff of time and distance
prepareSeasonalFit <- function(data) {
    data <- lapply(data, function(x) x[, -1])
    data <- lapply(data, as.matrix)
    lapply(data, function(x) apply(x, 2, diff))
}




## The negative log-likelihood for seasonal analysis data.
## input: theta, integrControl, numThreads: are the same as
##               other nllk function.
##        data: list have the *similar* format as the output from
##              'seasonFilter' after 'prepareSeasonalFit'.
## output: negative log-likelihood of seasonal filtered data.
##' @import foreach
nllk_seasonal_parallel <- function(theta, data,
                          integrControl, numThreads) {
    n.year <- length(data)

    ## create parallel backend
    cl = parallel::makeCluster(numThreads); on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)

    i = 1 #Dummy line for Rstudio warnings
    result <- foreach(i = 1:n.year) %dopar% {
        nllk_fwd_ths(theta, data[[i]], integrControl)
    }

    sum(unlist(result))
}

## another parallel version which does not work
## in win-server
## nllk_seasonal_parallel <- function(theta, data,
##                                    integrControl, numThreads) {
##     n.year <- length(data)
##     result <- numeric(n.year)

##     grainSize <- lapply(data, function(x) {ceiling(nrow(x) / numThreads)})
##     grainSize <- unlist(grainSize)

##     for (i in 1:n.year) {
##         result[i] <- nllk_fwd_ths_parallel(theta, data[[i]], integrControl, grainSize[i])
##     }

##     sum(result)
## }



##' Fit MRH, MR, BMME Models for Seasonal Analysis (seasonal analysis toolbox)
##'
##' Fit a MRH, MR or BMME models, that is a wrapper of \code{\link{fitMovResHan}},
##' \code{\link{fitMovRes}} and \code{\link{fitBmme}} separately for
##' seasonal analysis. The special structure of data is required,
##' that is the same as the return from \code{\link{seasonFilter}}.
##' During seasonal analysis, the data is actually the subsets of whole
##' dataset which come from some certain period of each year.
##'
##' @param data The dataset should be fitted, which have the same
##' format as the output of \code{seasonFilter}.
##' @param start The initial value for optimization.
##' @param lower,upper,numThreads,integrControl are the same as \code{\link{fitMovResHan}}.
##' @param method,... are the parameters passed to \code{optim}.
##' @param likelihood,logtr,optim.control are the same as \code{\link{fitMovRes}}.
##'
##' @return A list of estimation result with the same format as
##' corresponding normal functions.
##' 
##' @seealso \code{\link{fitMovResHan}} for normal fit MRH model,
##' \code{\link{fitMovRes}} for normal fit MR model,
##' \code{\link{fitBmme}} for normal fit BMME model,
##' \code{\link{transfData}} for transferring original GPS data to standard
##' format, \code{\link{seasonFilter}} for subsetting given date interval
##' from transferred data.
##' 
##' @author Chaoran Hu
##' @export
fitMovResHan.seasonal <- function(data, start,
                                  lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
                                  upper = c(10, 10, 10, 10, 0.999),
                                  numThreads = NULL,
                                  integrControl = integr.control()) {
    if (is.null(numThreads)) numThreads <- RcppParallel::defaultNumThreads() * 3 / 4
    
    dinc <- prepareSeasonalFit(data)
    integrControl <- unlist(integrControl)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_seasonal_parallel,
                          data = dinc,
                          integrControl = integrControl,
                          numThreads = numThreads,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_COBYLA",
                                      "print_level" = 3,
                                      "maxeval" = 0))

    result <- list(estimate    =  fit[[18]],
                   loglik      = -fit[[17]],
                   convergence =  fit[[13]])
    result
}

