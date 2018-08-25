## transfData: transfData from raw table to standard table. 
## dateFilter: subset data during given date interval.
## seasonFilter: subset data during given season for each year.



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
## output: a data.frame contains additional column 'BATCH' indicates
##         which observation within given time interval.
dateFilter <- function(data, startDate, endDate) {
    startDate <- as.Date(startDate)
    endDate <- as.Date(endDate)
    time.label <- which((data[, 1] >= startDate) & (data[, 1] <= endDate))
    data[time.label, 'BATCH'] <- as.numeric(format(startDate, "%Y%m"))
    data
}


##' Subsetting data during given season for each year (seasonal analysis toolbox)
##'
##' Return subsets of data from each year, which is in given
##' time interval between \code{startDate} and \code{endDate}.
##'
##' @param data The data be filtered, which has the same format
##' as the output from \code{\link{transfData}}.
##' @param startDate,endDate Start point and end point of
##' time interval during a year, which has the format "MM-DD".
##'
##' @return A \code{data.frame} with inputted data and additional
##' column 'BATCH' indicates which subset of inputted data is located
##' within given time interval. In column 'BATCH', different integers
##' stands for different segments and 0 stands for outside given time
##' interval.
##'
##' @author Chaoran Hu
##' @export
seasonFilter <- function(data, startDate, endDate) {
    year <- unique(format(data[, 1], "%Y"))
    year <- as.numeric(year)
    n.year <- length(year)
    startMth <- as.numeric(substr(startDate, 1, 2))
    endMth <- as.numeric(substr(endDate, 1, 2))

    if (startMth > endMth) {

        
        minyear <- min(year)
        maxyear <- max(year)
        myyear <- c(minyear - 1, minyear:maxyear)
        n_myyear <- length(myyear)
        
        for(i in 1:n_myyear) {
            startday <- paste(myyear[i], "-", startDate, sep = "")
            endday <- paste(myyear[i] + 1, "-", endDate, sep = "")
            data <- dateFilter(data, startday, endday)
        }

        
    } else {

        
        for(i in 1:n.year) {
            startday <- paste(year[i], "-", startDate, sep = "")
            endday <- paste(year[i], "-", endDate, sep = "")
            data <- dateFilter(data, startday, endday)
        }

        
    }
    
    data$BATCH[is.na(data$BATCH)] <- 0
    
    data[, -1]
}
