## nllk_bmme_seasonal:          nllk for seasonal analysis with bmme process.
## ncllk_m1_inc_seasonal:       composite nllk for seasonal analysis for MR model.
## nllk_inc_seasonal:           full nllk for seasonal analysis for MR model.
## nllk_mrme_seasonal:          nllk for seasonal analysis for MRME model.
## nllk_mrme_approx_seasonal:   nllk for seasonal analysis for MRME approx model.
## nllk_seasonal_parallel:      full nllk for seasonal analysis for MRH model.
## fitBMME_seasonal:            fit bmme model for seasonal analysis.
## fitMR_seasonal:              fit a moving-resting model for seasonal analysis.
## fitMRME_seasonal:            fit a moving-resting model with measurement error for seasonal analysis.
## fitMRMEapprox_seasonal:      fit approximate moving-resting model with measurement error for seasonal analyis.
## fitMRH_seasonal:             fit a moving-resting-handling model for seasonal analysis.







## transfer data to list according to given segments
## data: standard data.frame with first col is time,
##       other cols are coordinates, one of them is
##       indicate segments.
## segment: character variable indicate which col
##       should be used for indicate segments.
## return: a list with each element is a segment.
seg2list <- function(data, segment="BATCH") {
    if (segment %in% names(data) != TRUE) {
        stop("Cannot find segment variable in data.frame.")
    }

    seg.ncol <- which(names(data) == segment)
    seg.col <- data[, seg.ncol]
    
    if (length(unique(seg.col)) == 1) {
        if (unique(seg.col) != 0) {
            return(list(data))
        } else {
            return(NA)
        }
    }
    
    uni.seg <- unique(seg.col)
    uni.seg <- uni.seg[uni.seg != 0]
    n.seg <- length(uni.seg)
    
    result <- vector("list", n.seg)
    
    for (i in seq_len(n.seg)) {
        result[[i]] <- data[seg.col == uni.seg[i], ]
    }
        
    res_len <- sapply(result, nrow)
    result[which(res_len > 2)]
}


## delete the date column in the output of 'seasonFilter'.
## convert each element of data list to matrix
## generate diff of time and distance
prepareSeasonalFit <- function(data, segment) {
    seg.col <- which(names(data[[1]]) == segment)
    data <- lapply(data, function(x) x[, -seg.col])
    data <- lapply(data, as.matrix)
    lapply(data, function(x) apply(x, 2, diff))
}



## The negative log-likelihood of bmme for seasonal analysis data.
## input:
##        param: vector of (sigma, delta)
##        data:  list have the *similar* format as the output from
##               'seasonFilter' after 'prepareSeasonalFit'
## output:
##        negative log-likelihood of seasonal filtered data
nllk_bmme_seasonal <- function(param, data) {
    n.year <- length(data)
    result <- lapply(data, nllk.bmme, param = param)
    sum(unlist(result))
}

## The nllk of moving-resting model for seasonal analysis data.
## input:
##        theta: vector of (lambda1, lambda0, sigma)
##        data:  list have the *similar* format as the output from
##               'seasonFilter' after 'prepareSeasonalFit'
##        integrControl, logtr: see ncllk_m1_inc
## output:
##        negative log-likelihood of seasonal filtered data

## composite nllk
ncllk_m1_inc_seasonal <- function(theta, data,
                                  integrControl, logtr) {
    n.year <- length(data)
    result <- lapply(data, ncllk_m1_inc,
                     theta = theta, integrControl = integrControl,
                     logtr = logtr)
    sum(unlist(result))
}

## forward nllk
nllk_inc_seasonal <- function(theta, data,
                              integrControl, logtr) {
    n.year <- length(data)
    result <- lapply(data, nllk_inc,
                     theta = theta, integrControl = integrControl,
                     logtr = logtr)
    sum(unlist(result))
}


## The nllk of moving-resting model with measurement error
## for seasonal analysis data.
## input:
##        theta: vector of (lambda1, lambda0, sigma, sig_err)
##        data:  list have the *similar* format as the output from
##               'seasonFilter' after 'prepareSeasonalFit'
##        integrControl, logtr: see ncllk_m1_inc
## output:
##        negative log-likelihood of seasonal filtered data
nllk_mrme_seasonal <- function(theta, data, integrControl) {
    n.year <- length(data)
    result <- lapply(data, nllk_mrme,
                     theta = theta, integrControl = integrControl)
    sum(unlist(result))
}
## naive composite llk for MRME
nllk_mrme_naive_cmp_seasonal <- function(theta, data, integrControl) {
    n.year <- length(data)
    result <- lapply(data, nllk_mrme_naive_cmp,
                     theta = theta, integrControl = integrControl)
    sum(unlist(result))
}


## The nllk of moving-resting model with approximate measurement
## error for seasonal analysis data.
## input:
##        theta: vector of (lambda1, lambda0, sigma, sig_err)
##        data:  list have the *similar* format as the output from
##               'seasonFilter' after 'prepareSeasonalFit'
##        integrControl: see ncllk_m1_inc
##        approx_norm_even, approx_norm_odd: see comment in MRME_approx.cpp
## output:
##        negative log-likelihood of seasonal filtered data
nllk_mrme_approx_seasonal <- function(theta, data, integrControl,
                                      approx_norm_even, approx_norm_odd) {
    n.year <- length(data)
    result <- lapply(data, nllk_mrme_approx,
                     theta = theta, integrControl = integrControl,
                     approx_norm_even = approx_norm_even, approx_norm_odd = approx_norm_odd)
    sum(unlist(result))
}







####################################
### BMME Model seasonal analysis ###
####################################

## obtain initial value for sigma and delta by method of moment
## (the wrapper of 'bmme.start' for seasonal analysis data.)
## input:
##      dat: list with the same format as the output from 'seasonalFilter'
## output:
##      a vector containing rough initial value of sigma 
bmme.start.seasonal <- function(dat, segment) {
    dif <- prepareSeasonalFit(dat, segment)
    dim <- ncol(dif[[1]]) - 1

    numerator <- sum(unlist(lapply(dif, function(x) {sum(x[,-1]^2)})))
    denominator <- sum(unlist(lapply(dif, function(x) {sum(2 + x[,1])}))) * dim
    st <- sqrt(numerator / denominator)
    c(st, st)
}

## internal function for fitBMME seasonal
##' @importFrom stats optim
##' @importFrom methods is
fitBMME_seasonal <- function(data, segment, start, method, optim.control) {
    data <- seg2list(data, segment)
    if (is.null(start)) start <- bmme.start.seasonal(data, segment)
    dinc <- prepareSeasonalFit(data, segment)
    fit <- optim(start, nllk_bmme_seasonal, data = dinc, method=method, control = optim.control)
    
    ## get variance estimate
    varest <- matrix(NA_real_, 2, 2)
    estimate <- fit$par
    ## always do this without log transformation
    hess <- tryCatch(numDeriv::hessian(nllk_bmme_seasonal, estimate, data = dinc),
                     error = function(e) e)
    if (!is(hess, "error")) varest <- solve(hess)
    ## not -hess because objfun is negative llk already
    
    list(estimate    = estimate,
         varest      = varest,
         loglik      = -fit$value,
         convergence = fit$convergence)
}


##############################################
### Moving-Resting Model seasonal analysis ###
##############################################



## glue all elements in a list together as a matrix
## input:
##        x: list, whose elements are same format data.fram
##           or matrix.
## output:
##        unlist the input as a matrix.
glue_list <- function(x) {
    n <- length(x)
    m <- sapply(x, nrow)
    result <- matrix(NA_real_, sum(m), 3)
    m <- c(0, cumsum(m))
    for(i in 1:n) {
        result[(m[i]+1):m[i+1],] <- x[[i]]
    }
    result
}

## obtain initial value for Mov-Res model for seasonal analysis
## input:
##      dat: list with the same format as the output from 'seasonalFilter'
## output:
##      a vector containing rough initial value of sigma
##' @importFrom stats coef lm
movres.start.seasonal <- function(dat, segment) {
    dat <- prepareSeasonalFit(dat, segment)
    dinc <- glue_list(dinc)
    s1 <- coef(lm(dinc[, 2]^2 ~ dinc[, 1] - 1))
    s2 <- coef(lm(dinc[, 3]^2 ~ dinc[, 1] - 1))
    ss <- sqrt(mean(c(s1, s2)))
    c(0.5, 0.5, ss)
}



## internal function for fitMR seasonal
##' @importFrom stats optim
##' @importFrom methods is
fitMR_seasonal <- function(data, segment, start, likelihood,
                           logtr, method, optim.control, integrControl) {
    data <- seg2list(data, segment)
    if (is.null(start)) start <- movres.start.seasonal(data, segment)
    dinc <- prepareSeasonalFit(data, segment)
    objfun <- switch(likelihood,
                     composite = ncllk_m1_inc_seasonal,
                     full = nllk_inc_seasonal,
                     stop("Not valid likelihood type.")
                     )
    integrControl <- unlist(integrControl)
    fit <- optim(start, objfun, data = dinc, method = method,
                 control = optim.control, 
                 integrControl = integrControl, 
                 logtr = logtr)
    ## get variance estimate
    varest <- matrix(NA_real_, 3, 3)
    estimate <- if (logtr) exp(fit$par) else fit$par
    if (likelihood == "full") {
        ## always do this without log transformation
        hess <- tryCatch(numDeriv::hessian(
            objfun, estimate, data = dinc,
            integrControl = integrControl, logtr = FALSE),
                         error = function(e) e)
        if (!is(hess, "error")) varest <- solve(hess)
        ## not -hess because objfun is negative llk already
    }
    
    list(estimate    = estimate,
         varest      = varest,
         loglik      = -fit$value,
         convergence = fit$convergence,
         likelihood  = likelihood)
}


#####################################################################
### Moving-Resting with Measurement Error Model seasonal analysis ###
#####################################################################
## internal function for fitMR seasonal
##' @importFrom stats optim
##' @importFrom methods is
fitMRME_seasonal <- function(data, segment, start,
                             lower, upper,
                             print_level,
                             #method, optim.control,
                             integrControl) {
    data <- seg2list(data, segment)
    if (is.null(start)) start <- movres.start.seasonal(data, segment)
    dinc <- prepareSeasonalFit(data, segment)
    integrControl <- unlist(integrControl)


    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_mrme_seasonal,
                          data = dinc,
                          integrControl = integrControl,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_COBYLA",
                                      "print_level" = print_level,
                                      "maxeval" = -5))

    result <- list(estimate    =  fit[[18]],
                   loglik      = -fit[[17]],
                   convergence =  fit[[14]])

    return(result)
    
    ## fit <- optim(start, nllk_mrme_seasonal, data = dinc, method = method,
    ##              control = optim.control, integrControl = integrControl)
    
    ## estimate <- fit$par
    
    ## list(estimate    = estimate,
    ##      loglik      = -fit$value,
    ##      convergence = fit$convergence)
}


fitMRME_naive_seasonal <- function(data, segment, start,
                                   lower, upper,
                                   #method, optim.control,
                                   integrControl) {
    data <- seg2list(data, segment)
    if (is.null(start)) start <- movres.start.seasonal(data, segment)
    dinc <- prepareSeasonalFit(data, segment)
    integrControl <- unlist(integrControl)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_mrme_naive_cmp_seasonal,
                          data = dinc,
                          integrControl = integrControl,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_COBYLA",
                                      "print_level" = 3,
                                      "maxeval" = -5))

    result <- list(estimate    =  fit[[18]],
                   loglik      = -fit[[17]],
                   convergence =  fit[[14]])

    return(result)
    
    ## fit <- optim(start, nllk_mrme_naive_cmp_seasonal, data = dinc, method = method,
    ##              control = optim.control, integrControl = integrControl)
    
    ## estimate <- fit$par
    
    ## list(estimate    = estimate,
    ##      loglik      = -fit$value,
    ##      convergence = fit$convergence)
}


fitMRMEapprox_seasonal <- function(data, segment, start,
                                   approx_norm_even, approx_norm_odd,
                                   method, optim.control, integrControl) {
    data <- seg2list(data, segment)
    if (is.null(start)) start <- movres.start.seasonal(data, segment)
    dinc <- prepareSeasonalFit(data, segment)
    integrControl <- unlist(integrControl)
    
    fit <- optim(start, nllk_mrme_approx_seasonal, data = dinc, method = method,
                 control = optim.control, integrControl = integrControl,
                 approx_norm_even = approx_norm_even, approx_norm_odd = approx_norm_odd)
    
    
    list(estimate    = fit$par,
         loglik      = -fit$value,
         convergence = fit$convergence)
}






#######################################################
### Moving-Resting-Handling Model Seasonal Analysis ###
#######################################################

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


## internal function for fitMRH seasonal
fitMRH_seasonal <- function(data, segment, start,
                            lower, upper,
                            numThreads, integrControl, print_level) {
    data <- seg2list(data, segment)
    dinc <- prepareSeasonalFit(data, segment)

    integrControl <- unlist(integrControl)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_seasonal_parallel,
                          data = dinc,
                          integrControl = integrControl,
                          numThreads = numThreads,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_COBYLA",
                                      "print_level" = print_level,
                                      "maxeval" = -5))

    result <- list(estimate    =  fit[[18]],
                   loglik      = -fit[[17]],
                   convergence =  fit[[14]])
    result
}

