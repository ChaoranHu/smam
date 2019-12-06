##' Fit a Moving-Resting-Handling Model with Embedded Brownian Motion
##'
##' Fit a Moving-Resting-Handling Model with Embedded Brownian Motion with
##' animal movement data at discretely observation times by maximizing
##' a full likelihood. Using \code{segment} to fit part of observations to
##' the model. A practical application of this feature is seasonal analysis.
##'
##' @param data a data.frame whose first column is the observation time, and other
#'     columns are location coordinates. If \code{segment} is not \code{NULL},
#'     additional column with the same name given by \code{segment} should be
#'     included. This additional column is used to indicate which part of
#'     observations shoule be used to fit model. The value of this column can
#'     be any integer with 0 means discarding this observation and non-0 means
#'     using this obversvation. Using different non-zero numbers indicate different
#'     segments. (See vignette for more details.)
##' @param start The initial value for optimization, in the order of rate
##'    of moving, rate of resting, rate of handling, volatility and switching
##'    probability.
##' @param segment character variable, name of the column which indicates segments,
##'    in the given \code{data.frame}. The default value, \code{NULL}, means using
##'    whole dataset to fit the model.
##' @param numThreads int, the number of threads allocated for parallel
##'    computation. The default setup is 3/4 available threads. If this parameter
##'    is less or equal to 1, the serial computation will be processed.
##' @param lower,upper Lower and upper bound for optimization.
##' @param integrControl Integration control vector includes rel.tol,
##'    abs.tol, and subdivisions.
##'
##' @return A list of estimation result with following components:
##' \item{estimate}{the estimated parameter vector}
##' \item{loglik}{maximized loglikelihood or composite loglikelihood
##' evaluated at the estimate}
##' \item{convergence}{convergence code from \code{nloptr}}
##'
##' @references
##' Pozdnyakov, V., Elbroch, L.M., Hu, C., Meyer, T., and Yan, J. (2018+)
##' On estimation for Brownian motion governed by telegraph process with
##' multiple off states. <arXiv:1806.00849>
##'
##' @seealso \code{\link{rMRH}} for simulation.
##'
##' @examples
##' \donttest{
##' ## slow work, may take several hours
##' set.seed(06269)
##' tgrid <- seq(0, 400, by = 8)
##' dat <- rMRH(tgrid, 4, 0.5, 0.1, 5, 0.8, 'm')
##' fitMRH(dat, c(4, 0.5, 0.1, 5, 0.8)) ## parallel process
##' fitMRH(dat, c(4, 0.5, 0.1, 5, 0.8), numThreads = -1) ## serial process
##'
##' ## fit part of dataset to the MRH model
##' batch <- c(rep(0, 10), rep(1, 7), rep(0, 10), rep(2, 10), rep(0, 14))
##' dat.segment <- cbind(dat, batch)
##' fit.segment <- fitMRH(dat.segment, start = c(4, 0.5, 0.1, 5, 0.8), segment = "batch")
##' head(dat.segment)
##' fit.segment
##' }
##' 
##' @author Chaoran Hu
##' @export
fitMRH <- function(data, start, segment = NULL,
                   numThreads = RcppParallel::defaultNumThreads() * 3 / 4,
                   lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
                   upper = c(   10,    10,    10,    10, 0.999),
                   integrControl = integr.control()) {
    if (is.null(segment)) {
        

        if (numThreads <= 1) { ## serial normal

            
            if (!is.matrix(data)) data <- as.matrix(data)
            dinc <- apply(data, 2, diff)
            integrControl <- unlist(integrControl)

            fit <- nloptr::nloptr(x0 = start, eval_f = nllk_fwd_ths,
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

            ## fit <- optim(par = start, fn = nllk_fwd_ths,
            ##              data = dinc,
            ##              integrControl = integrControl,
            ##              method = "L-BFGS-B",
            ##              lower  = lower,
            ##              upper = upper,
            ##              control = list(maxit =  1000))

            ## result <- list(estimate = fit$par,
            ##                loglik = -fit$value,
            ##                convergence = fit$convergence)
            
            return(result)

            
        } else { ## parallel normal

            
            result <- fitMRH_parallel(data, start, lower, upper, numThreads, integrControl)
            return(result)

            
        }

        
    } else { ## parallel seasonal

        
        if (numThreads < 1) numThreads <- 1
        result <- fitMRH_seasonal(data, segment, start,
                                  lower, upper, numThreads, integrControl)
        return(result)

        
    }
}

## internal function
fitMRH_parallel <- function(data, start, lower, upper,
                            numThreads, integrControl) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    grainSize <- ceiling(nrow(dinc) / numThreads)

    ## allocate threads
    RcppParallel::setThreadOptions(numThreads = numThreads)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_fwd_ths_parallel,
                          data = dinc,
                          integrControl = integrControl,
                          grainSize = grainSize,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_COBYLA",
                                      "print_level" = 3,
                                      "maxeval" = -5))

    result <- list(estimate    =  fit[[18]],
                   loglik      = -fit[[17]],
                   convergence =  fit[[14]])


    ## fit <- optim(par = start, fn = nllk_fwd_ths_parallel,
    ##              data = dinc,
    ##              integrControl = integrControl,
    ##              grainSize = grainSize,
    ##              method = "L-BFGS-B",
    ##              lower  = lower,
    ##              upper = upper,
    ##              control = list(maxit =  1000))

    ## result <- list(estimate = fit$par,
    ##                loglik = -fit$value,
    ##                convergence = fit$convergence)
    
    result
}



##### testing code for getting hessian matrix of nllk
hessMRH <- function(estimate, data, integrControl, numThreads) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    if (numThreads <= 1) {
        hess <- tryCatch(numDeriv::hessian(nllk_fwd_ths,
                                           x = estimate,
                                           data = dinc,
                                           integrControl = integrControl),
                         error = function(e) {
                             print(e)
                             matrix(NA, ncol = 5, nrow = 5)
                         })
    } else {
        grainSize <- ceiling(nrow(dinc) / numThreads)

        ## allocate threads
        RcppParallel::setThreadOptions(numThreads = numThreads)

        hess <- tryCatch(numDeriv::hessian(nllk_fwd_ths_parallel,
                                           x = estimate,
                                           data = dinc,
                                           integrControl = integrControl,
                                           grainSize = grainSize),
                         error = function(e) {
                             print(e)
                             matrix(NA, ncol = 5, nrow = 5)
                         })
    }
    
    hess
}


##### ending of testing code


##### do not export composite llk for MovResHan,
##### because it is more time consuming.

## composite nllk
## nllk_composite_parallel <- function(theta, data, groupSize,
##                                     integrControl, numThreads) {
##     nGroup <- nrow(data) %/% groupSize

##     ## create parallel backend
##     cl = parallel::makeCluster(numThreads); on.exit(parallel::stopCluster(cl))
##     doParallel::registerDoParallel(cl)

##     i = 1 #Dummy line for RStudio warnings
##     result <- foreach(i = 1:nGroup) %dopar% {
##         dataCart <- data[((i - 1) * groupSize + 1):(i * groupSize), ]
##         nllk_fwd_ths(theta, dataCart, integrControl)
##     }

##     sum(unlist(result))
## }

## fit base on composite nllk
## fitMovResHan.composite.parallel <- function(data, start, lower, upper,
##                                             groupSize,
##                                             numThreads = RcppParallel::defaultNumThreads() * 3 / 4,
##                                             integrControl = integr.control()) {
##     if (!is.matrix(data)) data <- as.matrix(data)
##     dinc <- apply(data, 2, diff)
##     integrControl <- unlist(integrControl)

##     fit <- nloptr::nloptr(x0 = start, eval_f = nllk_composite_parallel,
##                           data = dinc,
##                           integrControl = integrControl,
##                           groupSize = groupSize,
##                           numThreads = numThreads,
##                           lb = lower,
##                           ub = upper,
##                           opts = list("algorithm"   = "NLOPT_LN_COBYLA",
##                                       "print_level" = 3,
##                                       "maxeval" = 0))

##     fit
## }



