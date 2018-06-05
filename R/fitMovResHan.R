##' Fit a Moving-Resting-Handling Model with Embedded Brownian Motion
##'
##' Fit a Moving-Resting-Handling Model with Embedded Brownian Motion with
##' animal movement data at discretely observation times by maximizing
##' a full likelihood. The parallel code is provided as \code{fitMovResHan.parallel},
##' which improves the code speed significantly.
##'
##' @param data a \code{data.frame} whose first column is the observation
##' time, and other columns are location coordinates.
##' @param start The initial value for optimization, in the order of rate
##' of moving, rate of resting, rate of handling, volatility and switching
##' probability.
##' @param lower,upper Lower and upper bound for optimization.
##' @param integrControl Integration control vector includes rel.tol,
##' abs.tol, and subdivisions.
##' @param numThreads int, the number of threads allocated for parallel
##' computation. The default setup is 3/4 available threads.
##'
##' @return A list of estimation result, with estimator, log-likelihood and
##' convergence code from \code{nloptr}.
##'
##' @references
##' Pozdnyakov, V., Elbroch, L.M., Hu, C., Meyer, T., and Yan, J. (2018+)
##' On estimation for Brownian motion governed by telegraph process with
##' multiple off states. <arXiv:1806.00849>
##'
##' @seealso \code{\link{rMovResHan}} for simulation.
##'
##' @examples
##' \donttest{
##' ## slow work, may take several hours
##' set.seed(06269)
##' tgrid <- seq(0, 400, by = 8)
##' dat <- rMovResHan(tgrid, 4, 0.5, 0.1, 5, 0.8, 'm')
##' fitMovResHan(dat, c(4, 0.5, 0.1, 5, 0.8))
##' fitMovResHan.parallel(dat, c(4, 0.5, 0.1, 5, 0.8))
##' }
##' 
##' @author Chaoran Hu
##' @export
fitMovResHan <- function(data, start,
                         lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
                         upper = c(   10,    10,    10,    10, 0.999),
                         integrControl = integr.control()) {
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
                                      "maxeval" = 0))

    result <- list(estimate    =  fit[[18]],
                   loglik      = -fit[[17]],
                   convergence =  fit[[13]])
    result
}

##' @rdname fitMovResHan
##' @export
fitMovResHan.parallel <- function(data, start,
                                  lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
                                  upper = c(   10,    10,    10,    10, 0.999),
                                  numThreads = NULL,
                                  integrControl = integr.control()) {
    if (is.null(numThreads)) numThreads <- RcppParallel::defaultNumThreads() * 3 / 4
    
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
                                      "maxeval" = 0))

    result <- list(estimate    =  fit[[18]],
                   loglik      = -fit[[17]],
                   convergence =  fit[[13]])
    result
}



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



