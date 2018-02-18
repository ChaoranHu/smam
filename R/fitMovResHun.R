##' Fit Three States Model (forward algrithm nllk)
##'
##' Fit Model using different optimation approaches. The fit1 bases on
##' Nelder-Mead. The fit5 bases on COBYLA. The fit6 bases on BOBYQA.
##'
##' @param data The dataset should be fitted.
##' @param start The initial value for optimization.
##' @param optim.control The control parameter for R function optim.
##' @param integrControl Integration control vector includes rel.tol,
##' abs.tol, and subdivisions.
##' @param grainSize Minimum chunk size for parallelization.
##' @param numThreads Allocate the number of threads.
##' @param lower Lower bound for optimization.
##' @param upper Upper bound for optimization.
##'
##' @return A list of estimation result.
##' 
##' @author Chaoran Hu
##' @export
fitMovResHun1 <- function(data, start,
                          optim.control = list(),
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    objfun <- function(theta) {
        if (theta[1] > 0 & theta[2] > 0 & theta[3] > 0 &
            theta[4] > 0 & theta[5] > 0 & theta[5] < 1) {
            return(nllk_fwd_ths(theta, dinc, integrControl))
        } else {
            return(NA)
        }
    }
    
    fit <- stats::optim(par = start, fn = objfun,
                        method = "Nelder-Mead",
                        control = optim.control)

    fit
}

##' @rdname fitMovResHun1
##' @export
fitMovResHun5 <- function(data, start, lower, upper,
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

    fit
}

##' @rdname fitMovResHun1
##' @export
fitMovResHun.parallel <- function(data, start, lower, upper,
                                  grainSize,
                                  numThreads = RcppParallel::defaultNumThreads() * 3 / 4,
                                  integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

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

    fit
}

##' @rdname fitMovResHun1
##' @export
fitMovResHun6 <- function(data, start, lower, upper,
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_fwd_ths,
                          data = dinc,
                          integrControl = integrControl,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_BOBYQA",
                                      "print_level" = 3,
                                      "maxeval" = 0))

    fit
}


##' Fit Three States Model (composite nllk with parallel feature)
##'
##' Fit Model using COBYLA optimation approaches for composite negetive
##' log likelihood.
##'
##' @param data The dataset should be fitted.
##' @param start The initial value for optimization.
##' @param lower Lower bound for optimization.
##' @param upper Upper bound for optimization.
##' @param groupSize int size of each group for composite likelihood.
##' @param integrControl Integration control vector includes rel.tol,
##' abs.tol, and subdivisions.
##' @param numThreads Allocate the number of threads.
##'
##' @return A list of estimation result.
##' 
##' @author Chaoran Hu
##' @export
fitMovResHun.composite.parallel <- function(data, start, lower, upper,
                                            groupSize,
                                            numThreads = RcppParallel::defaultNumThreads() * 3 / 4,
                                            integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_composite_parallel,
                          data = dinc,
                          integrControl = integrControl,
                          groupSize = groupSize,
                          numThreads = numThreads,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_COBYLA",
                                      "print_level" = 3,
                                      "maxeval" = 0))

    fit
}




#' Auxiliary for Controlling Numerical Integration
#'
#' Auxiliary function for the numerical integration used in the
#' likelihood and composite likelihood functions. Typically only
#' used internally by 'fitMovRes'. 
#'
#' @param rel.tol relative accuracy requested.
#' @param abs.tol absolute accuracy requested.
#' @param subdivisions the maximum number of subintervals.
#'
#' @details
#' The arguments are the same as \code{integrate}, but passed
#' down to the C API of Rdqags used by \code{integrate}.
#'
#' @return
#' A list with components named as the arguments.
#'
#' @export
integr.control <- function(rel.tol = .Machine$double.eps^.25,
                           abs.tol = rel.tol, subdivisions = 100L) {
    if (!is.numeric(rel.tol) || rel.tol <= 0) 
        stop("value of 'rel.tol' must be > 0")
    if (!is.numeric(abs.tol) || abs.tol <= 0) 
        stop("value of 'abs.tol' must be > 0")
    if (!is.numeric(subdivisions) || subdivisions <= 0) 
        stop("maximum number of subintervals must be > 0")
    list(rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions)
}



.onUnload <- function (libpath) {
  library.dynam.unload("thmam", libpath)
}
