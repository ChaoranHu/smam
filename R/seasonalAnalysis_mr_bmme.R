## nllk_bmme_seasonal:    nllk for seasonal analysis with bmme
##                        process.
## ncllk_m1_inc_seasonal: composite nllk for seasonal analysis.
## nllk_inc_seasonal:     forward algorithm nllk for seasonal analysis.
## fitBmme.seasonal:      fit bmme model for seasonal analysis.
## fitMovRes.seasonal:    fit a moving-resting model for seasonal analysis.



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

## The composite and forward algorithm nllk of two-states model
## for seasonal analysis data.
## input:
##        theta: vector of (lambda1, lambda0, sigma)
##        data:  list have the *similar* format as the output from
##               'seasonFilter' after 'prepareSeasonalFit'
##        integrControl, logtr: see ncllk_m1_inc
## output:
##        negative log-likelihood of seasonal filtered data
ncllk_m1_inc_seasonal <- function(theta, data,
                                  integrControl, logtr) {
    n.year <- length(data)
    result <- lapply(data, ncllk_m1_inc,
                     theta = theta, integrControl = integrControl,
                     logtr = logtr)
    sum(unlist(result))
}

nllk_inc_seasonal <- function(theta, data,
                              integrControl, logtr) {
    n.year <- length(data)
    result <- lapply(data, nllk_inc,
                     theta = theta, integrControl = integrControl,
                     logtr = logtr)
    sum(unlist(result))
}

## obtain initial value for sigma and delta by method of moment
## (the wrapper of 'bmme.start' for seasonal analysis data.)
## input:
##      dat: list with the same format as the output from 'seasonalFilter'
## output:
##      a vector containing rough initial value of sigma 
bmme.start.seasonal <- function(dat) {
    dif <- prepareSeasonalFit(dat)
    dim <- ncol(dif[[1]]) - 1

    numerator <- sum(unlist(lapply(dif, function(x) {sum(x[,-1]^2)})))
    denominator <- sum(unlist(lapply(dif, function(x) {sum(2 + x[,1])}))) * dim
    st <- sqrt(numerator / denominator)
    c(st, st)
}


##' @importFrom stats optim
##' @importFrom methods is
##' @rdname fitMovResHun.seasonal
##' @export
fitBmme.seasonal <- function(data, start = NULL, method = "Nelder-Mead", ...) {
    if (is.null(start)) start <- bmme.start.seasonal(data)
    dinc <- prepareSeasonalFit(data)
    fit <- optim(start, nllk_bmme_seasonal, data = dinc, method=method, ...)
    
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
movres.start.seasonal <- function(dat) {
    dat <- lapply(dat, function(x) {x[, -1]})
    dinc <- lapply(dat, function(x) {apply(x, 2, diff)})
    dinc <- glue_list(dinc)
    s1 <- coef(lm(dinc[, 2]^2 ~ dinc[, 1] - 1))
    s2 <- coef(lm(dinc[, 3]^2 ~ dinc[, 1] - 1))
    ss <- sqrt(mean(c(s1, s2)))
    c(0.5, 0.5, ss)
}


##' @importFrom stats optim
##' @importFrom methods is
##' @rdname fitMovResHun.seasonal
##' @export
fitMovRes.seasonal <- function(data, start = NULL, likelihood = c("full", "composite"),
                               logtr = FALSE,
                               method = "Nelder-Mead",
                               optim.control = list(),
                               integrControl = integr.control()) {
    if (is.null(start)) start <- movres.start.seasonal(data)
    dinc <- prepareSeasonalFit(data)
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
