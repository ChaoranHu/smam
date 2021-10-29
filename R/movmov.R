## simulation of breaking time points bbs
## for a 2 state telegraph process
sim1mm.times.bbz <- function(s, lam1, lam2) {
    tsum <- 0
    tt <- NULL
    state <- NULL
    while (TRUE) {
        tnew <- rexp(1, lam1) ## starting from state 1
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        state <- c(state, 0)
        if (tsum > s) break
        tnew <- rexp(1, lam2)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        state <- c(state, 1)
        if (tsum > s) break
    }
    cbind(tt, state)
}

## simulation of a realization given breaking times 
sim1mm.bbz <- function(s, sigma1, sigma2, time, brtimes, t0moving=TRUE) {
#### time: time points in [0, s]
    tt <- sort(unique(c(time, brtimes)))
    nt <- length(tt)
    nb <- length(brtimes) 
    x <- rep(NA, nt)
    x[1] <- 0
    status <- as.integer(t0moving)
    tend <- brtimes[1]
    j <- 1
    for (i in 2:nt) {
        if (tt[i] <= tend) { ## status unchanged
            if (status == 1)   ## moving1
              x[i] <- x[i - 1] + rnorm(1, sd = sigma1 * sqrt(tt[i] - tt[i - 1]))
            else               ## moving2
              x[i] <- x[i - 1] + rnorm(1, sd = sigma2 * sqrt(tt[i] - tt[i - 1]))
        }
        if (tt[i] == tend) {
            status <- 1 - status ## switch status
            j <- j + 1
            tend <- brtimes[j]
        }
    }
    x[tt %in% time]
}


## simulation a moving-moving path given a grid time

#' Sampling from a Moving-Moving Process with 2 Embedded Brownian Motion
#'
#' A moving-moving process consists of two states: moving (large) and moving (small).
#' The transition between the two states is modeled by an alternating
#' renewal process, with exponentially distributed duration. An animal
#' moves according to two Brownian motions with different volatility parameters.
#'
#' @param time time points at which observations are to be simulated
#' @param lamM1 rate parameter of the exponential duration while moving1
#' @param lamM2 rate parameter of the exponential duration while moving2
#' @param sigma1 volatility parameter of the Brownian motion while moving1
#' @param sigma2 volatility parameter of the Brownian motion while moving2
#' @param s0 the state at time 0, must be one of "m1" or "m2", for moving1 and
#' moving2, respectively
#' @param dim (integer) dimension of the Brownian motion
#'
#' @return
#' A \code{data.frame} whose first column is the time points and whose
#' other columns are coordinates of the locations.
#' @references
#' Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakov, V.,
#' Williams, S., and Meyer, T. (2014) A moving-resting process with an
#' embedded Brownian motion for animal movements.
#' Population Ecology. 56(2): 401--415.
#'
#' Pozdnyakov, V., Elbroch, L., Labarga, A., Meyer, T., and Yan, J.
#' (2017) Discretely observed Brownian motion governed by telegraph
#' process: estimation. Methodology and Computing in Applied Probability.
#' doi:10.1007/s11009-017-9547-6.
#' 
#' @examples
#' tgrid <- seq(0, 100, length=100)
#' 
#' dat <- rMM(tgrid, 1, 0.1, 1, 0.1, "m1")
#' plot(dat[,1], dat[,2], xlab="t", ylab="X(t)", type='l')
#' 
#' @export

rMM <- function(time, lamM1, lamM2, sigma1, sigma2, s0, dim = 2) {
    time <- time - time[1]
    stopifnot(s0 %in% c("m1", "m2"))
    t0moving <- (s0 == "m1")
    lam1 <- if (t0moving) lamM1 else lamM2
    lam2 <- if (t0moving) lamM2 else lamM1
    
    timeIND <- 0
    if (length(time) == 1) {
        timeIND <- 1
        time <- c(0, time)
    }
    tmax <- time[length(time)]
    brtimes <- sim1mm.times.bbz(tmax, lam1, lam2)
    coord <- replicate(dim, sim1mm.bbz(tmax, sigma1, sigma2, time, brtimes[, 1], t0moving))

    stateresult <- simmr.state(time, brtimes, t0moving = t0moving)

    if (timeIND == 1) {
        return(data.frame(time = time, coord)[-1, ])
    }
    
    data.frame(time = time, coord)
}









#' Fit a Moving-Moving Model with 2 Embedded Brownian Motion
#'
#' Fit a Moving-Moving Model with 2 Embedded Brownian Motion with animal
#' movement data at discretely observation times by maximizing a full
#' likelihood constructed from the marginal density of increment.
#' 'estVarMM' uses parametric bootstrap to obtain variance matrix
#' of estimators from 'fitMM'.
#'
#' @param data a data.frame whose first column is the observation time, and other
#'     columns are location coordinates.
#' @param start starting value of the model, a vector of four components
#'     in the order of rate for moving1, rate for moving2,
#'     and volatility1(larger), volatility2(smaller).
#' @param logtr logical, if TRUE parameters are estimated on the log scale.
#' @param method the method argument to feed \code{optim}.
#' @param optim.control a list of control to be passed to \code{optim}.
#' @param integrControl a list of control parameters for the \code{integrate}
#'     function: rel.tol, abs.tol, subdivision.
#' 
#' @return
#' a list of the following components:
#' \item{estimate}{the esimated parameter vector}
#' \item{loglik}{maximized loglikelihood or composite loglikelihood
#' evaluated at the estimate}
#' \item{convergence}{convergence code from \code{optim}}
#' @references
#' Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakov, V.,
#' Williams, S., and Meyer, T. (2014) A moving-resting process with an
#' embedded Brownian motion for animal movements.
#' Population Ecology. 56(2): 401--415.
#'
#' Pozdnyakov, V., Elbroch, L., Labarga, A., Meyer, T., and Yan, J.
#' (2017) Discretely observed Brownian motion governed by telegraph
#' process: estimation. Methodology and Computing in Applied Probability.
#' doi:10.1007/s11009-017-9547-6.
#'
#' @examples
#' \dontrun{
#' ## time consuming example
#' tgrid <- seq(0, 100, length=100)
#' set.seed(123)
#' dat <- rMM(tgrid, 1, 0.1, 1, 0.1, "m1")
#'
#' ## fit whole dataset to the MR model
#' fit <- fitMM(dat, start=c(1, 0.1, 1, 0.1))
#' fit
#'
#' var <- estVarMM(fit$estimate, dat, nBS = 10, numThreads = 6)
#' var
#' }
#' @export
fitMM <- function(data, start,
                  logtr = FALSE,
                  method = "Nelder-Mead",
                  optim.control = list(),
                  integrControl = integr.control()) {
    
    ## normal process
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    
    fit <- optim(start, nllk_inc_mm, data = dinc, method = method,
                 control = optim.control, 
                 integrControl = integrControl, 
                 logtr = logtr)
    ## ## get variance estimate
    ## varest <- matrix(NA_real_, 3, 3)
    estimate <- if (logtr) exp(fit$par) else fit$par
    ## if (likelihood == "full") {
    ##     ## always do this without log transformation
    ##     hess <- tryCatch(numDeriv::hessian(
    ##                                    objfun, estimate, data = dinc,
    ##                                    integrControl = integrControl, logtr = FALSE),
    ##                      error = function(e) e)
    ##     if (!is(hess, "error")) varest <- solve(hess)
    ##     ## not -hess because objfun is negative llk already
    ## }
    
    result <- list(estimate    = estimate,
                #varest      = varest,
                loglik      = -fit$value,
                convergence = fit$convergence,
                data = data)
    attr(result, "class") <- "smam_mm"
    return(result)

}


#' @param est_theta estimators of MRME model
#' @param data data used to process estimation
#' @param nBS number of bootstrap.
#' @param detailBS whether or not output estimation results during bootstrap,
#'     which can be used to generate bootstrap CI.
#' @param numThreads the number of threads for parallel computation. If its value
#'     is greater than 1, then parallel computation will be processed. Otherwise,
#'     serial computation will be processed.
#' @param integrControl a list of control parameters for the \code{integrate}
#'     function: rel.tol, abs.tol, subdivision.
#' @rdname fitMM
#' @export
estVarMM <- function(est_theta, data, nBS, detailBS = FALSE,
                     numThreads = 1,
                     integrControl = integr.control()) {
    tgrid <- data[, 1]
    dim <- ncol(data) - 1
    
    lam1 <- est_theta[1]
    lam2 <- est_theta[2]
    sig1 <- est_theta[3]
    sig2 <- est_theta[4]

    p_1 <- 1/lam1/(1/lam1 + 1/lam2)
    p_2 <- 1 - p_1

    if (numThreads <= 1) {
        
        result <- matrix(NA, ncol = 4, nrow = nBS)

        for (i in seq_len(nBS)) {
            start_state <- sample(c("m1", "m2"), size = 1, prob = c(p_1, p_2))
            datBS <- rMM(tgrid, lam1, lam2, sig1, sig2, s0 = start_state, dim = dim)
            result[i, ] <- fitMM(datBS, start = est_theta,
                                 integrControl = integrControl)$estimate
        }

    } else {

        ## create parallel backend
        cl = parallel::makeCluster(numThreads)
        on.exit(close(pb), add = TRUE)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        doSNOW::registerDoSNOW(cl)
        pb <- utils::txtProgressBar(max = nBS, style = 3)
        progress <- function(n) utils::setTxtProgressBar(pb, n)
        opts <- list(progress = progress)

        i = 1
        result <- foreach(i = 1:nBS, .combine = rbind, .options.snow = opts) %dopar% {
            start_state <- sample(c("m1", "m2"), size = 1, prob = c(p_1, p_2))
            datBS <- rMM(tgrid, lam1, lam2, sig1, sig2, s0 = start_state, dim = dim)
            fit <- fitMM(datBS, start = est_theta, integrControl = integrControl)$estimate
            fit
        }
    }
    
    
    if (detailBS) {
        return(list(cov = cov(na.omit(result)),
                    BS_detail = result))
    } else {
        return(cov(na.omit(result)))
    }
        
}
