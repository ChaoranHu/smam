#' @importFrom stats dnorm integrate optim rexp rnorm
#' @importFrom methods is
#' @importFrom numDeriv hessian
#' @importFrom Rcpp evalCpp
#' @useDynLib smam

## simulation of breaking time points bbs
## for a 2 state telegraph process
sim1mr.times.bbz <- function(s, lam1, lam2) {
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
sim1mr.bbz <- function(s, sigma, time, brtimes, t0moving=TRUE) {
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
            if (status == 1)   ## moving
              x[i] <- x[i - 1] + rnorm(1, sd = sigma * sqrt(tt[i] - tt[i - 1]))
            else               ## resting
              x[i] <- x[i - 1]
        }
        if (tt[i] == tend) {
            status <- 1 - status ## switch status
            j <- j + 1
            tend <- brtimes[j]
        }
    }
    x[tt %in% time]
}

## figure out the state given breaking times
simmr.state <- function(time, brtimes, t0moving) {
    brstate <- brtimes[, 2]
    brtimes <- brtimes[, 1]
    tt <- sort(unique(c(time, brtimes)))
    nt <- length(tt)
    nb <- length(brtimes) 
    x <- rep(NA, nt)
    tend <- brtimes[1]
    tend.state <- brstate[1]
    j <- 1
    for (i in 1:nt) {
        if (tt[i] <= tend) { ## status unchanged
            x[i] <- tend.state
        }
        if (tt[i] == tend) { ## switch status
            j <- j + 1
            tend <- brtimes[j]
            tend.state <- brstate[j]
        }
    }
    
    ## return result according to the start state
    if (t0moving) {
        return(1 - x[tt %in% time])
    } else {
        return(x[tt %in% time])
    }
}


## simulation a moving-resting path given a grid time

#' Sampling from a Moving-Resting Process with Embedded Brownian Motion
#'
#' A moving-resting process consists of two states: moving and resting.
#' The transition between the two states is modeled by an alternating
#' renewal process, with exponentially distributed duration. An animal
#' stays at the same location while resting, and moves according to a
#' Brownian motion while moving.
#'
#' @param time time points at which observations are to be simulated
#' @param lamM rate parameter of the exponential duration while moving
#' @param lamR rate parameter of the exponential duration while resting
#' @param sigma volatility parameter of the Brownian motion while moving
#' @param s0 the state at time 0, must be one of "m" or "r", for moving and
#' resting, respectively
#' @param dim (integer) dimension of the Brownian motion
#' @param state indicates whether the simulation show the states at given
#' time points.
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
#' tgrid <- seq(0, 10, length=1001)
#' ## make it irregularly spaced
#' tgrid <- sort(sample(tgrid, 800))
#' dat <- rMR(tgrid, 1, 1, 1, "m")
#' plot(dat[,1], dat[,2], xlab="t", ylab="X(t)", type='l')
#'
#' dat2 <- rMR(tgrid, 1, 1, 1, "m", state = TRUE)
#' head(dat2)
#' 
#' @export

rMR <- function(time, lamM, lamR, sigma, s0, dim = 2, state = FALSE) {
    stopifnot(s0 %in% c("m", "r"))
    t0moving <- (s0 == "m")
    lam1 <- if (t0moving) lamM else lamR
    lam2 <- if (t0moving) lamR else lamM
    
    timeIND <- 0
    if (length(time) == 1) {
        timeIND <- 1
        time <- c(0, time)
    }
    tmax <- time[length(time)]
    brtimes <- sim1mr.times.bbz(tmax, lam1, lam2)
    coord <- replicate(dim, sim1mr.bbz(tmax, sigma, time, brtimes[, 1], t0moving))

    stateresult <- simmr.state(time, brtimes, t0moving = t0moving)

    if (timeIND == 1) {
        if (state) {
            return(data.frame(time = time, state = stateresult, coord)[-1, ])
        }
        return(data.frame(time = time, coord)[-1, ])
    }
    if (state) {
        return(data.frame(time = time, state = stateresult, coord))
    }
    data.frame(time = time, coord)
}

#' 'rMovRes' is deprecated. Using new function 'rMR' instead.
#' @rdname rMR
#' @export
rMovRes <- function(time, lamM, lamR, sigma, s0, dim = 2) {
    .Deprecated("rMR")
    rMR(time, lamM, lamR, sigma, s0, dim)
}


## atom at w = t for dtm.m or dtm.r
dt.atm <- function(t, lam) {
    exp(- lam * t)
}

## dt{mr}.{mr}{mr}
## prob/dens of time spent in state (m or r) by time t,
## the starting state is m or r and ending state is m or r.
##
## dt{mr}.{mr}
## prob/dens of time spent in state (m or r) by time t,
## the starting state is m or r.
##
## using besselI function is much faster than my C code
dtm.mm.b <- function(w, t, lamM, lamR) {
    ## The formula of Zacks (2004, p.500, JAP:497-507) mistakenly has
    ## 2 inside the square root!
    lm <- lamM * w
    lr <- lamR * (t - w)
    elmr <- exp(-lm - lr)
    z <- 2 * sqrt(lm * lr)
    ## ifelse(w > t | w < 0, 0, 
    ##        ifelse(w == t, atm(t, lamM),
    ##               elmr * sqrt(lamM * lamR * w / (t - w)) * besselI(z, 1)))
    ifelse(w > t | w < 0, 0,  elmr * sqrt(lamM * lamR * w / (t - w)) * besselI(z, 1))
}

dtm.mr.b <- function(w, t, lamM, lamR) {
    lm <- lamM * w
    lr <- lamR * (t - w)
    elmr <- exp(-lm - lr)
    z <- 2 * sqrt(lm * lr)
    ifelse(w > t | w < 0, 0, lamM * elmr * besselI(z, 0))
}

## calling C code yields to using modified besselI (mb) as default
dtm.mm <- function(w, t, lamM, lamR, mb = TRUE) { 
    if (mb) return(dtm.mm.b(w, t, lamM, lamR))
    wlen <- length(w)
    t <- rep(t, length.out = wlen)
    dens <- .C("pmm", as.double(w), as.double(t), as.double(lamM), as.double(lamR), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
    ##           ifelse(w == t, atm(t, lamM), dens))
}

dtm.mr <- function(w, t, lamM, lamR, mb = TRUE) {
    if (mb) return(dtm.mr.b(w, t, lamM, lamR))
    wlen <- length(w)
    t <- rep(t, length.out = wlen)
    dens <- .C("pmr", as.double(w), as.double(t), as.double(lamM), as.double(lamR), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
}

dtr.rr <- function(w, t, lamM, lamR, mb = TRUE) {
    if (mb) return(dtm.mm.b(w, t, lamR, lamM))
    wlen <- length(w)
    t <- rep(t, length.out=wlen)
    dens <- .C("pmm", as.double(w), as.double(t), as.double(lamR), as.double(lamM), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
    ##           ifelse(w == t, atm(t, lamR), dens))
}

dtr.rm <- function(w, t, lamM, lamR, mb = TRUE) {
    if (mb) return(dtm.mr.b(w, t, lamR, lamM))
    wlen <- length(w)
    t <- rep(t, length.out=wlen)
    dens <- .C("pmr", as.double(w), as.double(t), as.double(lamR), as.double(lamM), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
}


dtm.m <- function(w, t, lamM, lamR) {
    dtm.mm(w, t, lamM, lamR) + dtm.mr(w, t, lamM, lamR)
}

dtm.r <- function(w, t, lamM, lamR) {
    dtr.rm(t - w, t, lamM, lamR) + dtr.rr(t - w, t, lamM, lamR)
}
    
dtr.m <- function(w, t, lamM, lamR) {
    dtm.mm(t - w, t, lamM, lamR) + dtm.mr(t - w, t, lamM, lamR)
}

dtr.r <- function(w, t, lamM, lamR) {
    dtr.rm(w, t, lamM, lamR) + dtr.rr(w, t, lamM, lamR)
}

#' Density for Time Spent in Moving or Resting
#'
#' Density for time spent in moving or resting in a time interval,
#' unconditional or conditional on the initial state.
#'
#' @param w time points at which the density is to be evaluated
#' @param t length of the time interval
#' @param lamM rate parameter of the exponentially distributed duration in moving
#' @param lamR rate parameter of the exponentially distributed duration in resting
#' @param s0 initial state. If \code{NULL}, the unconditional density is
#' returned; otherwise, it is one of "m" or "s", standing for moving and
#' resting, respectively, and the conditional density is returned given
#' the initial state.
#' @return
#' a vector of the density evaluated at \code{w}.
#' @details
#' \code{dtm} returns the density for time in moving;
#' \code{dtr} returns the density for time in resting.
#' @references
#' Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakov, V.,
#' Williams, S., and Meyer, T. (2014) A moving-resting process with an
#' embedded Brownian motion for animal movements.
#' Population Ecology. 56(2): 401--415.
#'
#' @examples
#' lamM <- 1
#' lamR <- c(1/2, 1, 2)
#' lr <- length(lamR)
#' totalT <- 10
#' old.par <- par(no.readonly=TRUE)
#' par(mfrow=c(1, 2), mar=c(2.5, 2.5, 1.1, 0.1), mgp=c(1.5, 0.5, 0), las=1)
#' curve(dtm(x, totalT, 1, 1/2, "m"), 0, totalT, lty=1, ylim=c(0, 0.34),
#'       xlab="M(10)", ylab="density")
#' curve(dtm(x, totalT, 1, 1, "m"), 0, totalT, lty=2, add=TRUE)
#' curve(dtm(x, totalT, 1, 2, "m"), 0, totalT, lty=3, add=TRUE)
#' mtext(expression("S(0) = 1"))
#' legend("topleft", legend = expression(lambda[r] == 1/2, lambda[r] == 1,
#'        lambda[r] == 2), lty = 1:lr)
#' curve(dtm(x, totalT, 1, 1/2, "r"), 0, totalT, lty=1, ylim=c(0, 0.34),
#'       xlab="M(10)", ylab="density")
#' curve(dtm(x, totalT, 1, 1, "r"), 0, totalT, lty=2, add=TRUE)
#' curve(dtm(x, totalT, 1, 2, "r"), 0, totalT, lty=3, add=TRUE)
#' mtext(expression("S(0) = 0"))
#' legend("topleft", legend = expression(lambda[r] == 1/2, lambda[r] == 1,
#'       lambda[r] == 2), lty = 1:lr)
#' par(old.par)
#'
#' @export

dtm <- function(w, t, lamM, lamR, s0 = NULL) {
    if (is.null(s0)) { ## unconditional
        pm <- lamR / (lamM + lamR)
        pr <- 1 - pm
        pm * dtm.m(w, t, lamM, lamR) + pr * dtm.r(w, t, lamM, lamR)
    }
    else { ## conditional on s0
        if (s0 == "m") dtm.m(w, t, lamM, lamR)
        else dtm.r(w, t, lamM, lamR)
    }
}

#' @describeIn dtm Density of time spent in resting
#' @export
dtr <- function(w, t, lamM, lamR, s0 = NULL) {
    if (is.null(s0)) { ## unconditional
        pm <- lamR / (lamM + lamR)
        pr <- 1 - pm
        pm * dtr.m(w, t, lamM, lamR) + pr * dtr.r(w, t, lamM, lamR)
    }
    else { ## conditional on s0
        if (s0 == "m") dtr.m(w, t, lamM, lamR)
        else dtr.r(w, t, lamM, lamR)
    }
}

## density of process at time s, in either state (m or r), given that
## the starting state is m or r

hmm <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    term1 <- exp(-lamM * s) * apply(dnorm(x, sd = sigma * sqrt(s)), 1, prod)
    term2 <- sapply(1:n,
                    function(i) {
                        f <- function(w) {
                            t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(w)))
                            t2 <- dtm.mm(w, s[i], lamM, lamR)
                            t1 * t2
                        }
                        vf <- Vectorize(f)
                        integrate(vf, 0, s[i])$value
                    })
    term1 + term2
}

hmr <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    sapply(1:n,
           function(i) {
               f <- function(w) {
                   t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(w)))
                   t2 <- dtm.mr(w, s[i], lamM, lamR)
                   t1 * t2
               }
               vf <- Vectorize(f)
               integrate(vf, 0, s[i])$value
           })
}


hrr <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    sapply(1:n,
           function(i) {
               f <- function(w) {
                   t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(s[i] - w)))
                   t2 <- dtr.rr(w, s[i], lamM, lamR)
                   t1 * t2
               }
               vf <- Vectorize(f)
               integrate(vf, 0, s[i])$value
           })
}

hrm <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    sapply(1:n,
           function(i) {
               f <- function(w) {
                   t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(s[i] - w)))
                   t2 <- dtr.rm(w, s[i], lamM, lamR)
                   t1 * t2
               }
               vf <- Vectorize(f)
               integrate(vf, 0, s[i])$value
           })
}

hm <- function(x, s, lamM, lamR, sigma) {
    ## integrate(hm, -infty, +infty) = 1 
    hmm(x, s, lamM, lamR, sigma) + hmr(x, s, lamM, lamR, sigma)
}

hr <- function(x, s, lamM, lamR, sigma) {
    ## integrate(hr, -infty, +infty) = 1 - exp(-lamR * s)
    hrm(x, s, lamM, lamR, sigma) + hrr(x, s, lamM, lamR, sigma)
}

## marginal likelihood for inferences
## marginal density of one increment
## x: the increment (not the original locations)
dmrbm.m1 <- function(x, tt, lamM, lamR, sigma) {
    if (lamM <= 0 || lamR <=0 || sigma <= 0) return(NA)
    if (!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    ## delta <- apply(x, 1, prod) == 0
    delta <- apply(x == 0, 1, all) ## staying
    pm <- 1 / lamM / (1 / lamM + 1 / lamR)
    pr <- 1 - pm
    tt <- rep(tt, length.out = n)
    dens <- pr * exp(-lamR * tt)
    if (!all(delta)) {
        x <- x[!delta, , drop=FALSE]
        tt <- tt[!delta]
        hm <- hmm(x, tt, lamM, lamR, sigma) + hmr(x, tt, lamM, lamR, sigma)
        hr <- hrm(x, tt, lamM, lamR, sigma) + hrr(x, tt, lamM, lamR, sigma)
        dens[!delta] <- pm * hm + pr * hr
    }
    dens
}

## likelihood for given one set of parameter value
cllk.m1.p <- function(x, t, lamM, lamR, sigma) { 
    cllk <- try(sum(log(dmrbm.m1(x, t, lamM, lamR, sigma))), silent=TRUE)
    if (inherits(cllk, "try-error")) NA else cllk
}


cllk.m1.p2v <- Vectorize(cllk.m1.p, c("lamM", "lamR", "sigma"))

cllk.m1 <- function(theta, data) {
    lamM <- theta[1]
    lamR <- theta[2]
    sigma <- theta[3]
    dinc <- apply(data, 2, diff)
    cllk.m1.p2v(dinc[,-1], dinc[,1], lamM, lamR, sigma)    
}

ncllk.m1.inc <- function(theta, data, logtr = FALSE) { ## data is increment already
    if (logtr) theta <- exp(theta)
    lamM <- theta[1]
    lamR <- theta[2]
    sigma <- theta[3]
    -cllk.m1.p2v(data[,-1], data[,1], lamM, lamR, sigma)    
}


#' Fit a Moving-Resting Model with Embedded Brownian Motion
#'
#' Fit a Moving-Resting Model with Embedded Brownian Motion with animal
#' movement data at discretely observation times by maximizing a composite
#' likelihood constructed from the marginal density of increment.
#' Using \code{segment} to fit part of observations to the model. A practical
#' application of this feature is seasonal analysis.
#'
#' @param data a data.frame whose first column is the observation time, and other
#'     columns are location coordinates. If \code{segment} is not \code{NULL},
#'     additional column with the same name given by \code{segment} should be
#'     included. This additional column is used to indicate which part of
#'     observations shoule be used to fit model. The value of this column can
#'     be any integer with 0 means discarding this observation and non-0 means
#'     using this obversvation. Using different non-zero numbers indicate different
#'     segments. (See vignette for more details.)
#' @param start starting value of the model, a vector of three components
#'     in the order of rate for moving, rate for resting, and volatility.
#' @param segment character variable, name of the column which indicates segments,
#'     in the given \code{data.frame}. The default value, \code{NULL}, means using
#'     whole dataset to fit the model.
#' @param likelihood a character string specifying the likelihood type to
#'     maximize in estimation. This can be "full" for full likelihood or
#'     "composite' for composite likelihood.
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
#' \item{likelihood}{likelihood type (full or composite) from the input}
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
#' \donttest{
#' tgrid <- seq(0, 10, length=500)
#' set.seed(123)
#' ## make it irregularly spaced
#' tgrid <- sort(sample(tgrid, 30)) # change to 400 for a larger sample
#' dat <- rMR(tgrid, 1, 2, 25, "m")
#'
#' ## fit whole dataset to the MR model
#' fit.fl <- fitMR(dat, start=c(2, 2, 20), likelihood = "full")
#' fit.fl
#' 
#' fit.cl <- fitMR(dat, start=c(2, 2, 20), likelihood = "composite")
#' fit.cl
#'
#' ## fit part of dataset to the MR model
#' batch <- c(rep(0, 5), rep(1, 7), rep(0, 4), rep(2, 10), rep(0, 4))
#' dat.segment <- cbind(dat, batch)
#' fit.segment <- fitMR(dat.segment, start = c(2, 2, 20), segment = "batch",
#'                      likelihood = "full")
#' head(dat.segment)
#' fit.segment
#' }
#' 
#' @export
fitMR <- function(data, start, segment = NULL,
                  likelihood = c("full", "composite"),
                  logtr = FALSE,
                  method = "Nelder-Mead",
                  optim.control = list(),
                  integrControl = integr.control()) {
    if (is.null(segment)) {

        
        ## normal process
        if (!is.matrix(data)) data <- as.matrix(data)
        dinc <- apply(data, 2, diff)
        objfun <- switch(likelihood,
                         composite = ncllk_m1_inc,
                         full = nllk_inc,
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
        
        return(list(estimate    = estimate,
                    varest      = varest,
                    loglik      = -fit$value,
                    convergence = fit$convergence,
                    likelihood  = likelihood))

        
    } else {

        
        ## seasonal process
        result <- fitMR_seasonal(data, segment, start, likelihood,
                                 logtr, method, optim.control, integrControl)
        return(result)

        
    }
}

#' 'fitMovRes' is deprecated. Using new function 'fitMR' instead.
#' @rdname fitMR
#' @export
fitMovRes <- function(data, start, likelihood = c("full", "composite"),
                      logtr = FALSE,
                      method = "Nelder-Mead",
                      optim.control = list(),
                      integrControl = integr.control()) {
    .Deprecated("fitMR")
    fitMR(data= data, start = start, likelihood = likelihood,
          logtr = logtr, method = method, optim.control = optim.control,
          integrControl = integrControl)
}

#' Auxiliary for Controlling Numerical Integration
#'
#' Auxiliary function for the numerical integration used in the
#' likelihood and composite likelihood functions. Typically only
#' used internally by 'fitMR' and 'fitMRH'. 
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
    


## The R version of composite likelihood estimation
## Kept only for comparison check with fitMR using cl = TRUE
fitMovRes.cl <- function(data, start, logtr = FALSE, method = "Nelder-Mead",
                         optim.control = list()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    fit <- optim(start, ncllk.m1.inc, data = dinc, method = method,
                 control = optim.control, logtr = logtr)
    list(estimate    = fit$par,
         loglik      = -fit$value,
         convergence = fit$convergence)
}


## dx <- function(x, t, lamM, lamR, Sigma) {
##     pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##     pr <- 1 - pm
##     hm <- dxmm(x, t, lamM, lamR, Sigma) + dxmr(x, t, lamM, lamR, Sigma)
##     hr <- dxrm(x, t, lamM, lamR, Sigma) + dxrr(x, t, lamM, lamR, Sigma)
##     pm * hm + pr * hr
## }

## cllk.1inc <- function(theta, xinc, t) {
##     Sigma <- diag(rep(theta[3], 2))
##     sum(log(dx(xinc, t, theta[1], theta[2], Sigma)))
## }






## #################################################################
## ## old code
## #################################################################
## dmnorm <- function(x, Sigma, w) {
##   ## x is a single observation vector
##   ## w is a vector
##   ## need to return a vector of the same length as w
##   if (!is.matrix(x)) x <- as.matrix(x)
##   if (!is.matrix(Sigma)) Sigma <- as.matrix(Sigma)
##   sapply(w, function(t) dmvnorm(x, sigma = t * Sigma))
## }
   
## dxmm <- function(x, s, lamM, lamR, Sigma) {
##   if (is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, w) * dtm.mm(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##     val[i] <- val[i] + exp(-lamM * s[i]) * dmnorm(matrix(x[i, ], 1, d), Sigma, s[i])
##   }
##   val
## }

## ## hmr(x, t)
## dxmr <- function(x, s, lamM, lamR, Sigma) {
##   if(is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, w) * dtm.mr(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##   }
    
##   val
## }

## ## hrr(x, t)
## dxrr <- function(x, s, lamM, lamR, Sigma) {
##   if(is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, s[i] - w) * dtr.rr(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##   }

##   val
## }

## ## hrm(x, t)
## dxrm <- function(x, s, lamM, lamR, Sigma) {
##   if(is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, s[i] - w) * dtr.rm(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##   }

##   val
## }


## ## conditional density
## ## Pr[ X_s in (x, x + dx), X_t - X_s in (y, y + dy) | s(0) = 1]
## dxxm <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##   pr <- 1 - pm
##   dxxm.m <- dxmm(x, s, lamM, lamR, Sigma) * (dxmm(y, t - s, lamM, lamR, Sigma) +
##                                             dxmr(y, t - s, lamM, lamR, Sigma))
##   dxxm.r <- dxmr(x, s, lamM, lamR, Sigma) * (dxrr(y, t - s, lamM, lamR, Sigma) +
##                                             dxrm(y, t - s, lamM, lamR, Sigma))
##   dxxm.m+dxxm.r 
## }

## ## Pr[ X_s in (x, x + dx), X_t - X_s in (y, y + dy) | s(0) = 0]
## dxxr <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##   pr <- 1 - pm
##   dxxr.m <- dxrm(x, s, lamM, lamR, Sigma) * (dxmm(y, t - s, lamM, lamR, Sigma) +
##                                             dxmr(y, t - s, lamM, lamR, Sigma))
##   dxxr.r <- dxrr(x, s, lamM, lamR, Sigma) * (dxrr(y, t - s, lamM, lamR, Sigma) +
##                                             dxrm(y, t - s, lamM, lamR, Sigma))
##   dxxr.m+dxxr.r  
## }

## ## Pr[ X_s in (x, x + dx), X_t - X_s in (y, y + dy) | s(0) = 1 or 0]
## dxx <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##   pr <- 1 - pm
##   dxxm(x, y, s, t, lamM, lamR, Sigma) * pm + dxxr(x, y, s, t, lamM, lamR, Sigma) * pr  
## }



## BDd <- function(y, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)
##   hm <- dxmm(y, t, lamM, lamR, Sigma) + dxmr(y, t, lamM, lamR, Sigma)
##   hr <- dxrm(y, t, lamM, lamR, Sigma) + dxrr(y, t, lamM, lamR, Sigma)
##   bdd <- pm * hm + pr * hr
##   bdd
## }

## ## pr[ X(u) in dx | X(t) = y] 
## BD <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)

##   if(is.vector(x)) x <- as.matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   zero <- ifelse(d == 1, 0, matrix(0, 1, d))

##   bd <- matrix(NA, length(s), xlen)
##   bdd <- BDd(y, t, lamM, lamR, Sigma)

##   for(j in 1:xlen) {  
    
##      # x != 0, y != 0, x != Y
##     if((!identical(matrix(as.numeric(x[j,]), 1, d), zero)) & (!identical(y, zero)) & (!identical(matrix(as.numeric(x[j,]), 1, d), y))) {
##       bdn <- sapply(1:length(s), function(i) dxx(matrix(x[j,], 1, d), y - matrix(x[j,], 1, d), s[i], t, lamM, lamR, Sigma))
##      # bdn <- sapply(1:length(s), function(i) pm * dxxm(x = matrix(x[j,], 1, d), y = y - matrix(x[j,], 1, d), s[i], t, lamM, lamR, Sigma) + pr * dxxr(x = matrix(x[j,], 1, d), y = y - x[j,], s[i], t, lamM, lamR, Sigma))
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, bdn[i] / bdd, 0))
##     }
##     # x = 0, y != 0
##     else if ((identical(matrix(as.numeric(x[j,]), 1, d), zero)) & (!identical(y, zero))) {  
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(s[i])) * (dxrm(y, t-s[i], lamM, lamR, Sigma)+dxrr(y, t-s[i], lamM, lamR, Sigma))) / bdd, 0))
##     }
##     # x = y, y != 0
##     else if ((identical(matrix(as.numeric(x[j,]), 1, d), y)) & (!identical(y, zero))) {  
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(t-s[i])) * (dxrm(y, s[i], lamM, lamR, Sigma)+dxrr(y, s[i], lamM, lamR, Sigma))) / bdd, 0))
##     }
##     # y = 0
##     else if (identical(y, zero)) {
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(identical(matrix(as.numeric(x[j,]), 1, d), zero), 1, 0))
##     }


##   }

##   bd
## }


## BDIni <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)

##   if(is.vector(x)) x <- as.matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   zero <- matrix(0, 1, d)

##   bd <- matrix(NA, length(s), xlen)
##   bdd <- BDd(y, t, lamM, lamR, Sigma)

##   for(j in 1:xlen) {
##     # x = 0, y != 0
##     bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(s[i])) * (dxrm(y, t-s[i], lamM, lamR, Sigma)+dxrr(y, t-s[i], lamM, lamR, Sigma))) / bdd, 0))

##   }

##   bd
## }

## BDEnd <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)

##   if(is.vector(x)) x <- as.matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   zero <- matrix(0, 1, d)

##   bd <- matrix(NA, length(s), xlen)
##   bdd <- BDd(y, t, lamM, lamR, Sigma)

##   for(j in 1:xlen) {
##     # x = y, y != 0
##     bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(t-s[i])) * (dxrm(y, s[i], lamM, lamR, Sigma)+dxrr(y, s[i], lamM, lamR, Sigma))) / bdd, 0))
  
##   }

##   bd
## }
