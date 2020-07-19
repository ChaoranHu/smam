#' @importFrom stats dnorm integrate optim rexp rnorm cov na.omit
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
#' dat3 <- rMRME(tgrid, 1, 1, 1, 0.01, "m", state = TRUE)
#' head(dat3)
#' plot(dat3[,1], dat3[,3], xlab="t", ylab="Z(t)=X(t)+GWN(0.01)", type="l")
#' 
#' @export

rMR <- function(time, lamM, lamR, sigma, s0, dim = 2, state = FALSE) {
    time <- time - time[1]
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


#' 'rMRME' samples from moving-resting process with Guassian measurement error
#' @param sig_err s.d. of Gaussian white noise
#' @rdname rMR
#' @export
rMRME <- function(time, lamM, lamR, sigma, sig_err, s0, dim = 2, state = FALSE){
    time <- time - time[1]
    dat <- rMR(time, lamM, lamR, sigma, s0, dim = dim, state)
    if (state) {
        for (i in 1:dim) {
            dat[, i+2] <- dat[, i+2] + rnorm(length(time), mean = 0, sd = sig_err)
        }
    } else {
        for (i in 1:dim) {
            dat[, i+1] <- dat[, i+1] + rnorm(length(time), mean = 0, sd = sig_err)
        }
    }
    dat
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
#' \dontrun{
#' ## time consuming example
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


#' Fit a Moving-Resting Model with Measurement Error
#'
#' 'fitMRME' fits a Moving-Resting Model with Measurement Error. The measurement
#' error is modeled by Guassian noise. Using \code{segment} to fit part
#' of observations to the model. A practical application of this feature
#' is seasonal analysis.
#'
#' @param data a data.frame whose first column is the observation time, and other
#'     columns are location coordinates. If \code{segment} is not \code{NULL},
#'     additional column with the same name given by \code{segment} should be
#'     included. This additional column is used to indicate which part of
#'     observations shoule be used to fit model. The value of this column can
#'     be any integer with 0 means discarding this observation and non-0 means
#'     using this obversvation. Using different non-zero numbers indicate different
#'     segments. (See vignette for more details.)
#' @param start starting value of the model, a vector of four components
#'     in the order of rate for moving, rate for resting, volatility, and
#'     s.d. of Guassian measurement error.
#' @param segment character variable, name of the column which indicates segments,
#'     in the given \code{data.frame}. The default value, \code{NULL}, means using
#'     whole dataset to fit the model.
#' @param lower,upper Lower and upper bound for optimization.
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
#' 
#' @references
#' Hu, C., Pozdnyakov, V., and Yan, J. Moving-resting model with measurement
#' error. In process.
#'
#' @author Chaoran Hu
#' 
#' @examples
#' ## time consuming example
#' #tgrid <- seq(0, 10*100, length=100)
#' #set.seed(123)
#' #dat <- rMRME(tgrid, 1, 0.5, 1, 0.01, "m")
#'
#' ## fit whole dataset to the MRME model
#' #fit <- fitMRME(dat, start=c(1, 0.5, 1, 0.01))
#' #fit
#'
#' ## fit whole dataset to the MRME model with naive composite likelihood
#' #fit.naive <- fitMRME_naive(dat, start=c(1, 0.5, 1, 0.01))
#' #fit.naive
#'
#' ## fit whole dataset to the MRME model with approximate error
#' #fit.approx <- fitMRMEapprox(dat, start=c(1, 0.5, 1, 0.01))
#' #fit.approx
#'
#' ## fit part of dataset to the MR model
#' #batch <- c(rep(0, 5), rep(1, 17), rep(0, 4), rep(2, 30), rep(0, 4), rep(3, 40))
#' #dat.segment <- cbind(dat, batch)
#' #fit.segment <- fitMRME(dat.segment, start = c(1, 0.5, 1, 0.01), segment = "batch")
#' #fit.segment.approx <- fitMRMEapprox(dat.segment, start = c(1, 0.5, 1, 0.01), segment = "batch")
#' #head(dat.segment)
#' #fit.segment
#' 
#' @export
fitMRME <- function(data, start, segment = NULL,
                    lower = c(0.000001, 0.000001, 0.000001, 0.000001),
                    upper = c(10, 10, 10, 10),
                    #method = "Nelder-Mead",
                    #optim.control = list(),
                    integrControl = integr.control()) {
    if (is.null(segment)) {

        
        if (!is.matrix(data)) data <- as.matrix(data)
        dinc <- apply(data, 2, diff)
        integrControl <- unlist(integrControl)

        fit <- nloptr::nloptr(x0 = start, eval_f = nllk_mrme,
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

        
        ## fit <- optim(start, nllk_mrme, data = dinc, method = method,
        ##              control = optim.control, integrControl = integrControl)

        ## estimate <- fit$par
        
        ## return(list(estimate    = estimate,
        ##             loglik      = -fit$value,
        ##             convergence = fit$convergence))

        
    } else {

        
        ## seasonal process
        result <- fitMRME_seasonal(data, segment, start,
                                   lower, upper,
                                   #method, optim.control,
                                   integrControl)
        return(result)

        
    }
}

#' Variance matrix of estimators from moving-resting process with measurement error
#'
#' 'estVarMRME_Godambe' uses Godambe information matrix to obtain variance matrix
#' of estimators from 'fitMRME'.
#' 'estVarMRME_pBootstrap' uses parametric bootstrap to obtain variance matrix
#' of estimators from 'fitMRME'.
#' 'estVarMRMEnaive_Godambe' use Godambe information matrix to obtain variance matrix
#' of estimators from 'fitMRME_naive'.
#' 'estVarMRMEnaive_pBootstrap' uses parametric bootstrap to obtain variance matrix
#' of estimators from 'fitMRME_naive'.
#'
#' @param est_theta estimators of MRME model
#' @param data data used to process estimation
#' @param nBS number of bootstrap.
#' @param numThreads the number of threads for parallel computation. If its value
#' is greater than 1, then parallel computation will be processed. Otherwise,
#' serial computation will be processed.
#' @param integrControl a list of control parameters for the \code{integrate}
#' function: rel.tol, abs.tol, subdivision.
#' @param gradMethod method used for numeric gradient (\code{numDeriv::grad}).
#'
#' @return variance-covariance matrix of estimators
#'
#' @author Chaoran Hu
#'
#' @examples
#' \dontrun{
#' ## time consuming example
#' tgrid <- seq(0, 10*100, length=100)
#' set.seed(123)
#' dat <- rMRME(tgrid, 1, 0.5, 1, 0.01, "m")
#'
#' estVarMRME_Godambe(c(1, 0.5, 1, 0.01), dat, nBS = 10)
#' estVarMRME_pBootstrap(c(1, 0.5, 1, 0.01), dat, nBS = 10)
#' estVarMRMEnaive_Godambe(c(1, 0.5, 1, 0.01), dat, nBS = 10)
#' estVarMRMEnaive_pBootstrap(c(1, 0.5, 1, 0.01), dat, nBS = 10)
#'
#' estVarMRME_Godambe(c(1, 0.5, 1, 0.01), dat, nBS = 10, numThreads = 6)
#' estVarMRME_pBootstrap(c(1, 0.5, 1, 0.01), dat, nBS = 10, numThreads = 6)
#' estVarMRMEnaive_Godambe(c(1, 0.5, 1, 0.01), dat, nBS = 10, numThreads = 6)
#' estVarMRMEnaive_pBootstrap(c(1, 0.5, 1, 0.01), dat, nBS = 10, numThreads = 6)
#' estVarMRMEnaive_pBootstrap(c(1, 0.5, 1, 0.01), dat, nBS = 10, numThreads = 6)
#' }
#' @export
estVarMRME_Godambe <- function(est_theta, data, nBS,
                               numThreads = 1,
                               gradMethod = "simple",
                               integrControl = integr.control()) {
    
    ## get J matrix in Godambe information matrix via bootstrap
    getJ_MRME <- function(est_theta, data, nBS, numThreads, gradMethod, integrControl) {
        
        tgrid <- data[, 1]
        dim <- ncol(data) - 1
        
        lamM <- est_theta[1]
        lamR <- est_theta[2]
        sigma <- est_theta[3]
        sig_err <- est_theta[4]

        p_m <- 1/lamM/(1/lamM + 1/lamR)
        p_r <- 1 - p_m

        integrControl <- unlist(integrControl)

        if (numThreads <= 1) {
            
            result <- matrix(NA, ncol = 4, nrow = nBS)

            for (i in seq_len(nBS)) {
                start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
                datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
                datBS <- as.matrix(datBS)
                dincBS <- apply(datBS, 2, diff)
                grad_cart <- numDeriv::grad(func = nllk_mrme, x = est_theta,
                                            method = gradMethod,
                                            data = dincBS,
                                            integrControl = integrControl)
                result[i, ] <- -grad_cart
            }
            
        } else {

            ## create parallel backend
            cl = parallel::makeCluster(numThreads)
            on.exit(close(pb), add = TRUE)
            on.exit(parallel::stopCluster(cl), add = TRUE)
            #doParallel::registerDoParallel(cl)
            doSNOW::registerDoSNOW(cl)
            pb <- utils::txtProgressBar(max = nBS, style = 3)
            progress <- function(n) utils::setTxtProgressBar(pb, n)
            opts <- list(progress = progress)

            i = 1 #Dummy line for Rstudio warnings

            result <- foreach(i = 1:nBS, .combine = rbind, .options.snow = opts) %dopar% {
                start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
                datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
                datBS <- as.matrix(datBS)
                dincBS <- apply(datBS, 2, diff)
                grad_cart <- numDeriv::grad(func = nllk_mrme, x = est_theta,
                                            method = gradMethod,
                                            data = dincBS,
                                            integrControl = integrControl)
                
                -grad_cart
            }

        }

        ## test code only ##
        ## print(result)
        ## test code only ##
        
        cov(na.omit(result))
    }

    ## get H matrix in Godambe information matrix
    getH_MRME <- function(est_theta, data, integrControl) {
        data <- as.matrix(data)
        dinc <- apply(data, 2, diff)
        integrControl <- unlist(integrControl)

        numDeriv::hessian(func = nllk_mrme, x = est_theta, data = dinc,
                          integrControl = integrControl)
    }

    Jmatrix <- getJ_MRME(est_theta, data, nBS, numThreads, gradMethod, integrControl)
    Hmatrix <- getH_MRME(est_theta, data, integrControl)

    solve(Hmatrix %*% solve(Jmatrix) %*% Hmatrix)

}


#' @param detailBS whether or not output estimation results during bootstrap,
#' which can be used to generate bootstrap CI.
#' @rdname estVarMRME_Godambe
#' @export
estVarMRME_pBootstrap <- function(est_theta, data, nBS, detailBS = FALSE,
                                  numThreads = 1,
                                  integrControl = integr.control()) {

    tgrid <- data[, 1]
    dim <- ncol(data) - 1
    
    lamM <- est_theta[1]
    lamR <- est_theta[2]
    sigma <- est_theta[3]
    sig_err <- est_theta[4]

    p_m <- 1/lamM/(1/lamM + 1/lamR)
    p_r <- 1 - p_m

    if (numThreads <= 1) {
        
        result <- matrix(NA, ncol = 4, nrow = nBS)

        for (i in seq_len(nBS)) {
            start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
            datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
            result[i, ] <- fitMRME(datBS, start = est_theta,
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
            start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
            datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
            fit <- fitMRME(datBS, start = est_theta, integrControl = integrControl)$estimate
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



#' 'fitMRME_naive' fits moving-resting model with measurement error
#' by MLE with a naive composite llk, that pretend two consecutive
#' increments are independent.
#' @rdname fitMRME
#' @export
fitMRME_naive <- function(data, start, segment = NULL,
                          lower = c(0.000001, 0.000001, 0.000001, 0.000001),
                          upper = c(10, 10, 10, 10),
                          #method = "Nelder-Mead",
                          #optim.control = list(),
                          integrControl = integr.control()) {
    if (is.null(segment)) {

        
        if (!is.matrix(data)) data <- as.matrix(data)
        dinc <- apply(data, 2, diff)
        integrControl <- unlist(integrControl)

        fit <- nloptr::nloptr(x0 = start, eval_f = nllk_mrme_naive_cmp,
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
        
        ## fit <- optim(start, nllk_mrme_naive_cmp, data = dinc, method = method,
        ##              control = optim.control, integrControl = integrControl)

        ## estimate <- fit$par
        
        ## return(list(estimate    = estimate,
        ##             loglik      = -fit$value,
        ##             convergence = fit$convergence))

        
    } else {

        
        ## seasonal process
        result <- fitMRME_naive_seasonal(data, segment, start,
                                         lower, upper,
                                         #method, optim.control,
                                         integrControl)
        return(result)

        
    }
}



#' @rdname estVarMRME_Godambe
#' @export
estVarMRMEnaive_Godambe <- function(est_theta, data, nBS,
                                    numThreads = 1,
                                    gradMethod = "simple",
                                    integrControl = integr.control()) {
    
    ## get J matrix in Godambe information matrix via bootstrap
    getJ_MRMEnaive <- function(est_theta, data, nBS, numThreads, gradMethod, integrControl) {
        
        tgrid <- data[, 1]
        dim <- ncol(data) - 1
        
        lamM <- est_theta[1]
        lamR <- est_theta[2]
        sigma <- est_theta[3]
        sig_err <- est_theta[4]

        p_m <- 1/lamM/(1/lamM + 1/lamR)
        p_r <- 1 - p_m

        integrControl <- unlist(integrControl)

        if (numThreads <= 1) {

            result <- matrix(NA, ncol = 4, nrow = nBS)

            for (i in seq_len(nBS)) {
                start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
                datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
                datBS <- as.matrix(datBS)
                dincBS <- apply(datBS, 2, diff)
                grad_cart <- numDeriv::grad(func = nllk_mrme_naive_cmp, x = est_theta,
                                            method = gradMethod,
                                            data = dincBS,
                                            integrControl = integrControl)
                result[i, ] <- -grad_cart
            }

        } else {

            ## create parallel backend
            cl = parallel::makeCluster(numThreads)
            on.exit(close(pb), add = TRUE)
            on.exit(parallel::stopCluster(cl), add = TRUE)
                                        #doParallel::registerDoParallel(cl)
            doSNOW::registerDoSNOW(cl)
            pb <- utils::txtProgressBar(max = nBS, style = 3)
            progress <- function(n) utils::setTxtProgressBar(pb, n)
            opts <- list(progress = progress)

            i = 1 #Dummy line for Rstudio warnings

            result <- foreach(i = 1:nBS, .combine = rbind, .options.snow = opts) %dopar% {
                start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
                datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
                datBS <- as.matrix(datBS)
                dincBS <- apply(datBS, 2, diff)
                grad_cart <- numDeriv::grad(func = nllk_mrme_naive_cmp, x = est_theta,
                                            method = gradMethod,
                                            data = dincBS,
                                            integrControl = integrControl)
                -grad_cart
            }


            
        }
        
        cov(na.omit(result))
    }

    ## get H matrix in Godambe information matrix
    getH_MRMEnaive <- function(est_theta, data, integrControl) {
        data <- as.matrix(data)
        dinc <- apply(data, 2, diff)
        integrControl <- unlist(integrControl)

        numDeriv::hessian(func = nllk_mrme_naive_cmp, x = est_theta, data = dinc,
                          integrControl = integrControl)
    }

    Jmatrix <- getJ_MRMEnaive(est_theta, data, nBS, numThreads, gradMethod, integrControl)
    Hmatrix <- getH_MRMEnaive(est_theta, data, integrControl)

    solve(Hmatrix %*% solve(Jmatrix) %*% Hmatrix)

}



#' @rdname estVarMRME_Godambe
#' @export
estVarMRMEnaive_pBootstrap <- function(est_theta, data, nBS, detailBS = FALSE,
                                       numThreads = 1,
                                       integrControl = integr.control()) {

    tgrid <- data[, 1]
    dim <- ncol(data) - 1
    
    lamM <- est_theta[1]
    lamR <- est_theta[2]
    sigma <- est_theta[3]
    sig_err <- est_theta[4]

    p_m <- 1/lamM/(1/lamM + 1/lamR)
    p_r <- 1 - p_m


    if (numThreads <= 1) {
        result <- matrix(NA, ncol = 4, nrow = nBS)

        for (i in seq_len(nBS)) {
            start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
            datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
            result[i, ] <- fitMRME_naive(datBS, start = est_theta,
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
            start_state <- sample(c("m", "r"), size = 1, prob = c(p_m, p_r))
            datBS <- rMRME(tgrid, lamM, lamR, sigma, sig_err, s0 = start_state, dim = dim)
            fit <- fitMRME_naive(datBS, start = est_theta,
                                 integrControl = integrControl)$estimate
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





#' 'fitMRMEapprox' also fits moving-resting model. However, in this function,
#' the gaussian error is approximated with two discrete distributions.
#'
#' @param approx_norm_even,approx_norm_odd numeric matrixes specify the
#'     discrete distributions used to approximate standard normal distribution.
#'     The first column is support of discrete distribution and the second
#'     column is probability mass function. \code{approx_norm_even} is used to
#'     approximate even step error and \code{approx_norm_odd} is used to
#'     approximate odd step error. We mention that the supports of these two
#'     discrete distributions should not have any common elements.
#' @rdname fitMRME
#' @export
fitMRMEapprox <- function(data, start, segment = NULL,
                          approx_norm_even = approxNormalOrder(5),
                          approx_norm_odd  = approxNormalOrder(6),
                          method = "Nelder-Mead",
                          optim.control = list(),
                          integrControl = integr.control()) {
    if (is.null(segment)) {

        
        if (!is.matrix(data)) data <- as.matrix(data)
        dinc <- apply(data, 2, diff)
        integrControl <- unlist(integrControl)
        
        fit <- optim(start, nllk_mrme_approx, data = dinc, method = method,
                     approx_norm_odd = approx_norm_odd,
                     approx_norm_even = approx_norm_even,
                     control = optim.control, integrControl = integrControl)
        
        return(list(estimate    = fit$par ,
                    loglik      = -fit$value,
                    convergence = fit$convergence))

        
    } else {

        
        ## seasonal process
        result <- fitMRMEapprox_seasonal(data, segment, start,
                                         approx_norm_odd, approx_norm_even,
                                         method, optim.control, integrControl)
        return(result)

        
    }
}






##############################################################
## the following code is for testing purpose only ############
nllk_mrme_approx_fixed_sig_err <- function(theta, sig_err, data, integrControl,
                                           approx_norm_even, approx_norm_odd) {
    nllk_mrme_approx(c(theta, sig_err), data, integrControl,
                     approx_norm_even, approx_norm_odd)
}

fitMRMEapprox_fixedSigErr <- function(data, start, sig_err,
                                      approx_norm_even = approxNormalOrder(5),
                                      approx_norm_odd  = approxNormalOrder(6),
                                      method = "Nelder-Mead",
                                      optim.control = list(),
                                      integrControl = integr.control()) {
    
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    
    fit <- optim(start, nllk_mrme_approx_fixed_sig_err,
                 data = dinc, sig_err = sig_err, method = method,
                 approx_norm_odd = approx_norm_odd,
                 approx_norm_even = approx_norm_even,
                 control = optim.control, integrControl = integrControl)
    
    return(list(estimate    = fit$par ,
                loglik      = -fit$value,
                convergence = fit$convergence))
    
}

## test code ends here #####################################
############################################################

##############################################################
## the following code is for testing purpose only ############
## fitMRME_fixed_sig_err <- function(data, start, sig_err,
##                                   method = "Nelder-Mead",
##                                   optim.control = list(),
##                                   integrControl = integr.control()){
##     ## start here contains lam1, lam0, sigma
##     ## sig_err should be given as known
##     if (!is.matrix(data)) data <- as.matrix(data)
##     dinc <- apply(data, 2, diff)
##     integrControl <- unlist(integrControl)
    
##     fit <- optim(start, nllk_mrme_fixed_sig_err, sig_err = sig_err,
##                  data = dinc, method = method,
##                  control = optim.control, integrControl = integrControl)

##     estimate <- fit$par
    
##     return(list(estimate    = estimate,
##                 loglik      = -fit$value,
##                 convergence = fit$convergence))
## }


## fitMRME_one_chain <- function(data, start,
##                               method = "Nelder-Mead",
##                               optim.control = list(),
##                               integrControl = integr.control()){
##     ## start here contains lam1, lam0, sigma, sig_err
##     if (!is.matrix(data)) data <- as.matrix(data)
##     dinc <- apply(data, 2, diff)
##     integrControl <- unlist(integrControl)
    
##     fit <- optim(start, nllk_mrme_one_chain,
##                  data = dinc, method = method,
##                  control = optim.control, integrControl = integrControl)

##     estimate <- fit$par
    
##     return(list(estimate    = estimate,
##                 loglik      = -fit$value,
##                 convergence = fit$convergence))
## }

## fitMRME_one_chain_fixed_sig_err <- function(data, start, sig_err,
##                                             method = "Nelder-Mead",
##                                             optim.control = list(),
##                                             integrControl = integr.control()){
##     ## start here contains lam1, lam0,m sigma
##     ## sig_err should be given as known
##     if (!is.matrix(data)) data <- as.matrix(data)
##     dinc <- apply(data, 2, diff)
##     integrControl <- unlist(integrControl)
    
##     fit <- optim(start, nllk_mrme_one_chain_fixed_sig_err,
##                  sig_err = sig_err,
##                  data = dinc, method = method,
##                  control = optim.control, integrControl = integrControl)

##     estimate <- fit$par
    
##     return(list(estimate    = estimate,
##                 loglik      = -fit$value,
##                 convergence = fit$convergence))
## }
## test code ends here #####################################
############################################################



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




#' Auxiliary for Preparing Discrete Distribution
#' used to approximating Standard Normal Distribution
#'
#' Auxiliary for preparing discrete distribution used to
#' approximate standard normal. This function generates
#' order statistics of standard normal with same probability
#' assigned. Then, the discrete distribution is standardized
#' to variance one and mean zero.
#'
#' @param m int, the number of order statistics used
#'
#' @details
#' This function use \code{EnvStats::evNormOrdStats} to get
#' the order statisics of standard normal distribution. The
#' same probability is assigned for each order statistics.
#'
#'
#' @return
#' A numeric matrix with first column is support of discrete
#' distribution and second column is corresponding p.m.f..
#'
#' @seealso \code{EnvStats::evNormOrdStats} for order
#' statisics of standard normal. \code{\link{fitMRMEapprox}}
#' for fit MRME with approximated measurement error.
#'
#' @author Chaoran Hu
#'
#' @export
approxNormalOrder <- function(m){
    result <- matrix(1/m, ncol = 2, nrow = m)
    result[, 1] <- EnvStats::evNormOrdStats(m)
    var <- sum((result[,1]-mean(result[,1]))^2)/m
    result[, 1] <- result[, 1] / sqrt(var)
    result
}


#' 'approxNormalOrder2' generates a even space grid first.
#' Then, the probability is calculated as normal kernal.
#' Finally, it is also standardized to variance one and
#' mean zero.
#'
#' @param width the width between two consecutive grid points.
#'
#' @rdname approxNormalOrder
#' @export
approxNormalOrder2 <- function(m, width) {
    ## generate grid with even space
    if (m %% 2 == 0) {
        grid <- seq(width/2, (m/2-1)*width + width/2, by = width)
        grid <- c(-rev(grid), grid)
    } else {
        grid <- seq(0, ((m-1)/2)*width, by = width)
        grid <- c(-rev(grid), grid[-1])
    }
    ## generate probability with norm kernel
    result <- cbind(grid, dnorm(grid)/sum(dnorm(grid)))
    ## standardize it to var 1 and mean 0
    var <- sum(((result[,1]-mean(result[,1]))^2) * result[, 2])
    result[, 1] <- result[, 1] / sqrt(var)
    result
}







##############################################################
## the following code is for testing purpose only ############
## The R version of composite likelihood estimation
## Kept only for comparison check with fitMR using cl = TRUE
## fitMovRes.cl <- function(data, start, logtr = FALSE, method = "Nelder-Mead",
##                          optim.control = list()) {
##     if (!is.matrix(data)) data <- as.matrix(data)
##     dinc <- apply(data, 2, diff)
##     fit <- optim(start, ncllk.m1.inc, data = dinc, method = method,
##                  control = optim.control, logtr = logtr)
##     list(estimate    = fit$par,
##          loglik      = -fit$value,
##          convergence = fit$convergence)
## }


## ### not public function ###
## ## llk for moving-resting model with given state
## ## param data: time state locations
## llk_mr <- function(data, theta, integrControl = integr.control()) {
##     if (!all(data[, 2] == 0 | data[, 2] == 1)) stop("state must be 0 or 1")
##     if (!is.matrix(data)) data <- as.matrix(data)
##     dinc <- apply(data[, -2], 2, diff)
##     integrControl <- unlist(integrControl)
##     mrllk_state(theta, dinc, data[, 2], integrControl)
## }
## test code ends here #####################################
############################################################






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
