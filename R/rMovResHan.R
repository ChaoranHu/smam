## simulation of breaking time points
## for a 3 states telegraph process
## the start state is move
sim1.times.bbz.mov <- function(s, lam1, lam2, p) {
    tsum <- 0
    tt <- NULL
    state <- NULL
    while (TRUE) {
        tnew <- stats::rexp(1, rate=lam1)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        state <- c(state, 0)
        if (tsum > s) break
        ind <- stats::rbinom(1, size=1, prob=p) ## the prob of choosing lam2[1] is p
        tnew <- stats::rexp(1, rate=lam2[1])*ind + stats::rexp(1, rate=lam2[2])*(1-ind)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        state <- c(state, (2 - ind))
        if (tsum > s) break
    }
    cbind(tt, state)
}

## simulation of breaking time points
## for a 3 states telegraph process
## the start state is stay or rest
## ini.state = 1 or 2
sim1.times.bbz.sta <- function(s, lam1, lam2, p, ini.state) {
    if (ini.state == 1) {
        tsum <- rexp(1, rate = lam2[1])
        tt <- tsum
        state <- 1
        if (tsum > s) return(cbind(tt, state))
    } else {
        tsum <- rexp(1, rate = lam2[2])
        tt <- tsum
        state <- 2
        if (tsum > s) return(cbind(tt, state))
    }
    
    while (TRUE) {
        
        tnew <- stats::rexp(1, rate=lam1)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        state <- c(state, 0)
        if (tsum > s) break

        ind <- stats::rbinom(1, size=1, prob=p) ## the prob of choosing lam2 is p
        tnew <- stats::rexp(1, rate=lam2[1])*ind + stats::rexp(1, rate=lam2[2])*(1-ind)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        state <- c(state, (2 - ind))
        if (tsum > s) break
        
    }
    cbind(tt, state)
}


## simulation of a realization given breaking times
sim1.bbz <- function(s, sigma, time, brtimes, t0moving=TRUE) {
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
              x[i] <- x[i - 1] + stats::rnorm(1, sd = sigma * sqrt(tt[i] - tt[i - 1]))
            else               ## resting or staying
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
sim.state <- function(time, brtimes) {
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
    x[tt %in% time]
}

### simulation a moving-resting-handling path given a grid time #######################

##' Sampling from a Moving-Resting-Handling Process with Embedded Brownian Motion
##'
##' A moving-resting-handling process consists of three states: moving, resting and handling.
##' The transition between the three states is modeled by an alternating
##' renewal process, with expenentially distributed duration.
##' An animal stays at the same location while resting and handling
##' (the choice of resting and handling depends on Bernoulli distribution),
##' and moves according to a Brownian motion while moving state.
##' The sequence of states is moving, resting or staying, moving, resting or staying ...
##' or versus
##'
##' @param time time points at which observations are to be simulated
##' @param lamM rate parameter of the exponential duration while moving
##' @param lamR rate parameter of the exponential duration while resting
##' @param lamH rate parameter of the exponential duration while handling
##' @param sigma volatility parameter of the Brownian motion while moving
##' @param p probability of choosing resting,
##' and 1-p is probability of choosing handling
##' @param s0 the state at time 0, must be one of "m" (moving), "r" (resting) or "h" (handling).
##' @param dim (integer) dimension of the Brownian motion
##' @param state indicates whether the simulation show the states at given
##' time points.
##'
##' @return
##' A \code{data.frame} whose first column is the time points and whose
##' other columns are coordinates of the locations. If \code{state} is
##' \code{TRUE}, the second column will be the simulation state.
##' 
##' @references
##' Pozdnyakov, V., Elbroch, L.M., Hu, C., Meyer, T., and Yan, J. (2018+)
##' On estimation for Brownian motion governed by telegraph process with
##' multiple off states. <arXiv:1806.00849>
##'
##' @seealso \code{\link{fitMRH}} for fitting model.
##' 
##' @examples
##' set.seed(06269)
##' tgrid <- seq(0, 8000, length.out=1001)
##' dat <- rMRH(time=tgrid, lamM=4, lamR=0.04, lamH=0.2,
##'             sigma=1000, p=0.5, s0="m", dim=2)
##' plot(dat$time, dat$X1, type='l')
##' plot(dat$time, dat$X2, type='l')
##' plot(dat$X1,   dat$X2, type='l')
##'
##' set.seed(06269) ## show the usage of state
##' dat2 <- rMRH(time=tgrid, lamM=4, lamR=0.04, lamH=0.2,
##'              sigma=1000, p=0.5, s0="m", dim=2, state=TRUE)
##' head(dat)
##' head(dat2)
##'
##' @author Chaoran Hu
##' @export
rMRH <- function(time, lamM, lamR, lamH, sigma, p, s0, dim = 2, state = FALSE) {
    stopifnot(s0 %in% c("m", "r", "h"))
    t0moving <- as.integer(s0 == "m")
    lam1 <- lamM
    lam2 <- c(lamR, lamH)

    time <- time - time[1]

    timeIND <- 0
    if (length(time) == 1) {
        timeIND <- 1
        time <- c(0, time)
    }
    tmax <- time[length(time)]
    if (t0moving == 1) {
        brtimes <- sim1.times.bbz.mov(tmax, lam1, lam2, p)
    } else {
        if (s0 == "r") {
            brtimes <- sim1.times.bbz.sta(tmax, lam1, lam2, p, ini.state = 1)
        } else {
            brtimes <- sim1.times.bbz.sta(tmax, lam1, lam2, p, ini.state = 2)
        }
    }
    coord <- replicate(dim, sim1.bbz(tmax, sigma, time, brtimes[, 1], t0moving))

    stateresult <- sim.state(time, brtimes)

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
