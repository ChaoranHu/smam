## simulation of breaking time points
## for a 3 states telegraph process
## the start state is move
sim1.times.bbz.mov <- function(s, lam1, lam2, p) {
    tsum <- 0
    tt <- NULL
    while (TRUE) {
        tnew <- stats::rexp(1, rate=lam1)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        if (tsum > s) break
        ind <- stats::rbinom(1, size=1, prob=p) ## the prob of choosing lam2[1] is p
        tnew <- stats::rexp(1, rate=lam2[1])*ind + stats::rexp(1, rate=lam2[2])*(1-ind)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        if (tsum > s) break
    }
    tt
}

## simulation of breaking time points
## for a 3 states telegraph process
## the start state is stay or rest
sim1.times.bbz.sta <- function(s, lam1, lam2, p) {
    tsum <- 0
    tt <- NULL
    while (TRUE) {
        ind <- stats::rbinom(1, size=1, prob=p) ## the prob of choosing lam2 is p
        tnew <- stats::rexp(1, rate=lam2[1])*ind + stats::rexp(1, rate=lam2[2])*(1-ind)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        if (tsum > s) break
        tnew <- stats::rexp(1, rate=lam1)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        if (tsum > s) break
    }
    tt
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

### simulation a moving-resting-staying path given a grid time #######################

##' Sampling from a Moving-Resting-Staying Process with Embedded brownian Motion
##'
##' A moving-resting-staying process consists of three states: moving, resting and staying.
##' The transition between the three states is modeled by an alternating
##' renewal process, with expenentially distributed duration.
##' An animal stays at the same location while resting and staying
##' (the choice of resting and staying depends on bernoulli distribution),
##' and moves according to a Brownian motion while moving.
##' The sequence of states is moving, resting or staying, moving, resting or staying ...
##' or versus
##'
##' @param time time points at which observations are to be simulated
##' @param lamM rate parameter of the exponential duration while moving
##' @param lamR rate parameter of the exponential duration while resting
##' @param lamS rate parameter of the exponential duration while staying
##' @param sigma volatility parameter of the Brownian motion while moving
##' @param p probability of choosing resting,
##' and 1-p is probability of choosing staying
##' @param s0 the state at time 0, must be one of "m" or "r", for moving and
##' resting or staying, respectively
##' @param dim (integer) dimension of the Brownian motion
##'
##' @return
##' A \code{data.frame} whose first column is the time points and whose
##' other columns are coordinates of the locations.
##' @references
##' Jun Yan and Vladimir Pozdnyakov (2016). smam: Statistical Modeling of Animal
##' Movements. R package version 0.3-0. https://CRAN.R-project.org/package=smam
##' @author Chaoran Hu
##'
##' @examples
##' tgrid <- seq(0, 8000, length.out=1001)
##' dat <- rMovResHun(time=tgrid, lamM=4, lamR=0.04, lamS=0.2,
##'                   sigma=1000, p=0.5, s0="m", dim=2)
##' plot(dat$time, dat$X1, type='l')
##' plot(dat$time, dat$X2, type='l')
##' plot(dat$X1,   dat$X2, type='l')
##' 
##' @export
rMovResHun <- function(time, lamM, lamR, lamS, sigma, p, s0, dim = 2) {
    stopifnot(s0 %in% c("m", "r"))
    t0moving <- as.integer(s0 == "m")
    lam1 <- lamM
    lam2 <- c(lamR, lamS)
    tmax <- time[length(time)]
    brtimes <- sim1.times.bbz.mov(tmax, lam1, lam2, p)*t0moving
             + sim1.times.bbz.sta(tmax, lam1, lam2, p)*(1 - t0moving)
    coord <- replicate(dim, sim1.bbz(tmax, sigma, time, brtimes, t0moving))
    data.frame(time = time, coord)
}
