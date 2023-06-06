#' Sampling from a Brownian bridge path give a grid time
#'
#' A Brownian bridge is a continuous-time stochastic process
#' B(t) whose probability distribution is the conditional
#' probability distribution of a standard Wiener process W(t)
#' subject to the condition (when standardized) that W(T) = 0,
#' so that the process is pinned to the same value at both t = 0 and t = T.
#' The implementation here is a generalized Brownian bridge
#' that allows start point and end point at different locations.
#'
#' @param time time points at which observations are to be simulated
#' @param start_pt the start point location of Brownian bridge
#' @param end_pt the end point location of Brownian brige
#' @param sigma volatility parameter of the Brownian motion
#'
#' @return
#' A \code{data.frame} whose first column is the time points
#' and second column is coordinate of the locations.
#' 
#' @examples
#' ## Brownian bridge starting from location 0 and ending at location 1
#' ## with sigma 0.1 from time 0 to time 10
#' plot(rBB(seq(0, 10, length.out = 100), 0, 1, 0.1), type = "l")
#' 
#' @export
rBB <- function(time, start_pt, end_pt, sigma) {
  time <- time - time[1]
  ntime <- length(time)
  u <- time[1]
  t <- time[ntime]
  alpha <- start_pt
  beta <- end_pt
  
  sim_BB <- numeric(ntime)
  for (i in seq_len(ntime)) {
    s <- time[i]
    cart_mean <- ((t-s)*alpha + (s-u)*beta)/(t-u)
    cart_sd <- sqrt(sigma^2 * (s-u)*(t-s)/(t-u))
    sim_BB[i] <- rnorm(1, cart_mean, cart_sd)
  }
  
  data.frame(time = time, BB = sim_BB)
}




sim1mrb.bbz <- function(s, sigma, time, brtimes, t0moving=TRUE,
                        start_pt, end_pt) {
  #### time: time points in [0, s]
  tt <- sort(unique(c(time, brtimes)))
  if (time[length(time)] < brtimes[length(brtimes)]) {
      tt <- tt[seq_len(length(tt)-1)]
  }
  nt <- length(tt)
  nb <- length(brtimes)
  status <- as.integer(t0moving)
  tend <- brtimes[1]
  j <- 1
  
  mov_segments <- numeric(0)
  for (i in 2:nt) {
    if (tt[i] <= tend) { ## status unchanged
      if (status == 1)   ## moving
        mov_segments <- c(mov_segments, tt[i] - tt[i - 1])
    }
    if (tt[i] == tend) {
      status <- 1 - status ## switch status
      j <- j + 1
      tend <- brtimes[j]
    }
  }
  
  sim_BB <- rBB(c(0, cumsum(mov_segments)), start_pt, end_pt, sigma)[, 2]
  
  x <- rep(NA, nt)
  x[1] <- start_pt
  status <- as.integer(t0moving)
  tend <- brtimes[1]
  j <- 1
  k <- 2
  for (i in 2:nt) {
    if (tt[i] <= tend) { ## status unchanged
      if (status == 1)   ## moving
        {
        x[i] <- sim_BB[k]
        k <- k + 1
        }
      else               ## resting
        {x[i] <- x[i - 1]}
    }
    if (tt[i] == tend) {
      status <- 1 - status ## switch status
      j <- j + 1
      tend <- brtimes[j]
    }
  }
  x[tt %in% time]
}


#' Sampling from a Moving-Resting bridge
#'
#' A moving-resting process consists of two states: moving and resting.
#' The transition between the two states is modeled by an alternating
#' renewal process, with exponentially distributed duration. An animal
#' stays at the same location while resting, and moves according to a
#' Brownian motion while moving. A moving-resting bridge is an extension
#' of Brownian bridge that uses moving-resting process to bridge given
#' starting point and ending point.
#'
#' @param time time points at which observations are to be simulated
#' @param start_pt the start point location of Brownian bridge
#' (numeric for one dimension, vector for multiple dimension)
#' @param end_pt the end point location of Brownian brige
#' (numeric for one dimension, vector for multiple dimension)
#' @param lamM rate parameter of the exponential duration while moving
#' @param lamR rate parameter of the exponential duration while resting
#' @param sigma volatility parameter of the Brownian motion while moving
#' @param s0 the state at time 0, must be one of "m" or "r", for moving and
#' resting, respectively
#'
#' @return
#' A \code{data.frame} whose first column is the time points and whose
#' other columns are coordinates of the locations.
#' 
#' @examples
#' tgrid <- seq(0, 10, length=1001)
#' ## make it irregularly spaced
#' tgrid <- sort(sample(tgrid, 800))
#' 
#' ## a 2-dim moving-resting bridge starting from c(0, 0)
#' ## and ending at c(10, -10)
#' dat <- rMRB(tgrid, start_pt = c(0, 0), end_pt = c(10, -10),
#'             lamM = 1, lamR = 1, sigma = 1, "m")
#' par(mfrow = c(2, 1))
#' plot(dat[,1], dat[,2], xlab="t", ylab="X(t)", type='l')
#' plot(dat[,1], dat[,3], xlab="t", ylab="X(t)", type='l')
#' par(mfrow = c(1, 1))
#' 
#' @export
rMRB <- function(time, start_pt, end_pt,
                 lamM, lamR, sigma, s0) {
  time <- time - time[1]
  stopifnot(s0 %in% c("m", "r"))
  t0moving <- (s0 == "m")
  
  stopifnot(length(start_pt) == length(end_pt))
  dim = length(start_pt)
  
  lam1 <- if (t0moving) lamM else lamR
  lam2 <- if (t0moving) lamR else lamM
  
  tmax <- time[length(time)]
  brtimes <- sim1mr.times.bbz(tmax, lam1, lam2)

  coord <- matrix(NA, ncol = dim, nrow = length(time))
  for (i in seq_len(dim)) {
    coord[, i] <- sim1mrb.bbz(tmax, sigma, time, brtimes[, 1], t0moving,
                              start_pt[i], end_pt[i])
  }
  
  data.frame(time = time, coord)
}







