### estimation of whole state path in moving-resting model

predM <- function(p.m, p.r) { ## internal func. for ratio of prob: fb
    length <- length(p.m)
    result <- numeric(length)

    for (i in 1:length) {
        if (is.infinite(p.m[i]) & p.m[i] < 0) {
            result[i] <- 0
            next
        }
        if (is.infinite(p.r[i]) & p.r[i] < 0) {
            result[i] <- 1
            next
        }
        result[i] <- p.m[i] / (p.m[i] + p.r[i])
    }
    result
}

predMV <- function(p.m, p.r) { ## internal func. for ratio of prob: viterbi
    length <- length(p.m)
    result <- numeric(length)

    for (i in 1:length) {
        if (is.infinite(p.m[i]) & p.m[i] < 0) {
            result[i] <- 0
            next
        }
        if (is.infinite(p.r[i]) & p.r[i] < 0) {
            result[i] <- 1
            next
        }
        result[i] <- 1 / (1 + exp(p.r[i] - p.m[i]))
    }
    result
}

##' Estimation of states at each time point with Moving-Resting Process
##'
##' Estimate the state at each time point under the Moving-Resting
##' process with Embedded Brownian Motion with animal movement data at
##' discretely time points. See the difference between \code{fitStateMR}
##' and \code{fitViterbiMR} in detail part. Using \code{fitPartialViterbiMR}
##' to estimate the state within a small piece of time interval.
##'
##' \code{fitStateMR} estimates the most likely state by maximizing
##' the probability of \eqn{Pr(S(t = t_k) = s_k | X)}, where X is the whole
##' data and \eqn{s_k} is the possible sates at \eqn{t_k} (moving, resting).
##'
##' \code{fitViterbiMR} estimates the most likely state path by maximizing
##' \eqn{Pr(S(t = t_0) = s_0, S(t = t_1) = s_1, ..., S(t = t_n) = s_n | X)}, where
##' X is the whole data and \eqn{s_0, s_1, ..., s_n} is the possible
##' state path.
##'
##' \code{fitPartialViterbiMR} estimates the most likely state path of
##' a small peice of time interval, by maximizing the probability of
##' \eqn{Pr(S(t = t_k) = s_k, ..., S(t = t_{k+q-1}) = s_{k+q-1} | X)},
##' where \eqn{k} is the start time point and \eqn{q} is the length of interested
##' time interval.
##'
##' @param data a \code{data.frame} whose first column is the observation
##' time, and other columns are location coordinates.
##' @param theta the parameters for Moving-Resting model, in the
##' order of rate of moving, rate of resting, volatility.
##' @param cutoff the cut-off point for prediction.
##' @param integrControl Integration control vector includes rel.tol,
##' abs.tol, and subdivisions.
##' @param startpoint Start time point of interested time interval.
##' @param pathlength the length of interested time interval.
##'
##' @return A \code{data.frame} contains estimated results, with elements:
##' \itemize{
##'  \item original data be estimated.
##'  \item conditional probability of moving, resting (\code{p.m},
##' \code{p.r}), which is \eqn{Pr(S(t = t_k) = s_k | X)} for
##' \code{fitStateMR}; \eqn{log-Pr(s_0, ..., s_k | X_k)} for
##' \code{fitViterbiMR}, where \eqn{X_k} is \eqn{(X_0, ..., X_k)};
##' and \eqn{log-Pr(s_k, ..., s_{k+q-1}|X)} for \code{fitPartialViterbiMR}.
##'  \item estimated states with 1-moving, 0-resting.
##' }
##'
##'
##' @seealso \code{\link{rMR}} for simulation.
##' \code{\link{fitMR}} for estimation of parameters.
##'
##' @examples
##' set.seed(06269)
##' tgrid <- seq(0, 400, by = 8)
##' dat <- rMR(tgrid, 4, 3.8, 5, 'm')
##' fitStateMR(dat, c(4, 3.8, 5), cutoff = 0.5)
##' fitViterbiMR(dat, c(4, 3.8, 5), cutoff = 0.5)
##' fitPartialViterbiMR(dat, c(4, 3.8, 5), cutoff = 0.5, 20, 10)
##'
##' @author Chaoran Hu
##' @export
fitStateMR <- function(data, theta, cutoff = 0.5,
                       integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    ncol_data <- ncol(data)

    cart_result <- fwd_bwd_mr(theta, dinc, integrControl)
    cart_fwd <- cart_result[, 1:2]
    cart_bwd <- cart_result[, 3:4]
    cart_result <- cart_fwd * cart_bwd
    predM <- predM(cart_result[, 1], cart_result[, 2])
    cart_state <- ifelse(predM > cutoff, 1, 0)
    result <- cbind(data, cart_result, cart_state)
    colnames(result)[(ncol_data+1):(ncol_data+3)] <- c("p.m", "p.r", "states")
    as.data.frame(result)
}


##' @rdname fitStateMR
##' @export
fitViterbiMR <- function(data, theta, cutoff = 0.5,
                         integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    ncol_data <- ncol(data)

    cart_result <- viterbi_mr(theta, dinc, integrControl)
    predM <- predMV(cart_result[, 1], cart_result[, 2])
    cart_state <- ifelse(predM > cutoff, 1, 0)
    result <- cbind(data, cart_result, cart_state)
    colnames(result)[(ncol_data+1):(ncol_data+3)] <- c("p.m", "p.r", "states")
    as.data.frame(result)
}

##' @rdname fitStateMR
##' @export
fitPartialViterbiMR <- function(data, theta, cutoff = 0.5,
                                startpoint, pathlength,
                                integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    nrow_data <- nrow(data)
    if (startpoint < 1 | startpoint > nrow_data) stop("start time point should be within data time interval.")
    if ((startpoint + pathlength - 1) > nrow_data) stop("end time point should be within data time interval.")
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    ncol_data <- ncol(data)

    cart_result <- partial_viterbi_mr(theta, dinc, integrControl, startpoint - 1, pathlength)
    predM <- predMV(cart_result[, 1], cart_result[, 2])
    cart_state <- ifelse(predM > cutoff, 1, 0)
    result <- cbind(data[startpoint:(startpoint+pathlength-1), ], cart_result, cart_state)
    colnames(result)[(ncol_data+1):(ncol_data+3)] <- c("p.m", "p.r", "states")
    as.data.frame(result)
}



## test code
## ## judge state by odds ratio
## giveResult.mr <- function(result, theta) {
##     lam1 <- theta[1]
##     lam0 <- theta[2]
##     odds <- lam0 / lam1 ## odds = p_m / p_r
    
##     myodds <- result[, 1] / result[, 2] ## estimated odds

##     ifelse(myodds > odds, 1, 0)
## }

## fitStateMR.beta <- function(data, theta,
##                        integrControl = integr.control()) {
##     if (!is.matrix(data)) data <- as.matrix(data)
##     dinc <- apply(data, 2, diff)
##     integrControl <- unlist(integrControl)
##     ncol_data <- ncol(data)

##     cart_result <- fwd_bwd_mr(theta, dinc, integrControl)
##     cart_fwd <- cart_result[, 1:2]
##     cart_bwd <- cart_result[, 3:4]
##     cart_result <- cart_fwd * cart_bwd
##     cart_state <- giveResult.mr(cart_result, theta)
##     result <- cbind(data, cart_result, cart_state)
##     colnames(result)[(ncol_data+1):(ncol_data+3)] <- c("p.m", "p.r", "states")
##     as.data.frame(result)
## }

