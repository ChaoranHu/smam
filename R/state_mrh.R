### estimation of whole state path in moving-resting-handling model

##' Estimation of states at each time point with Moving-Resting-Handling Process
##'
##' Estimate the state at each time point under the Moving-Resting-Handling
##' process with Embedded Brownian Motion with animal movement data at
##' discretely time points. See the difference between \code{fitStateMRH}
##' and \code{fitViterbiMRH} in detail part. Using \code{fitPartialViterbiMRH}
##' to estimate the state during a small piece of time interval.
##'
##' \code{fitStateMRH} estimates the most likely state by maximizing
##' the probability of \eqn{Pr(S(t = t_k) = s_k | X)}, where X is the whole
##' data and \eqn{s_k} is the possible sates at \eqn{t_k} (moving, resting
##' or handling).
##'
##' \code{fitViterbiMRH} estimates the most likely state path by maximizing
##' \eqn{Pr(S(t = t_0) = s_0, S(t = t_1) = s_1, ..., S(t = t_n) = s_n | X)}, where
##' X is the whole data and \eqn{s_0, s_1, ..., s_n} is the possible
##' state path.
##'
##' \code{fitPartialViterbiMRH} estimates the most likely state path of
##' a small peice of time interval, by maximizing the probability of
##' \eqn{Pr(S(t = t_k) = s_k, ..., S(t = t_{k+q-1}) = s_{k+q-1} | X)},
##' where \eqn{k} is the start time point and \eqn{q} is the length of interested
##' time interval.
##'
##' @param data a \code{data.frame} whose first column is the observation
##' time, and other columns are location coordinates.
##' @param theta the parameters for Moving-Resting-Handling model, in the
##' order of rate of moving, rate of resting, rate of handling, volatility
##' and switching probability.
##' @param integrControl Integration control vector includes rel.tol,
##' abs.tol, and subdivisions.
##' @param startpoint Start time point of interested time interval.
##' @param pathlength the length of interested time interval.
##'
##' @return A \code{data.frame} contains estimated results, with elements:
##' \itemize{
##'  \item original data be estimated.
##'  \item conditional probability of moving, resting, handling (\code{p.m},
##' \code{p.r}, \code{p.h}), which is \eqn{Pr(S(t = t_k) = s_k | X)} for
##' \code{fitStateMRH}; \eqn{log-Pr(s_0, ..., s_k | X_k)} for
##' \code{fitViterbiMRH}, where \eqn{X_k} is \eqn{(X_0, ..., X_k)};
##' and \eqn{log-Pr(s_k, ..., s_{k+q-1}|X)} for \code{fitPartialViterbiMRH}.
##'  \item estimated states with 0-moving, 1-resting, 2-handling.
##' }
##'
##' @references
##' Pozdnyakov, V., Elbroch, L.M., Hu, C., Meyer, T., and Yan, J. (2018+)
##' On estimation for Brownian motion governed by telegraph process with
##' multiple off states. <arXiv:1806.00849>
##'
##' @seealso \code{\link{rMRH}} for simulation.
##' \code{\link{fitMRH}} for estimation of parameters.
##'
##' @examples
##' \dontrun{
##' ## time consuming example
##' set.seed(06269)
##' tgrid <- seq(0, 400, by = 8)
##' dat <- rMRH(tgrid, 4, 0.5, 0.1, 5, 0.8, 'm')
##' fitStateMRH(dat, c(4, 0.5, 0.1, 5, 0.8))
##' fitViterbiMRH(dat, c(4, 0.5, 0.1, 5, 0.8))
##' fitPartialViterbiMRH(dat, c(4, 0.5, 0.1, 5, 0.8), 20, 10)
##' }
##' 
##' @author Chaoran Hu
##' @export
fitStateMRH <- function(data, theta,
                        integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    ncol_data <- ncol(data)

    cart_result <- fwd_bwd_ths(theta, dinc, integrControl)
    cart_fwd <- cart_result[, 1:3]
    cart_bwd <- cart_result[, 4:6]
    cart_result <- cart_fwd * cart_bwd
    cart_state <- apply(cart_result, 1, which.max) - 1
    result <- cbind(data, cart_result, cart_state)
    colnames(result)[(ncol_data+1):(ncol_data+4)] <- c("p.m", "p.r", "p.h", "states")
    as.data.frame(result)
}


##' @rdname fitStateMRH
##' @export
fitViterbiMRH <- function(data, theta,
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    ncol_data <- ncol(data)

    cart_result <- viterbi_ths(theta, dinc, integrControl)
    cart_state <- apply(cart_result, 1, which.max) - 1
    result <- cbind(data, cart_result, cart_state)
    colnames(result)[(ncol_data+1):(ncol_data+4)] <- c("p.m", "p.r", "p.h", "states")
    as.data.frame(result)
}

##' @rdname fitStateMRH
##' @export
fitPartialViterbiMRH <- function(data, theta, startpoint, pathlength,
                                 integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    nrow_data <- nrow(data)
    if (startpoint < 1 | startpoint > nrow_data) stop("start time point should be within data time interval.")
    if ((startpoint + pathlength - 1) > nrow_data) stop("end time point should be within data time interval.")
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    ncol_data <- ncol(data)

    cart_result <- partial_viterbi_ths(theta, dinc, integrControl, startpoint - 1, pathlength)
    cart_state <- apply(cart_result, 1, which.max) - 1
    result <- cbind(data[startpoint:(startpoint+pathlength-1), ], cart_result, cart_state)
    colnames(result)[(ncol_data+1):(ncol_data+4)] <- c("p.m", "p.r", "p.h", "states")
    as.data.frame(result)
}
