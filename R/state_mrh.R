### estimation of whole state path in moving-resting-handling model

##' Estimation of states at each time point with Moving-Resting-Handling Process
##'
##' Estimate the state at each time point under the Moving-Resting-Handling
##' process with Embedded Brownian Motion with animal movement data at
##' discretely time points. See the difference between \code{fitStateMRH}
##' and \code{fitViterbiMRH} in detail part.
##'
##' \code{fitStateMRH} estimates the most likely state by maximizing
##' the probability of \eqn{Pr(S(t = t_k) = s_k | X)}, where X is the whole
##' data and \eqn{s_k} is the possible sates at \eqn{t_k} (moving, resting
##' or handling).
##'
##' \code{fitViterbiMRH} estimates the most likely state path by maximizing
##' \eqn{Pr(S(t = 0) = s_0, S(t = 1) = s_1, ..., S(t = n) = s_n | X)}, where
##' X is the whole data and \eqn{s_0, s_1, ..., s_n} is the possible
##' state path.
##'
##' @param data a \code{data.frame} whose first column is the observation
##' time, and other columns are location coordinates.
##' @param theta the parameters for Moving-Resting-Handling model, in the
##' order of rate of moving, rate of resting, rate of handling, volatility
##' and switching probability.
##' @param integrControl Integration control vector includes rel.tol,
##' abs.tol, and subdivisions.
##'
##' @return A \code{data.frame} contains estimated results, with elements:
##' \itemize{
##'  \item original data be estimated.
##'  \item conditional probability of moving, resting, handling (\code{p.m},
##' \code{p.r}, \code{p.h}), which is \eqn{Pr(S(t = t_k) = s_k | X)} for
##' \code{fitStateMRH} and \eqn{log-Pr(s_0, ..., s_k | X_k)} for
##' \code{fitViterbiMRH}, where \eqn{X_k} is \eqn{(X_0, ..., X_k)}.
##'  \item estimated states with 0-moving, 1-resting, 2-handling.
##' }
##'
##' @references
##' Pozdnyakov, V., Elbroch, L.M., Hu, C., Meyer, T., and Yan, J. (2018+)
##' On estimation for Brownian motion governed by telegraph process with
##' multiple off states. <arXiv:1806.00849>
##'
##' @seealso \code{\link{rMovResHan}} for simulation.
##' \code{\link{fitMovResHan}} for estimation of parameters.
##'
##' @examples
##' \donttest{
##' set.seed(06269)
##' tgrid <- seq(0, 400, by = 8)
##' dat <- rMovResHan(tgrid, 4, 0.5, 0.1, 5, 0.8, 'm')
##' fitStateMRH(dat, c(4, 0.5, 0.1, 5, 0.8))
##' fitViterbiMRH(dat, c(4, 0.5, 0.1, 5, 0.8))
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
