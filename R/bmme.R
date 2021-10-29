##### need the Matrix package for sparse matrix operation
## require(Matrix)

#### generate data from BMME
#### Input:
####      tgrid: time grid, vector of time points
####      dim  : dimension of Brownian Motion (Wiener process) 
####      sigma: vector of dim, Brownian motion sd, recyclable
####      delta: vector of dim, measurement error sd, recyclable
#### Output:
####      a matrix of dim + 1 columns: time, location coordinates


#' Sampling from Brown Motion with Measurement Error
#'
#' Given the volatility parameters of a Brownian motion and normally
#' distributed measurement errors, generate the process at discretely
#' observed time points of a given dimension.
#'
#' @param time vector of time points at which observations are to be sampled
#' @param dim (integer) dimension of the Brownian motion
#' @param sigma volatility parameter (sd) of the Brownian motion
#' @param delta sd parameter of measurement error
#'
#' @return
#'   A \code{data.frame} whose first column is the time points and whose
#'   other columns are coordinates of the locations.
#'
#' @references
#' Pozdnyakov V., Meyer, TH., Wang, Y., and Yan, J. (2013)
#' On modeling animal movements using Brownian motion with measurement
#' error. Ecology 95(2): p247--253. doi:doi:10.1890/13-0532.1.
#'
#' @examples
#' tgrid <- seq(0, 10, length = 1001)
#' ## make it irregularly spaced
#' tgrid <- sort(sample(tgrid, 800))
#' dat <- rBMME(tgrid, 1, 1)
#' plot(dat[,1], dat[,2], xlab="t", ylab="X(t)", type="l")
#'
#' @export
rBMME <- function(time, dim = 2,  sigma = 1, delta = 1) {  
    n <- length(time)
    dat <- matrix(NA_real_, n, dim)
    t.inc.root <- sqrt(diff(time))
    sigma <- rep(sigma, length = dim)
    delta <- rep(delta, length = dim)
    for (i in 1:dim) {
        bm <- c(0, cumsum(rnorm(n - 1, 0, sd = t.inc.root * sigma[i])))
        err <- rnorm(n, 0, sd = delta[i])
        dat[,i] <- bm + err
    }
    as.data.frame(cbind(time = time, dat))
}

#' 'rBmme' is deprecated. Using new function 'rBMME' instead.
#' @rdname rBMME
#' @export
rBmme <- function(time, dim = 2,  sigma = 1, delta = 1) {
    .Deprecated("rBMME")
    rBMME(time, dim, sigma, delta)
}

#### generate the covariance matrix in sparse structure
#### Input:
####      tinc:  vector of time increment
####      param: vector of (sigma, delta)
#### Output:
####      a banded sparse covariance matrix of the increments

getSparseSigma <- function(tinc, param) {
    n <- length(tinc)
    d0 <- c(tinc * param[1]^2 + 2 * param[2]^2)   
    d1 <- rep(- param[2]^2, n - 1)
    Matrix::bandSparse(n, k = c(0, 1), diagonals = list(d0, d1), symmetric = TRUE)
}

#### multivariate density using sparse covariance matrix
#### Input:
####      x:     normal vector at which the density needs to be evaluated
####      sigl:  sigl = t(chol(Sigma)), Sigma is sparse covariance matrix
####      log:   if TRUE, return log density
#### Output:
####      log density of density
dmvnormSparse <- function(x, sigl, log = TRUE) {
    n <- NROW(x)
    ## sigl <- t(chol(Sigma))
    part1 <- - n / 2 * log (2 * pi)  - sum(log(diag(sigl)))
    y <- forwardsolve(sigl, x)
    part2 <- - sum(y^2) / 2
    if (log) return(part1 + part2)
    else return(exp(part1 + part2))
}

#### negative loglikelihood which can be fed to optim
#### Input:
####      param: vector of (sigma, delta)
####      dinc:  increment matrix of three columns: time, location_x, location_y
#### Output:
####      loglikelihood
nllk.bmme <- function(param, dinc) {
    if (any(param <= 0)) return(NaN)
    dim <- ncol(dinc) - 1
    Sigma <- getSparseSigma(dinc[,1], c(param[1], param[2]))
    sigl <- t(chol(Sigma))
    llk <- 0
    for (i in 1:dim) {
        llk <- llk + dmvnormSparse(dinc[, i + 1], sigl, log=TRUE)
    }
    return(-llk)

}

#### obtain initial value for sigma and delta by method of moment
#### Input:
####      dat: generated data matrix from gendata ignoring measurement error
####           assuming equal sigma in two directions
#### Output
####      a vector containing rough initial value of sigma 
bmme.start <- function(dat) {
    dif <- apply(dat, 2, diff)
    dim <- ncol(dat) - 1
    st <- sqrt(sum(dif[,-1]^2) / (dim * sum(2 + dif[,1])))
    c(st, st)
}

#### find the MLE by feeding negloglik.full to optim
#### Input:
####      dat:   generated data matrix from gendata
####      start: starting value of (sigma, delta)
#### Output
####      vector containing parameter estimate, standard error, and convergence code


#' Fit a Brownian Motion with Measurement Error
#'
#' Given discretely observed animal movement locations, fit a Brownian
#' motion model with measurement errors. Using \code{segment} to fit
#' part of observations to the model. A practical application of this
#' feature is seasonal analysis.
#'
#' @param data a data.frame whose first column is the observation time, and other
#'     columns are location coordinates. If \code{segment} is not \code{NULL},
#'     additional column with the same name given by \code{segment} should be
#'     included. This additional column is used to indicate which part of
#'     observations shoule be used to fit model. The value of this column can
#'     be any integer with 0 means discarding this observation and non-0 means
#'     using this obversvation. Using different non-zero numbers indicate different
#'     segments. (See vignette for more details.)
#'     
#' @param start starting value of the model, a vector of two component, one for
#'     sigma (sd of BM) and the other for delta (sd for measurement error).
#'     If unspecified (NULL), a moment estimator will be used assuming equal
#'     sigma and delta.
#' @param segment character variable, name of the column which indicates segments,
#'     in the given \code{data.frame}. The default value, \code{NULL}, means using
#'     whole dataset to fit the model.
#' @param method the method argument to feed \code{optim}.
#' @param optim.control a list of control that is passed down to \code{optim}.
#'
#' @details
#'   The joint density of the increment data is multivariate normal with a
#'   sparse (tri-diagonal) covariance matrix. Sparse matrix operation from
#'   package Matrix is used for computing efficiency in handling large data.
#'
#' @return
#'   A list of the following components:
#'   \item{estimate }{the estimated parameter vector}
#'   \item{var.est }{variance matrix of the estimator}
#'   \item{loglik }{loglikelihood evaluated at the estimate}
#'   \item{convergence}{convergence code from optim}
#'
#' @references
#' Pozdnyakov V., Meyer, TH., Wang, Y., and Yan, J. (2013)
#' On modeling animal movements using Brownian motion with measurement
#' error. Ecology 95(2): p247--253. doi:doi:10.1890/13-0532.1.
#' @seealso
#'   \code{\link{fitMR}}
#' @examples
#' set.seed(123)
#' tgrid <- seq(0, 500, by = 1)
#' dat <- rBMME(tgrid, sigma = 1, delta = 0.5)
#'
#' ## using whole dataset to fit BMME
#' fit <- fitBMME(dat)
#' fit
#'
#' ## using part of dataset to fit BMME
#' batch <- c(rep(0, 100), rep(1, 200), rep(0, 50), rep(2, 100), rep(0, 51))
#' dat.segment <- cbind(dat, batch)
#' fit.segment <- fitBMME(dat.segment, segment = "batch")
#' head(dat.segment)
#' fit.segment
#' 
#' @export
fitBMME <- function(data, start = NULL, segment = NULL,
                    method = "Nelder-Mead",
                    optim.control = list()) {
    if (is.null(segment)) {


        ## normal process
        if (is.null(start)) start <- bmme.start(data)
        dinc <- apply(data, 2, diff)
        fit <- optim(start, nllk.bmme, dinc = dinc, hessian = TRUE,
                     method=method, control = optim.control)
        ## Sigma <- getSparseSigma(data[,1], fit$par)
        ans <- list(estimate = fit$par,
                    var.est = solve(fit$hessian),
                    loglik = - fit$value,
                    convergence = fit$convergence
                    ## vmat.l = t(chol(Sigma))
                    )
        ans$data <- data
        attr(ans, "class") <- "smam_bmme"
        return(ans)

        
    } else {


        ## seasonal process
        result <- fitBMME_seasonal(data = data, segment = segment,
                                   start = start, method = method,
                                   optim.control = optim.control)
        result$data <- data
        attr(result, "class") <- "smam_bmme"
        return(result)

        
    }
}


#' 'fitBmme' is deprecated. Using new function 'fitBMME' instead.
#' @rdname fitBMME
#' @export
fitBmme <- function(data, start = NULL, method = "Nelder-Mead",
                    optim.control = list()) {
    .Deprecated("fitBMME")
    fitBMME(data = data, start = start, segment = NULL,
            method = method, optim.control = optim.control)
}

## remind Vladmir to check linear transformation of (B0, ... Bn, xi0, ..., xin)
getVarObs <- function(tgrid, param) {
    m <- length(tgrid)
    sigma2 <- param[1]^2
    delta2 <- param[2]^2
    diags <- sigma2 * tgrid + delta2
    vmat <- diag(diags)
    offdiags <- sigma2 * rep(tgrid[-m], (m-1):1)
    vmat[row(vmat) > col(vmat)] <- offdiags
    vmat <- vmat + t(vmat) - diag(diags)
    vmat
}
#### density at point x at time t, Prob(X(t) \in dx)
#### Input:
####      times: vector of times at which the density is needed
####      x:     one point at which the density is needed
####      param: vector of (sigma, delta)
####      dat:   data matrix
#### Output:
####      a vector containing the density evaluation at times for point x

dbmme.t.x <- function(times, x, param, SigmaZ.l, dat) {
    dim <- ncol(dat) - 1
    stopifnot(length(x) == dim)
    sigma2 <- param[1]^2
    delta2 <- param[2]^2
    tgrid <- dat[,1]
    sig12 <- sigma2 * outer(times, tgrid, pmin) 
    SigmaZ.l.sig21 <- forwardsolve(SigmaZ.l, t(sig12))
    dens <- 1
    v.t <- sigma2 * times - colSums(SigmaZ.l.sig21^2)
    
    for (i in 1:dim) {
        mu.t <- colSums(SigmaZ.l.sig21 * forwardsolve(SigmaZ.l, dat[, i + 1]))
        dens <- dens * dnorm(x[i], mu.t, sqrt(v.t))
    }
    dens
}

dbb.t.x <- function(times, x, param, dat) {
    dim <- ncol(dat) - 1
    sigma2 <- param[1]^2
    tgrid <- dat[,1]
    t0 <- tgrid[1]
    t1 <- tgrid[2]
    v.t <- (t1 - times) * (times - t0) / (t1 - t0) * sigma2
    dens <- 1
    for (i in 1:dim) {
        mu.t <- dat[1, i + 1] + (times - t0) / (t1 - t0) * diff(dat[, i + 1])
        dens <- dens * dnorm(x[i], mu.t, sqrt(v.t))
    }
    dens
}

## x is a single point
dbmme.x <- function(x, tint, param, SigmaZ.l, dat) {
    tall <- tint[2] - tint[1]
    integrand <- function(u) 100 * dbmme.t.x(u, x, param, SigmaZ.l, dat)
    integrate(integrand, tint[1], tint[2])$value / tall / 100
}

dbmme.x.my <- function(x, tint, param, SigmaZ.l, dat, nstep=100) {
    tall <- tint[2] - tint[1]
    tstep <- tall / nstep
    times <- seq(tint[1]+ tstep / 2, tint[2], by = tstep)
    mean(dbmme.t.x(times, x, param, SigmaZ.l, dat))
}
    
dbb.x <- function(x, tint, param, dat) {
    tall <- tint[2] - tint[1]
    integrand <- function(u) 100 * dbb.t.x(u, x, param, dat)
    integrate(integrand, tint[1], tint[2])$value / tall / 100
}

dbb.x.my <- function(x, tint, param, dat, nstep=100) {
    tall <- tint[2] - tint[1]
    tstep <- tall / nstep
    times <- seq(tint[1]+ tstep / 2, tint[2], by = tstep)
    mean(dbb.t.x(times, x, param, dat))
}

dbmme.x.v <- function(x, tint, param, SigmaZ.l, dat) {
    if (!is.matrix(x)) x <- as.matrix(x)
    apply(x, 1, dbmme.x, tint=tint, param=param, SigmaZ.l=SigmaZ.l, dat=dat)
}

doccuTime <- function(x, tint, param, dat) {
    if (!is.matrix(x)) x <- as.matrix(x)
    SigmaZ <- getVarObs(dat[,1], param)
    SigmaZ.l <- t(chol(SigmaZ))
    doccu <- dbmme.x.v(x, tint, param, SigmaZ.l, dat)
    cbind(x, doccu)
}

