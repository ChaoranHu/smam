### vcov part

#' Variance-Covariance Matrix of smam Estimators
#'
#' This function calculates variance covariance matrix for
#' estimators from smam package. Different methods will
#' be used for different `smam` models.
#'
#' @name vcov
#' @param object a fitted object from one of `smam::fitXXXX` functions
#' @param ... Optional arguments that are not used
#' @examples
#' ## time consuming example
#' #tgrid <- seq(0, 100, length=100)
#' #set.seed(123)
#' #dat <- rMRME(tgrid, 1, 0.5, 1, 0.01, "m")
#'
#' ## fit whole dataset to the MRME model
#' #fit <- fitMRME(dat, start=c(1, 0.5, 1, 0.01))
#' #fit
#'
#' ## get covariance matrix of estimators
#' #vcov(fit)
#'
#' @importFrom stats vcov
NULL

#' @param nBS number of bootstrap.
#' @param detailBS whether or not output estimation results of bootstrap,
#'  which can be used to generate bootstrap CI. Required
#'  when `vcovMethod=='pBootstrap'`.
#' @param numThreads the number of threads for parallel computation.
#'  If its value is greater than 1, then parallel computation will be
#'  processed. Otherwise, serial computation will be processed.
#' @param gradMethod method used for numeric gradient (\code{numDeriv::grad}).
#'  Required when `vcovMethod=='Godambe'`.
#' @param vcovMethod method of calculating variance covariance matrix.
#'  This should be one of `pBootstrap` (default) and `Godambe`.
#' @param integrControl a list of control parameters for the \code{integrate}
#' function: rel.tol, abs.tol, subdivision.
#' @rdname vcov
#' @export
vcov.smam_mrme <- function(object, nBS = 25, detailBS = TRUE, numThreads = 5,
                           gradMethod = "simple",
                           vcovMethod = "pBootstrap",
                           integrControl = integr.control(), ...) {
    stopifnot(vcovMethod %in% c("pBootstrap", "Godambe"))
    if (vcovMethod == "pBootstrap") {
        result <- estVarMRME_pBootstrap(object$estimate, object$data, nBS = nBS, detailBS = detailBS,
                              numThreads = numThreads, integrControl = integrControl)
    }
    if (vcovMethod == "Godambe") {
        result <- estVarMRME_Godambe(object$estimate, object$data, nBS = nBS, numThreads = numThreads,
                           gradMethod = gradMethod, integrControl = integrControl)
    }
    result
}

#' @rdname vcov
#' @export
vcov.smam_mm <- function(object, nBS = 25, detailBS = TRUE, numThreads = 5,
                         integrControl = integr.control(), ...) {
    result <- estVarMM(object$estimate, object$data, nBS = nBS, detailBS = detailBS,
                       numThreads = numThreads, integrControl = integrControl)

    result
}

#' @rdname vcov
#' @export
vcov.smam_mrh <- function(object, numThreads = 5,
                          integrControl = integr.control(), ...) {
    result <- hessMRH(object$estimate, object$data, integrControl = integrControl,
                      numThreads = numThreads)

    solve(result)
}


#' @rdname vcov
#' @export
vcov.smam_mr <- function(object, ...) {
    object$varest
}

#' @rdname vcov
#' @export
vcov.smam_bmme <- function(object, ...) {
    object$var.est
}


### estimate part

#' Estimate Result of smam Estimators
#'
#' `estimate` function returns the estimate result
#' of `smam::fitXXXX` from smam package.
#'
#' @param x a fitted object from one of `smam::fitXXXX` functions
#' @param ... other arguments
#' @examples
#' ## time consuming example
#' #tgrid <- seq(0, 100, length=100)
#' #set.seed(123)
#' #dat <- rMRME(tgrid, 1, 0.5, 1, 0.01, "m")
#'
#' ## fit whole dataset to the MRME model
#' #fit <- fitMRME(dat, start=c(1, 0.5, 1, 0.01))
#' #fit
#' 
#' ## get covariance matrix of estimators
#' #estimate(fit)
#' 
#' @export
estimate <- function(x, ...) UseMethod("estimate")

#' @rdname estimate
#' @method estimate smam_mrme
#' @export
estimate.smam_mrme <- function(x, ...) {
    result <- x$estimate
    names(result) <- c("lamM", "lamR", "sigma", "sig_err")
    result
}

#' @rdname estimate
#' @method estimate smam_mr
#' @export
estimate.smam_mr <- function(x, ...) {
    result <- x$estimate
    names(result) <- c("lamM", "lamR", "sigma")
    result
}

#' @rdname estimate
#' @method estimate smam_mm
#' @export
estimate.smam_mm <- function(x, ...) {
    result <- x$estimate
    names(result) <- c("lam1", "lam2", "sigma1", "sigma2")
    result
}


#' @rdname estimate
#' @method estimate smam_mrh
#' @export
estimate.smam_mrh <- function(x, ...) {
    result <- x$estimate
    names(result) <- c("lamM", "lamR", "lamH", "sigma", "p")
    result
}


#' @rdname estimate
#' @method estimate smam_bmme
#' @export
estimate.smam_bmme <- function(x, ...) {
    result <- x$estimate
    names(result) <- c("sigma", "delta")
    result
}