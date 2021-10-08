#' Variance-Covariance Matrix of smam Estimators
#'
#' This function calculates variance covariance matrix for
#' estimators from smam package. Different methods will
#' be used for different `smam` models.
#'
#' @name vcov
#' @param x a fitted object from one of `smam::fitXXXX` functions
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
vcov.smam_mrme <- function(x, nBS = 25, detailBS = TRUE, numThreads = 5,
                           gradMethod = "simple",
                           vcovMethod = "pBootstrap",
                           integrControl = integr.control(), ...) {
    stopifnot(vcovMethod %in% c("pBootstrap", "Godambe"))
    if (vcovMethod == "pBootstrap") {
        result <- estVarMRME_pBootstrap(x$estimate, x$data, nBS = nBS, detailBS = detailBS,
                              numThreads = numThreads, integrControl = integrControl)
    }
    if (vcovMethod == "Godambe") {
        result <- estVarMRME_Godambe(x$estimate, x$data, nBS = nBS, numThreads = numThreads,
                           gradMethod = gradMethod, integrControl = integrControl)
    }
    result
}
