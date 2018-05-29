##' Sampling from a Statistical Model of Animal Movement
##'
##' Sampling from a statistcal model of animal movement with
##' a given grid time. The specific model can be chosen from
##' Brownian motion with measurement error (\code{"bmme"}),
##' moving-resting process with embedded Brownian motion
##' (\code{"MovRes"}), and moving-resting-handling process
##' with embedded Brownian motion (\code{"MovResHan"}).
##'
##' This function \code{rSmam} is a wrapper funtion of
##' \code{rbmme(time, dim, sigma, delta)},
##' \code{rMovRes(time, lamM, lamR, sigma, s0, dim)},
##' and \code{rMovResHan(time, lamM, lamR, lamH, sigma, p, s0, dim)}.
##' The parameters should be chosen according to the corresponding
##' model.
##'
##' @param time time points at which observations are to be simulated
##' @param model specify the statistical model be used, which should be
##' within \code{c("bmme", "MovRes", "MovResHan")}.
##' @param lamM rate parameter of the exponential duration while moving
##' @param lamR rate parameter of the exponential duration while resting
##' @param lamH rate parameter of the exponential duration while handling
##' @param sigma volatility parameter of the Brownian motion while moving
##' @param delta sd parameter of measurement error
##' @param p probability of choosing resting, and 1-p is probability of choosing handling
##' @param s0 the state at time 0, must be one of "m" (moving) or "r" (resting/handling).
##' @param dim (integer) dimension of the Brownian motion
##'
##' @return
##' A \code{data.frame} whose first column is the time points and whose
##' other columns are coordinates of the locations.
##' 
##' @references
##' Pozdnyakov V., Meyer, TH., Wang, Y., and Yan, J. (2013)
##' On modeling animal movements using Brownian motion with measurement
##' error. Ecology 95(2): p247--253. doi:doi:10.1890/13-0532.1.
##' 
##' Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakov, V.,
##' Williams, S., and Meyer, T. (2014) A moving-resting process with an
##' embedded Brownian motion for animal movements.
##' Population Ecology. 56(2): 401--415.
##'
##' Pozdnyakov, V., Elbroch, L., Labarga, A., Meyer, T., and Yan, J.
##' (2017) Discretely observed Brownian motion governed by telegraph
##' process: estimation. Methodology and Computing in Applied Probability.
##' doi:10.1007/s11009-017-9547-6.
##' 
##' Pozdnyakov, V., Elbroch, L.M., Hu, C., Meyer, T., and Yan, J. (2018+)
##' On estimation for Brownian motion governed by telegraph process with
##' multiple off states. (Under Review)
##'
##' @seealso \code{\link{rbmme}}, \code{\link{rMovRes}},
##' and \code{\link{rMovResHan}} for original functions.
##'
##' @examples
##' tgrid <- seq(0, 10, length = 1001)
##' ## make it irregularly spaced
##' tgrid <- sort(sample(tgrid, 800))
##' 
##' dat.bmme <- rSmam(tgrid, model = "bmme", sigma = 0.9, delta = 0.1)
##' ## equivalent to rbmme(tgrid, sigma = 0.9, delta = 0.1)
##' 
##' dat.MovRes <- rSmam(tgrid, model = "MovRes",
##'                     lamM = 1, lamR = 0.1,
##'                     sigma = 1000, s0 = "m")
##' ## equivalent to rMovRes(tgrid, lamM = 1, lamR = 0.1,
##' ##                       sigma = 1000, s0 = "m")
##' 
##' dat.MovResHan <- rSmam(tgrid, model = "MovResHan",
##'                        lamM = 4, lamR = 0.4, lamH = 0.1,
##'                        sigma = 1000, p = 0.8, s0 = "m")
##' ## equivalent to rMovResHan(tgrid, lamM = 4, lamR = 0.4, lamH = 0.1,
##' ##                          sigma = 1000, p = 0.8, s0 = "m")
##'
##' 
##' @export
rSmam <- function(time, model = c("bmme", "MovRes", "MovResHan"), lamM, lamR, lamH,
                  sigma, delta, p, s0, dim = 2) {
    stopifnot(model %in% c("bmme", "MovRes", "MovResHan"))

    result <- switch(model,
                     bmme = as.data.frame(rbmme(time = time, dim = dim,
                                                sigma = sigma, delta = delta)),
                     MovRes = rMovRes(time = time, lamM = lamM, lamR = lamR,
                                      sigma = sigma, s0 = s0, dim = dim),
                     MovResHan = rMovResHan(time = time, lamM = lamM, lamR = lamR, lamH = lamH,
                                            sigma = sigma, p = p, s0 = s0, dim = dim))
    result
}





##' Fit a Statistical Model of Animal Movement
##'
##' Fit a statistcal model of animal movement, which can be
##' Brownian motion with measurement error (\code{"bmme"}),
##' moving-resting process with embedded Brownian motion
##' (\code{"MovRes"}), or moving-resting-handling process
##' with embedded Brownian motion (\code{"MovResHan"}).
##'
##' This function \code{fitSmam} is a wrapper funtion of
##' \code{fitBmme}, \code{fitMovRes}, \code{fitMovResHan},
##' \code{fitMovResHan.parallel}, \code{fitBmme.seasonal},
##' \code{fitMovRes.seasonal}, and \code{fitMovResHan.seasonal}.
##' The parameters should be chosen according to the corresponding
##' function.
##' 
##' @param data The dataset should be fitted, which is a \code{data.frame} whose
##' first column is the observation time, and other columns are location coordinates.
##' If it is a \code{List} with the same format as the output of \code{seasonFilter},
##' the seasonal analysis will be processed.
##' @param model specify the statistical model be used, which should be
##' within \code{c("bmme", "MovRes", "MovResHan")}.
##' @param start starting value of the corresponding model. For bmme model, it should
##' be a vector of two component, one for sigma (sd of BM) and the other for delta
##' (sd for measurement error). If unspecified (NULL), a moment estimator will be used
##' assuming equal sigma and delta. For moving-resting model, it should be a vector of
##' three components in the order of rate for moving, rate for resting, and volatility.
##' For moving-resting-handling model, it should be a vector in order of rate
##' of moving, rate of resting, rate of handling, volatility and switching
##' probability.
##' @param likelihood a character string specifying the likelihood type to
##' maximize in estimation. This can be "full" for full likelihood or
##' "composite' for composite likelihood.
##' @param logtr logical, if TRUE parameters are estimated on the log scale.
##' @param lower,upper Lower and upper bound for optimization.
##' @param method the method argument to feed \code{optim}.
##' @param optim.control,... a list of control to be passed to \code{optim}.
##' @param integrControl a list of control parameters for the \code{integrate}
##' function: rel.tol, abs.tol, subdivision.
##' @param numThreads int, the number of threads allocated for parallel
##' computation (for moving-resting-handling model). The default setup
##' is 3/4 available threads. (If this parameter is set to 1, serial algorithm
##' will be processed.)
##'
##' @return A list of estimation result with elements corresponding to selected
##' model.
##'
##' @references
##' Pozdnyakov V., Meyer, TH., Wang, Y., and Yan, J. (2013)
##' On modeling animal movements using Brownian motion with measurement
##' error. Ecology 95(2): p247--253. doi:doi:10.1890/13-0532.1.
##' 
##' Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakov, V.,
##' Williams, S., and Meyer, T. (2014) A moving-resting process with an
##' embedded Brownian motion for animal movements.
##' Population Ecology. 56(2): 401--415.
##'
##' Pozdnyakov, V., Elbroch, L., Labarga, A., Meyer, T., and Yan, J.
##' (2017) Discretely observed Brownian motion governed by telegraph
##' process: estimation. Methodology and Computing in Applied Probability.
##' doi:10.1007/s11009-017-9547-6.
##' 
##' Pozdnyakov, V., Elbroch, L.M., Hu, C., Meyer, T., and Yan, J. (2018+)
##' On estimation for Brownian motion governed by telegraph process with
##' multiple off states. (Under Review)
##'
##' @seealso \code{\link{fitBmme}}, \code{\link{fitMovRes}},
##' \code{\link{fitMovResHan}}, \code{\link{fitMovResHan.parallel}},
##' \code{\link{fitBmme.seasonal}}, \code{\link{fitMovRes.seasonal}},
##' and \code{\link{fitMovResHan.seasonal}} for original functions.
##'
##' @examples
##' set.seed(06269)
##' tgrid <- seq(0, 400, by = 8)
##'
##' ## fit bmme model
##' dat.bmme <- rSmam(tgrid, model = "bmme", sigma = 0.9, delta = 0.1)
##' fitSmam(dat.bmme, model = "bmme")
##'
##' ## fit moving-resting model
##' dat.movres <- rSmam(tgrid, model = "MovRes",
##'                     lamM = 1, lamR = 0.1,
##'                     sigma = 1000, s0 = "m")
##' fitSmam(dat.movres, model = "MovRes", start = c(1, 0.1, 1000),
##'         likelihood = "full")
##'
##' ## fit moving-resting-handling model
##' \donttest{
##' ## slow work, may take several hours
##' dat.movreshan <- rSmam(tgrid, model = "MovResHan",
##'              lamM = 4, lamR = 0.5, lamH = 0.1,
##'              sigma = 5, p = 0.8, s0 = "m")
##' fitSmam(dat.movreshan, model = "MovResHan", start = c(4, 0.5, 0.1, 5, 0.8))
##' }
##'
##' @export
fitSmam <- function(data, model = c("bmme", "MovRes", "MovResHan"),
                    start = NULL, likelihood = c("full", "composite"),
                    logtr = FALSE,
                    lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
                    upper = c(   10,    10,    10,    10, 0.999),
                    method = "Nelder-Mead",
                    optim.control = list(),
                    ...,
                    integrControl = integr.control(),
                    numThreads = NULL) {
    stopifnot(model %in% c("bmme", "MovRes", "MovResHan"))
    if (is.null(numThreads)) numThreads <- RcppParallel::defaultNumThreads() * 3 / 4
    
    if (is.data.frame(data)) {
        result <- switch(model,
                         bmme = fitBmme(data = data, start = start,
                                        method = method, optim.control = optim.control),
                         MovRes = fitMovRes(data = data, start = start,
                                            likelihood = likelihood, logtr = logtr,
                                            method = method, optim.control = optim.control,
                                            integrControl = integrControl),
                         MovResHan = ifelse(numThreads == 1,
                                            fitMovResHan(data = data, start = start,
                                                         lower = lower, upper = upper,
                                                         integrControl = integrControl),
                                            fitMovResHan.parallel(data = data, start = start,
                                                                  lower = lower, upper = upper,
                                                                  numThreads = numThreads,
                                                                  integrControl = integrControl)))
    } else {
        result <- switch(model,
                         bmme = fitBmme.seasonal(data = data, start = start,
                                                 method = method, ...),
                         MovRes = fitMovRes.seasonal(data = data, start = start,
                                                     likelihood = likelihood, logtr = logtr,
                                                     method = method, optim.control = optim.control,
                                                     integrControl = integrControl),
                         MovResHan = fitMovResHan.seasonal(data = data, start = start,
                                                           lower = lower, upper = upper,
                                                           numThreads = numThreads,
                                                           integrControl = integrControl))
    }

    result
}



