##' @import foreach
nllk_composite_parallel <- function(theta, data, groupSize,
                                    integrControl, numThreads) {
    nGroup <- nrow(data) %/% groupSize

    ## create parallel backend
    cl = parallel::makeCluster(numThreads); on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)

    i = 1 #Dummy line for RStudio warnings
    result <- foreach(i = 1:nGroup) %dopar% {
        dataCart <- data[((i - 1) * groupSize + 1):(i * groupSize), ]
        nllk_fwd_ths(theta, dataCart, integrControl)
    }

    sum(unlist(result))
}

