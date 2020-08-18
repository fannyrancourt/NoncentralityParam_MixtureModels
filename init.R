# Initialize functions

# Beta prime distribution
setMethod("initialize", "BetaPrime",
          function(.Object, shape1 = 1, shape2 = 1, ncp = 0) {
            .Object@img <- new("Reals")
            .Object@param <- new("BetaPrimeParameter", shape1 = shape1, 
                                 shape2 = shape2)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
              { sapply(rbeta(n, shape1 = shape1Sub, shape2 = shape2Sub, 
                      ncp = ncpSub), function(t) t/(1-t)) },
              list(shape1Sub = shape1, shape2Sub = shape2, 
                   ncpSub = ncp)
            )
            body(.Object@d) <- substitute(
              { dbeta(x/(x+1), shape1 = shape1Sub, shape2 = shape2Sub, 
                      ncp = ncpSub, log = log)/(x+1)^2 },
              list(shape1Sub = shape1, shape2Sub = shape2, 
                   ncpSub = ncp)
            )
            body(.Object@p) <- substitute(
              { pbeta(q/(q+1), shape1 = shape1Sub, shape2 = shape2Sub, 
                      ncp = ncpSub, lower.tail = lower.tail, 
                      log.p = log.p) },
              list(shape1Sub = shape1, shape2Sub = shape2, 
                   ncpSub = ncp)
            )
            body(.Object@q) <- substitute(
              { qbeta(p, shape1 = shape1Sub, shape2 = shape2Sub, 
                      ncp = ncpSub, lower.tail = lower.tail, 
                      log.p = log.p) / (1 - qbeta(p, shape1 = shape1Sub, shape2 = shape2Sub, 
                                                  ncp = ncpSub, lower.tail = lower.tail, 
                                                  log.p = log.p)) },
              list(shape1Sub = shape1, shape2Sub = shape2, 
                   ncpSub = ncp)
            )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- FALSE
            .Object
          })