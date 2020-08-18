# Author: Fanny Rancourt




################################
##
## Class: BetaPrimeParameter
##
################################


## Access methods

setMethod("shape1", "BetaPrimeParameter", function(object)
  object@shape1)
setMethod("shape2", "BetaPrimeParameter", function(object)
  object@shape2)


## Replace Method

setReplaceMethod("shape1", "BetaPrimeParameter",
                 function(object, value) {
                   object@shape1 <- value
                   object
                 })
setReplaceMethod("shape2", "BetaPrimeParameter",
                 function(object, value) {
                   object@shape2 <- value
                   object
                 })


setValidity("BetaPrimeParameter", function(object) {
  if (length(shape1(object)) != 1)
    stop("shape1 has to be a numeric of length 1")
  if (shape1(object) <= 0)
    stop("shape1 has to be positive")
  if (length(shape2(object)) != 1)
    stop("shape2 has to be a numeric of length 1")
  if (shape2(object) <= 0)
    stop("shape2 has to be positive")
  else
    return(TRUE)
})


################################
##
## Class: BetaPrimeDistribution
##
################################

BetaPrime <- function(shape1 = 1,
                      shape2 = 1)
  new("BetaPrime",
      shape1 = shape1,
      shape2 = shape2)

## wrapped access methods
setMethod("shape1", "BetaPrime", function(object)
  shape1(param(object)))
setMethod("shape2", "BetaPrime", function(object)
  shape2(param(object)))


## wrapped replace methods
setMethod("shape1<-", "BetaPrime",
          function(object, value)
            new("BetaPrime",
                shape1 = value,
                shape2 = shape2(object)))
setMethod("shape2<-", "BetaPrime",
          function(object, value)
            new("BetaPrime",
                shape1 = shape1(object),
                shape2 = value))

setMethod("/", c("numeric", "BetaPrime"),
          function(e1, e2) {
            if (isTRUE(all.equal(e1, 1)) &&
                isTRUE(all.equal(ncp(e2), 0)))
              return(BetaPrime(shape1 = shape2(e2), shape2 = shape1(e2)))
            else
              e1 / as(e2, "AbscontDistribution")
          })


# psi_lowerbound --------------------------------------------------------------------

## Moments
setMethod("expected_value", "BetaPrime", function(object) {
  if (shape2(param(object)) > 1) {
    return(shape1(param(object)) /
             (shape2(param(object)) - 1))
  } else {
    return(NaN)
  }
})

setMethod("variance", "BetaPrime", function(object) {
  if (shape2(param(object)) > 2) {
    return((shape1(param(object)) * (shape1(param(object)) + shape2(param(object)) - 1)) /
      ((shape2(param(
        object
      )) - 1) ^ 2 * (shape2(param(object)) - 2)))
  } else {
    return(NaN)
  }
})

## Unbiased shape1 estimator
setMethod("unbiased",
          signature = c("BetaPrime", "numeric"),
          function(object, obs) {
            if (shape2(param(object)) > 1) {
              return((shape2(param(object)) - 1) * obs)
            } else {
              return(NaN)
            }
          })

## psi lowerbound
setMethod("psi_lowerbound",
          signature = c("BetaPrime", "numeric"),
          function(object, obs) {
            # with beta representation
            # s ~ Beta(a,b), then s/(s+1) ~ BetaPrime(a,b), a, b > 0
            
            -(dbeta(obs / (obs + 1), shape1(object) + 1, shape2(object) - 1) * (1 - obs / (obs + 1)) / 
                ((shape2(object) - 1) * pbeta(obs / (obs + 1), shape1(object), shape2(object), lower.tail = T)))
          })

## Plot psi lowerbound
setMethod("plot_psi", signature = "BetaPrime",
          function(object, interval = c(0,10), step = 0.01, main = NULL) {
            if (step > 0 & !is.na(max(interval))) {
              
              min <- min(interval)
              if (min < 0) {
                min <- 0
                warning("Interval lowerbound sets to 0")
              }
              
              pts <- matrix(data = NA,
                            nrow = ((max(interval) - min) / step) + 1,
                            ncol = 2)
              
              
              for (i in 1:((max(interval) - min) / step)) {
                t <- min + i * step
                pts[i + 1, 1] <- t
                pts[i + 1, 2] <- psi_lowerbound(object, t)
              }
              
              if (is.null(main)){
                return(plot(
                  pts[, 1],
                  pts[, 2],
                  type = "l",
                  xlab = "t",
                  ylab = "psi",
                  main = paste("psi Beta Prime shape2 =", shape2(object), "t_0 =", shape1(object)),
                  col = 1
                ))
              } else {
                return(plot(
                  pts[, 1],
                  pts[, 2],
                  type = "l",
                  xlab = "t",
                  ylab = "psi",
                  main = main,
                  col = 1
                ))
              }
              
            } else if (step <= 0) {
              stop("Needs a strictly positive step.")
            } else {
              stop("Incorrect interval")
            }
          })