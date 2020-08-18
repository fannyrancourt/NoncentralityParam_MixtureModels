# require(distr)

################################
##
## Class: ChisqDistribution
##
################################

# ?Chisq for more info

# psi_lowerbound --------------------------------------------------------------------

## Moments
setMethod("expected_value", "Chisq", function(object) {
  return(df(object) + ncp(object))
})

setMethod("variance", "Chisq", function(object) {
  return(2 * df(object) + 4 * ncp(object))
})

## Unbiased ncp estimator
setMethod("unbiased",
          signature = c("Chisq", "numeric"),
          function(object, obs) {
            return(obs - df(object))
          })

## Lowerbound
setMethod("psi_lowerbound",
          signature = c("Chisq", "numeric"),
          function(object, obs) {
            if (ncp(object) > 0) {
              num_TruncExpValue <- integrate(function(t) {t * d(object)(t)}, 0, obs, subdivisions = 10000L)
              
              return((num_TruncExpValue$value - ncp(object) * p(Chisq(df = df(object)+2, ncp = ncp(object)))(obs)) 
                     / p(object)(obs) - df(object))
            }
            else{
              return(pgamma(obs, df(object) / 2 + 1,scale = 2) / pgamma(obs, df(object) / 2, scale = 2) * df(object) -
                       df(object))
            }
          })

## Plot psi lowerbound
setMethod("plot_psi", signature = "Chisq",
          function(object, interval = c(0,20), step = 0.01, main = NULL) {
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
                  main = paste("psi Noncentral ChiSq df =", df(object), "t_0 =", ncp(object)),
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