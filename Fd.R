# require(distr)

################################
##
## Class: FdDistribution
##
################################

# ?Fd for more info

# psi_lowerbound --------------------------------------------------------------------

## Moments
setMethod("expected_value", "Fd", function(object) {
  return((df1(object) + ncp(object)) / (df2(object) - 2) * df2(object) / df1(object))
})

setMethod("variance", "Fd", function(object) {
  return(2 * ((df1(object) + 2 * ncp(object)) * (df2(object) - 2) + (df1(object) + ncp(object))^2) / 
           ((df2(object) - 2)^2 * (df2(object) - 4)) * (df2(object) / df1(object))^2)
})

## Unbiased ncp estimator
setMethod("unbiased",
          signature = c("Fd", "numeric"),
          function(object, obs) {
            return((df2(object) - 2) * df1(object) / df2(object) * obs - df1(object))
          })

## Lowerbound
setMethod("psi_lowerbound",
          signature = c("Fd", "numeric"),
          function(object, obs) {
              # qte <- function(t) {t * d(object)(t)}
              num_TruncExpValue <- integrate(function(t) {t * d(object)(t)}, 0 , obs)
              if (ncp(object) > 0) {
                return(((df2(object) - 2) * df1(object) / df2(object) * num_TruncExpValue$value - 
                         ncp(object) * p(Fd(df1(object) + 2, df2(object), ncp(object)))(obs)) / p(object)(obs) - 
                         df1(object))
              } else {
                  return(- 2 * (df1(object) * obs / df2(object))^(df1(object / 2) * 
                            (1 + df1(object) * obs / df2(object))^(-df1(object) + df2(object)) / 2 + 1) / 
                            (pbeta(((df1(object) * obs / df2(object)) / 
                              (1 + df1(object) * obs / df2(object))), df1(object)/2,df_2/2) * 
                                beta(df1(object) / 2, df2(object) / 2)))
                } 
          })


## Plot psi lowerbound
setMethod("plot_psi", signature = "Fd",
          function(object, interval = c(0,40), step = 0.01, main = NULL) {
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
                  main = paste("psi Noncentral Fisher df1 =", df1(object), "df2 =", df2(object), "t_0 =", ncp(object)),
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