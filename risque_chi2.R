# Graphiques loi chi-deux décentrée

# Estimation du paramètre de décentralité
# MLE based on Saxena & Alam's paper for unrestrited ncp
# MLE with lowerbound restriction by Éric Marchand & Fanny Rancourt

# Perte quadratique

# Auteur : Fanny Rancourt

library(rootSolve)
library(distr)

# Additional functions
source("generic.R")
source("ChiSq.R")

theta_0 <- 1 #lowerbound theta (ncp)

chi2 <- Chisq(df = 6, ncp = theta_0)

nb_ec <- 6 #integral upperbound
pas <- 0.1

max <- theta_0 + 10 # Graph param

# Maximum likelihood estimator --------------------------------------------

# lowerbound for 
bound <- df(chi2) 
if (ncp(chi2) > 0) {
  evm <-
    function(x) {
      sqrt(x / ncp(chi2)) * besselI(sqrt(x * ncp(chi2)), df(chi2) / 2) /
        besselI(sqrt(x * ncp(chi2)), df(chi2) / 2 - 1) - 1
    }
  bound <- uniroot.all(evm, interval = c(0.1, (10 * ncp(chi2))))
}

mle_SA <- function(w) {
  result <- vector(mode = "numeric", length = length(w))
  for (i in 1:length(w)) {
    if (w[i] <= bound) {
      result[i] <- ncp(chi2)
    }
    else{
      evm <- function(x) {
        sqrt(w[i] / x) * besselI(sqrt(w[i] * x), df(chi2) / 2) / besselI(sqrt(w[i] * x), df(chi2) / 2 - 1) - 1
      }
      result[i] <- uniroot.all(evm, interval = c(0.1, (10 * w[i])))
    }
  }
  return(result)
}


# psi_lowerbound ----------------------------------------------------------

plot_psi(chi2)

# graphique risque --------------------------------------------------------

# Comparaisons des riques: ESB, partie positive de l'ESB, psi*, EVM

Pts <- matrix(data = NA,
              nrow = ((max - ncp(chi2)) / pas + 1),
              ncol = 5)

for (i in 0:((max - ncp(chi2)) / pas)) {
  theta <- ncp(chi2) + i * pas
  
  chi <- Chisq(df = df(chi2), ncp = theta)
  
  total_d0 <-
    function(w) {
      (unbiased(chi, w) - theta) ^ 2 * d(chi)(w)
    }
  total_psi <-
    function(w) {
      (unbiased(chi, w) - psi_lowerbound(chi2, w) - theta) ^ 2 * d(chi)(w)
    }
  total_sa <-
    function(w) {
      (mle_SA(w) - theta) ^ 2 * d(chi)(w)
    }
  
  E_total_d0 <-
    integrate(total_d0, 0, expected_value(chi) + nb_ec * sqrt(variance(chi)))
  E_total_max <-
    integrate(total_d0,
              theta_0 + df(chi),
              expected_value(chi) + nb_ec * sqrt(variance(chi)))
  
  # Vectorize helps integration
  E_total_psi <-
    integrate(Vectorize(total_psi), 0, expected_value(chi) + nb_ec * sqrt(variance(chi)), subdivisions = 10000L)
  E_total_sa <-
    integrate(total_sa,
              bound + .1,
              expected_value(chi) + nb_ec * sqrt(variance(chi)))
  
  Pts[i + 1, 1] <- theta
  Pts[i + 1, 2] <- E_total_psi$value #psi
  Pts[i + 1, 3] <- E_total_max$value + (ncp(chi) - theta) ^ 2 * p(chi)(ncp(chi) + df(chi)) #max of unbiased estimator
  Pts[i + 1, 4] <- E_total_d0$value #unbiased estimator
  Pts[i + 1, 5] <- E_total_sa$value + (ncp(chi) - theta) ^ 2 * pchisq(bound, df(chi), ncp = theta) #MLE
}

plot(
  Pts[, 1],
  Pts[, 2] / Pts[, 4],
  type = "l",
  xlab = "ncp",
  ylab = "ratio",
  main = paste("Noncentral chisq risk ratio, lambda_0 =", theta_0, " df =", df(chi)),
  col = 1,
  cex.main = 1
)

# Rapports risques evm/max, psi/max ---------------------------------------

leg <- c("psi/l_+", "mle/l_+")
plot(
  Pts[, 1],
  Pts[, 2] / Pts[, 3],
  type = "l",
  xlab = "ncp",
  ylab = "R(mle) / R(lambda_+)",
  main = paste("Noncentral chisq risk ratio, lambda_0 =", theta_0, " df =", df(chi)),
  col = 1,
  cex.main = 1
)
lines(Pts[, 1], Pts[, 5] / Pts[, 3], lty = 2, col = 2)

legend('topright',
       leg,
       lty = c(1, 2),
       col = c(1, 2))
