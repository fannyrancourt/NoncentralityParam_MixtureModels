# Graphiques loi beta prime 
# Estimation du 1er paramètre de forme
# Perte quadratique        

# Auteur : Fanny Rancourt

source("class.R")
source("init.R")
source("generic.R")
source("BetaPrime.R")

theta_0 <- 1 #borne inferieure espace parametrique de theta (shape1)
shape2 <- 10 #2e parametre de forme

nb_ec <- 6 #integral upperbound
pas <- 0.005 #step for theta

max <- 10 # Graph param

# psi_lowerbound --------------------------------------------------------------------

bii <- BetaPrime(theta_0, shape2)

plot_psi(bii)

# Graphique risque --------------------------------------------------------

Pts <- matrix(data = NA,
              nrow = ((max - theta_0) / pas + 1),
              ncol = 3)

for (i in 0:((max - theta_0) / pas)) {
  theta <- theta_0 + i * pas
  
  bii <- BetaPrime(theta, shape2)
  
  total_d0 <-
    function(s) {
      (unbiased(bii, s) - theta) ^ 2 * d(bii)(s)
    }
  total_psi <-
    function(s) {
      (unbiased(bii, s) - psi_lowerbound(bii, s) - theta) ^ 2 * d(bii)(s)
    }
  
  E_total_d0 <- integrate(total_d0, 0, expected_value(bii) + nb_ec * sqrt(variance(bii)))
  
  # Vectorize helps integration
  E_total_psi <- integrate(Vectorize(total_psi), 0, expected_value(bii) + nb_ec * sqrt(variance(bii)))
  
  Pts[i + 1, 1] <- theta
  Pts[i + 1, 2] <- E_total_psi$value
  Pts[i + 1, 3] <- E_total_d0$value
}

plot(
  Pts[, 1],
  Pts[, 2] / Pts[, 3],
  type = "l",
  xlab = "theta",
  ylab = "ratio",
  main = paste("Risk ratio for Beta prime", nb_ec, "std, shape2 =", shape2, " t_0 =", theta_0),
  col = 1
)
