# Graphiques loi Fisher décentrée
# Estimation du paramètre de décentralité

# Perte quadratique

# Auteur : Fanny Rancourt

library(distr)

# Additional functions
source("generic.R")
source("Fd.R")

#######################################################
## Graphiques loi Fisher d?centr?e perte quadratique ##
#######################################################

theta_0 <- 5  #lowerbound theta (ncp)

df_1 <- 7
df_2 <- 5 # must be >4 !

pas <- 0.05
nb_ec <- 6 #nombre d'?cart-types pour la borne sup?rieure

max <- 35 #Graph param


# psi_lowerbound ----------------------------------------------------------

fisher <- Fd(df1 = df_1, df2 = df_2, ncp = theta_0)

plot_psi(fisher)

# graphique risque --------------------------------------------------------

Pts <- matrix(data = NA,
              nrow = ((max - ncp(fisher)) / pas + 1),
              ncol = 3)

for (i in 0:((max - ncp(fisher)) / pas)) {
  theta <- ncp(fisher) + i * pas
  
  fd <- Fd(df1 = df1(fisher), df2 = df2(fisher), ncp = theta)
  
  total_d0 <-
    function(w) {
      (unbiased(fd, w) - theta) ^ 2 * d(fd)(w)
    }
  total_psi <-
    function(w) {
      (unbiased(fd, w) - psi_lowerbound(fisher, w) - theta) ^ 2 * d(fd)(w)
    }
  
  E_total_d0 <-
    integrate(total_d0, 0, expected_value(fd) + nb_ec * sqrt(variance(fd)))
  
  # Vectorize helps integration
  E_total_psi <-
    integrate(Vectorize(total_psi), 0, expected_value(fd) + nb_ec * sqrt(variance(fd)), subdivisions = 10000L)
  
  Pts[i + 1, 1] <- theta
  Pts[i + 1, 2] <- E_total_psi$value #psi
  Pts[i + 1, 3] <- E_total_d0$value #unbiased estimator
}

plot(
  Pts[, 1],
  Pts[, 2] / Pts[, 3],
  type = "l",
  xlab = "ncp",
  ylab = "ratio",
  main = paste("Noncentral Fisher risk ratios, t_0 =", theta_0, " df1 =", df1(fd), "df2 =", df2(fd)),
  col = 1,
  cex.main = 1
)
