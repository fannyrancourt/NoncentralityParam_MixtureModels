require(distr)



# Classe : distribution
# Slots : paramètres **validity
# Méthodes : estimateur sans biais,  mean, std, evm (s'il existe), param (acces & modif), densité, cdf, quantile, rng 

# Sous-classe : famille bornée de la distribution (héritage)
# Méthodes : psi*, risque (L1, L2) de psi* et de d0

## Beta Prime Parameter
setClass("BetaPrimeParameter", 
         representation = representation(shape1 = "numeric", 
                                         shape2 = "numeric"
         ), 
         prototype = prototype(shape1 = 1, shape2 = 1, ncp = 0, name = 
                                 gettext("Parameter of a Beta prime distribution")
         ), 
         contains = "Parameter"
)

## Beta Prime Distribution
setClass("BetaPrime",
         prototype = prototype(
           r = function(n){
             x <- rbeta(n,  shape1 = 1, shape2 = 1, ncp = 0)
             return(sapply(x, function(x) x/(1+x)))
           },
           d = function(x, log = FALSE){
             dbeta(x/(x+1),  shape1 = 1, shape2 = 1, ncp = 0,
                   log = log)/(x+1)^2
           },
           p = function(q, lower.tail = TRUE, log.p = FALSE ){
             pbeta(q/(q+1),  shape1 = 1, shape2 = 1, ncp = 0,
                   lower.tail = lower.tail, log.p = log.p)
           },
           q = function(p, lower.tail = TRUE, log.p = FALSE ){
             qb <- qbeta(p,  shape1 = 1, shape2 = 1, ncp = 0,
                   lower.tail = lower.tail, log.p = log.p)
             return(qb/(1-qb))
           },
           param = new("BetaPrimeParameter"),
           .logExact = TRUE,
           .lowerExact = TRUE
         ),
         contains = "AbscontDistribution"
)

## Beta Prime Family Parameter
setClass("BetaPrimeFamilyParameter", 
         representation = representation(shape1 = "numeric", 
                                         theta_0 = "numeric", 
                                         ncp = "numeric"
         ), 
         prototype = prototype(shape1 = 1, theta_0 = .Machine$double.eps, ncp = 0, name = 
                                 gettext("Parameter of a Beta prime family")
         ), 
         contains = "BetaPrimeParameter"
)

## Beta prime family distribution
setClass("BetaPrimeFamily",
         prototype = prototype(
           r = function(n){
             x <- rbeta(n,  shape1 = 1, theta_0 = 1, ncp = 0)
             return(sapply(x, function(x) x/(1+x)))
           },
           d = function(x, log = FALSE){
             dbeta(x/(x+1),  shape1 = 1, theta_0 = 1, ncp = 0,
                   log = log)/(x+1)^2
           },
           p = function(q, lower.tail = TRUE, log.p = FALSE ){
             pbeta(q/(q+1),  shape1 = 1, theta_0 = 1, ncp = 0,
                   lower.tail = lower.tail, log.p = log.p)
           },
           q = function(p, lower.tail = TRUE, log.p = FALSE ){
             qb <- qbeta(p,  shape1 = 1, theta_0 = 1, ncp = 0,
                         lower.tail = lower.tail, log.p = log.p)
             return(qb/(1-qb))
           },
           param = new("BetaPrimeParameter"),
           .logExact = TRUE,
           .lowerExact = TRUE
         ),
         contains = "BetaPrime"
)
