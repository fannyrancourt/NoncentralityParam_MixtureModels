# generic functions

if (!isGeneric("unbiased_estimator")) {
  setGeneric("unbiased_estimator", function(object) standardGeneric("unbiased_estimator"))
}

if (!isGeneric("expected_value")) {
  setGeneric("expected_value", function(object) standardGeneric("expected_value"))
}

if (!isGeneric("variance")) {
  setGeneric("variance", function(object) standardGeneric("variance"))
}

if (!isGeneric("psi_lowerbound")) {
  setGeneric("psi_lowerbound", function(object, obs) standardGeneric("psi_lowerbound"))
}

if (!isGeneric("plot_psi")) {
  setGeneric("plot_psi", function(object, ...) standardGeneric("plot_psi"))
}

if (!isGeneric("unbiased")) {
  setGeneric("unbiased", function(object, obs) standardGeneric("unbiased"))
}