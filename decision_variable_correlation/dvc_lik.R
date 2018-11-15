f_lik_dvc <- function(gamma_s, dprime_s, rho, dvr) {
  # Decision Variation Corrleation Likelihood Function
  #
  # Args:
  #   gamma_s: Subject bais.
  #   dprime_s: Decision Variable of the subject
  #   rho:  Decision Variable Correlation
  #   dprime_m: Dprime of the model
  #   dvr:  Decision Variable Response
  #
  # Returns:
  #   Likelihood P(A|DVR)
  x <- (gamma_s + dprime_s/2 - rho * dvr)/sqrt(1 - rho^2)
  
  pnorm(x)
}
