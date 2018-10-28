# Criterion for unequal variance
get_criterion <- function(meanA, meanB, sdA, sdB) {
  mu_A <- meanA
  mu_B <- meanB
  sigma_A <- sdA
  sigma_B <- sdB
  
  a <- -1/(2*sigma_B^2) + 1/(2*sigma_A^2)
  b <- (mu_B/sigma_B^2 - mu_A/sigma_A^2)
  c <- -mu_B^2/(2*sigma_B^2) + mu_A^2/(2*sigma_A^2) + log(sigma_A/sigma_B)
  
  x_1 <- (-b - sqrt(b^2 - 4*a*c))/(2*a)
  x_2 <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
  
  data.frame(criterion_left = x_1, criterion_right = x_2)
}
