#' Regular placement of points on surface of a sphere
#'
#' @param r 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
regular.placement <- function(r, N) {
  r <- r
  num <- N
  
  a <- 4.0 * pi * r ^ 2 / num
  d <- sqrt(a)
  m_theta  <- round(pi / d)
  d_theta <- pi / m_theta
  d_phi   <- a / d_theta
  
  vals <- lapply(0:(m_theta -1), FUN = function(m) {
    theta <- pi * (m + 0.5) / m_theta
    m_phi <- round(2.0 * pi * sin(theta) / d_phi)
    
    lapply(0:(m_phi - 1), FUN = function(n) {
      phi <- 2.0 * pi * n / m_phi
      x <- r * sin(theta) * cos(phi)
      y <- r * sin(theta) * sin(phi)
      z <- r * cos(theta)
      c(theta, phi, pi/m_theta, 2.0 * pi / m_phi)
    }) %>% do.call(rbind,.)
  }) %>% do.call(rbind, .)
  
  return(vals)
}
