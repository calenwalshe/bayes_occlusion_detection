#' Compute roots of log(pr(X|A) * PR(A)) - log(Pr(X|B) * PR(B)) = 0
#'
#' @param theta 
#' @param gamma 
#'
#' @return
#' @export
#'
#' @examples
roots.polar <- function(cov1 = diag(3), cov2 = diag(3), mean1 = c(0,0,0), mean2 = c(0,0,0), theta = 0, gamma = 0, pr_1 = .5, pr_2 = .5) {
  #A <- matrix(runif(9, 0, 3), 3, 3)
  #rCov1 <- t(A) %*% A
  
  #B <- matrix(runif(9, 0, 3), 3, 3)
  #rCov2 <- t(B) %*% B
  
  #cov1 <- rCov1
  #cov2 <- rCov2
  
  cov1_inv <- solve(cov1)
  cov2_inv <- solve(cov2)
  
  #mean1 <- runif(3, -3, 3)
  #mean2 <- runif(3, -3, 3)
  
  #pr_1 <- .5
  #pr_2 <- .5
  
  x_p <- c(sin(theta) * cos(gamma), sin(theta) * sin(gamma), cos(theta))
  
  a <- (t(x_p) %*% (-1 / 2 * cov2_inv) %*% x_p) -
    (t(x_p) %*% (-1 / 2 * cov1_inv) %*% x_p)
  
  b <- t(cov2_inv %*% mean2) %*% x_p - t(cov1_inv %*% mean1) %*% x_p
      
      c <- (-1 / 2 * t(mean2) %*% cov2_inv %*% mean2 -
              (1 / 2 * log(det(cov2))) +
              log(pr_2)) - ((-1 / 2 * t(mean1) %*% cov1_inv %*% mean1) -
                              (1 / 2 * log(det(cov1))) +
                              log(pr_1))
      
      
      
      if (a == 0) {
        r_1 <- (-2 * c) / (b + sqrt(b ^ 2 - 4 * a * c))
        r_2 <- (-2 * c) / (b - sqrt(b ^ 2 - 4 * a * c))
      } else {
        r_1 <- (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
        r_2 <- (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
      }
      
      x_1 <- r_1 * sin(theta) * cos(gamma)
      y_1 <- r_1 * sin(theta) * sin(gamma)
      z_1 <- r_1 * cos(theta)
      
      x_2 <- r_2 * sin(theta) * cos(gamma)
      y_2 <- r_2 * sin(theta) * sin(gamma)
      z_2 <- r_2 * cos(theta)
      
      roots.val <- c(r_1, r_2)
}

