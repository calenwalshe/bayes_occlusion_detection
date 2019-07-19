#' Compute roots of log(pr(X|A) * PR(A)) - log(Pr(X|B) * PR(B)) = 0
#'
#' @param theta 
#' @param gamma 
#'
#' @return
#' @export
#'
#' @examples
roots.polar <- function(cov1 = diag(3), cov2 = diag(3), cov1_inv = diag(3), cov2_inv = diag(3), mean1 = c(0,0,0), mean2 = c(0,0,0), pr_1 = .5, pr_2 = .5) {
  const <- (-1 / 2 * t(mean2) %*% cov2_inv %*% mean2 -
              (1 / 2 * log(det(cov2))) +
              log(pr_2)) - ((-1 / 2 * t(mean1) %*% cov1_inv %*% mean1) -
                              (1 / 2 * log(det(cov1))) +
                              log(pr_1))
  
  b_const <- (t(cov2_inv %*% mean2) - t(cov1_inv %*% mean1))
  
  A_2 <- (-1 / 2 * cov2_inv)
  A_1 <- (-1 / 2 * cov1_inv)

  f <- function(theta, gamma) {
    x_p <- c(sin(theta) * cos(gamma), sin(theta) * sin(gamma), cos(theta))
    
    a <- t(x_p) %*% A_2 %*% x_p -
      t(x_p) %*% A_1 %*% x_p
    
    b <- (t(cov2_inv %*% mean2) %*% x_p) - (t(cov1_inv %*% mean1) %*% x_p)

    if (a == 0) {
      r_1 <- (-2 * const) / (b + sqrt(b ^ 2 - 4 * a * const))
      r_2 <- (-2 * const) / (b - sqrt(b ^ 2 - 4 * a * const))
    } else {
      r_1 <- (-b + sqrt(b ^ 2 - 4 * a * const)) / (2 * a)
      r_2 <- (-b - sqrt(b ^ 2 - 4 * a * const)) / (2 * a)
    }
    
    # case when discriminant == 0
    if((b ^ 2 - 4 * a * const) == 0) {
      r_2 <- Inf
    }
    
    roots.val <- c(theta, gamma, r_1, r_2)  
  }
}

