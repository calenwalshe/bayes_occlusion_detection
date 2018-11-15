dvc <- function(x, y, z) {
  library(purrr)
  
  response_label <-
    matrix(c('cr', 'miss', 'fa', 'hit'),
           ncol = 2,
           byrow = TRUE)
  
  response_data <- data.frame(response = x,
                              category = y,
                              dvr = z)
  
  response_data$response_type <-
    unlist(map2(x, y, function(x, y) {
      response_label[x + 1, y + 1]
    }))
  
  n_present <- sum(response_data$category == 1)
  n_absent  <- sum(response_data$category == 0)
  
  pr_hit <- sum(x[y==1 & x== 1]) / sum(y == 1)
  pr_fa <- sum(x[y == 0 & x == 1]) / sum(y == 0)
  
  dprime_h <- qnorm(pr_hit) - qnorm(pr_fa)
  browser()
  gamma_h    <- -(dprime_h / 2 + qnorm(pr_fa))
  
  z[y == 1] <- (z[y == 1] - mean(z[y == 1])) / sd(z[y == 1])
  z[y == 0] <- (z[y == 0] - mean(z[y == 0])) / sd(z[y == 0])
  
  dvr_AA <- z[y == 0 & x == 0]
  dvr_AB <- z[y == 0 & x == 1]
  dvr_BA <- z[y == 1 & x == 0]
  dvr_BB <- z[y == 1 & x == 1]
  
  f_A <- function(rho_hat) {
    lik <-
      -sum(
        log(f_lik_dvc(gamma_h, dprime_h, rho = rho_hat, dvr = dvr_AA)),
        log(1 - f_lik_dvc(
          gamma_h, dprime_h, rho = rho_hat, dvr = dvr_AB
        )))
  }
  
  dvcA <- optimize(f_A, lower = -1, upper = 1)$minimum
  
  f_B <- function(rho_hat) {
    lik <-
      -sum(
        log(1 - f_lik_dvc(gamma_h, -dprime_h, rho = rho_hat, dvr = dvr_BB)),
        log(f_lik_dvc(
          gamma_h, -dprime_h, rho = rho_hat, dvr = dvr_BA
        )))
  }
  
  dvcB <- optimize(f_B, lower = -1, upper = 1)$minimum
  
  tibble(dprime_h = dprime_h, gamma_h = gamma_h, dvcA = dvcA, dvcB = dvcB)
}
  