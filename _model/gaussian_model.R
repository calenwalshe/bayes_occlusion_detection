get_scaled_model <- function(model.wide, scale.val = 1) {
  library(dplyr)
  library(broom)
  library(bbmle)
  library(R2Cuba)
  library(multidplyr)
  library(mvtnorm)
  library(purrr)
  
  get_discriminant <- function(mean1, mean2, cov1, cov2) {
    g1 <- function(distance) {
      x <- mean1 + distance * (mean2 - mean1)
      d1 <- t(x) %*% (-1/2 * solve(cov1)) %*% x + t(solve(cov1) %*% mean1) %*% x + 
        (-1/2 * t(mean1) %*% solve(cov1) %*% mean1) - 1/2*det(cov1) + log(.5)
      
      d2 <- t(x) %*% (-1/2 * solve(cov2)) %*% x + t(solve(cov2) %*% mean2) %*% x + 
        (-1/2 * t(mean2) %*% solve(cov2) %*% mean2) - 1/2*det(cov2) + log(.5)
      
      return(abs(d1 - d2))
    }
  }
  
  integrate_distribution <- function(scale = 1, center.val, mean1, mean2, cov1, cov2) {
    inv_cov1 <- solve(cov1)
    inv_cov2 <- solve(cov2)
    
    scalar.dist1 <- 1/ sqrt(det(2 * pi * cov1))
    scalar.dist2 <- 1/ sqrt(det(2 * pi * cov2))
    f <- function(x) {
      gauss_f1 <- scalar.dist1 * exp(-.5*(t(x - mean1) %*% inv_cov1 %*% (x - mean1)) - scale)
      gauss_f2 <- scalar.dist2 * exp(-.5*(t(x - mean2) %*% inv_cov2 %*% (x - mean2)) - scale)
      min(gauss_f1,gauss_f2)
    }
    int.p <- 2
    results <- vegas(3, 1, f, 
            lower = center.val - int.p * sqrt((cov1[diag(3)==1] + cov2[diag(3)==1])/2), 
            upper = center.val + int.p * sqrt((cov1[diag(3)==1] + cov2[diag(3)==1])/2),
            flags = list(verbose = 0),
            min.eval = 10^6,
            max.eval = 2*10^6)
    return(results$value)
  }
  
  cluster <- create_cluster(cores = 16)
  cluster_assign_value(cluster, 'mle2', mle2)
  cluster_assign_value(cluster, 'integrate_distribution', integrate_distribution)
  cluster_library(cluster, 'purrr')
  cluster_library(cluster, 'R2Cuba')
  cluster_library(cluster, 'bbmle')
  set_default_cluster(cluster)

  scaled.model <- model.wide %>% 
    arrange(TARGET, BIN, eccentricity) %>%
    rowwise() %>%
    mutate(n_obs = length(data_response_vec_0[,1]),
           SUBJECT = "model",
           find_bound        = list(optim(0, get_discriminant(data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1), lower = 0, upper = 1, method = "Brent")),
           max_LLR    = list(data_mean_vec_0 + find_bound$par * (data_mean_vec_1 - data_mean_vec_0)),
           scale_gauss = dmvnorm(max_LLR, data_mean_vec_1, data_cov_mat_1, log = T)) %>%
    select(TARGET, BIN, scale_gauss, max_LLR, eccentricity, data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1, SUBJECT) %>%
    group_by(TARGET, BIN, eccentricity) %>%
    nest() %>%
    partition(TARGET, BIN, eccentricity) %>%
    mutate(integrate_min = 
             map(data, function(x) 
               {
               integrate_distribution(x$scale_gauss[[1]], x$max_LLR[[1]], x$data_mean_vec_0[[1]], x$data_mean_vec_1[[1]], x$data_cov_mat_0[[1]], x$data_cov_mat_1[[1]])
             }
             )) %>%
    collect() %>%
    unnest() %>%
    mutate(dprime = -2*qnorm(log(integrate_min) + scale_gauss, 0, 1, log.p = T)) %>%
    arrange(TARGET, BIN, eccentricity)

  return(scaled.model)
}

#' Computes a negative log likelihood.
get_NLL <- function(ECC, template_response, human_pc_exp) {
  f <- get_model_results(ECC, template_response, human_pc_exp)
  NLL <- function(a, b, c0, d) {
    if (a <= 0 || b <= 0 || c0 <= 0) {
      return(1e+06)
    }
    
    model_response <- f(a, b, c0, d)
    NLL_val <- sum(model_response$NLL)
    return(NLL_val)
  }
}
