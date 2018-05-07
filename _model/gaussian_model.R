# Compute performance of the Bayes' optimal detector
get_optimal_model <- function(model.wide) {
  library(dplyr)
  library(broom)
  library(bbmle)
  library(R2Cuba)
  library(multidplyr)
  library(mvtnorm)
  library(purrr)
  library(tidyr)
  library(bbmle)
  library(DEoptim)
  
  integrate_distribution <- function(scale = 1, mean1, mean2, cov1, cov2, int.region = c(0,0,0), box.width = 1) {
    inv_cov1 <- solve(cov1)
    inv_cov2 <- solve(cov2)
    
    scalar.dist1 <- 1 / sqrt(det(2 * pi * cov1))
    scalar.dist2 <- 1 / sqrt(det(2 * pi * cov2))
    
    f <- function(x) {
      gauss_f1 <- scalar.dist1 * exp(-.5*(t(x - mean1) %*% inv_cov1 %*% (x - mean1) - scale))
      gauss_f2 <- scalar.dist2 * exp(-.5*(t(x - mean2) %*% inv_cov2 %*% (x - mean2) - scale))
      min(gauss_f1,gauss_f2, na.rm = T)
    }
    
    box_max <- apply(rbind(cov1[diag(3) == 1], cov2[diag(3) == 1]), 2, max)
    
    results <- cuhre(3, 1, f, 
                     lower = int.region - sqrt(box_max),
                     upper = int.region + sqrt(box_max),
                     flags = list(verbose = 0), rel.tol = 10e-4, max.eval = 10e7)
    
    return(results)
  }
    
  cluster <- create_cluster(cores = 16)
  cluster_assign_value(cluster, 'integrate_distribution', integrate_distribution)
  cluster_library(cluster, 'bbmle')
  cluster_library(cluster, 'purrr')
  cluster_library(cluster, 'R2Cuba')
  cluster_library(cluster, 'bbmle')
  cluster_library(cluster, 'tidyr')
  cluster_library(cluster, 'dplyr')
  cluster_library(cluster, 'broom')
  cluster_library(cluster, 'mvtnorm')
  set_default_cluster(cluster)
  
  #box.width.df <- data.frame(box.width = (seq(.01,3,length.out = 20)))

  scaled.model <- model.wide %>% 
    arrange(TARGET, BIN, eccentricity) %>%
    rowwise() %>%
    mutate(n_obs = length(data_response_vec_0[,1]),
           SUBJECT     = "model",
           roots       = list(find_root(data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1)),
           roots_mle   = list(find_roots_optim(roots[1:3], data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1)),
           coord       = list(roots_mle[1:3]),
           scale_gauss = list((coord - data_mean_vec_0) %*% solve(data_cov_mat_0) %*% (coord - data_mean_vec_0))) %>%
    select(TARGET, BIN, scale_gauss, coord, roots, roots_mle, eccentricity, data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1, SUBJECT) 

    scaled.integrate <- scaled.model %>%
    group_by(TARGET, BIN, eccentricity) %>%
    as_tibble() %>%
    nest() %>%
    partition(TARGET, BIN, eccentricity, cluster = cluster) %>%
    mutate(integrate_min = 
             map(data, function(x) 
               {
               integrate_distribution(x$scale_gauss[[1]], x$data_mean_vec_0[[1]], x$data_mean_vec_1[[1]], x$data_cov_mat_0[[1]], x$data_cov_mat_1[[1]], x$coord[[1]])
               }
             )) %>%
      collect() %>%
      unnest(data) %>%
      mutate(integral = map(integrate_min, function(x) { x$value }))
    
    integrate.result <- scaled.integrate %>% mutate(dprime = -2*qnorm(log(integral[[1]]) - .5*scale_gauss[[1]], 0, 1, log.p = T),
           dprime_compare = -2*qnorm(-.5*scale_gauss[[1]], 0, 1, log.p = T)) %>%
    arrange(TARGET, BIN, eccentricity)

  return(integrate.result)
}

# Compute the model based on dprime summation with assumptions of independence and equal variance.
get_independent_model <- function(model.wide) {
  optimal.df <- model.wide %>%
    select(BIN, TARGET, eccentricity, data_mean_vec_0, data_mean_vec_1, data_sd_vec_0, data_sd_vec_1)
  
  optimal.dprime <- by_row(optimal.df, function(row) {
    dprime <- sqrt(
      sum(
        (
          abs(row$data_mean_vec_1[[1]] - row$data_mean_vec_0[[1]]) /
            sqrt(row$data_sd_vec_0[[1]]^2 + row$data_sd_vec_1[[1]]^2)
          )^2
        )) # sqrt sum dprime^2
  }, .collate = "row", .to = "dprime")
  
  optimal.dprime$SUBJECT <- "independent_cue"
  
  return(optimal.dprime)
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

#' Compute the efficiency of the optimal and human responses
get_efficiency <- function(human.psychometric, optimal.observer) {
  library(purrrlyr)
  
  model.dat <- optimal.observer %>%
    select(SUBJECT, BIN, TARGET, eccentricity, dprime)
  
  human.dat <- human.psychometric %>%
    select(SUBJECT, BIN, TARGET, d0, e0, b, gamma) %>%
    mutate(e0 = as.numeric(e0))
  
  optim.human.dat <- merge(model.dat, human.dat, by = c("TARGET", "BIN"))
  
  efficiency <- by_row(optim.human.dat, function(row) {
    e0  <- row$e0
    b   <- row$b
    d0  <- row$d0
    ecc <- row$eccentricity
    
    dprime.human <- d0 * e0^b/(e0^b + ecc^b)
  }, .to = "dprime.human") %>%
    unnest(dprime.human) %>%
    mutate(efficiency = dprime.human/dprime) %>%
    select(BIN, TARGET, eccentricity, efficiency) %>%
    arrange(BIN, TARGET, eccentricity)
  
  return(efficiency)
}
