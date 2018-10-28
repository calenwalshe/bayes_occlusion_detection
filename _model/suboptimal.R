# Compute performance of the Bayes' optimal detector
get_optimal_model_suboptimal <- function(model.wide) {
  library(dplyr)
  library(Rmpfr)
  library(broom)
  library(bbmle)
  library(R2Cuba)
  library(multidplyr)
  library(mvtnorm)
  library(purrr)
  library(tidyr)
  library(bbmle)
  library(DEoptim)
  library(parallel)
  
  integrate_distribution <- function(scale = 1, mod_mean1, mod_mean2, mod_cov1, mod_cov2, target_mean, target_cov, int.region = c(0,0,0), box.width = 1) {
    
    inv_cov1 <- solve(target_cov)
    
    mod_inv_cov1 <- solve(mod_cov1)
    mod_inv_cov2 <- solve(mod_cov2)
    
    scalar.dist1 <- 1 / sqrt(det(2 * pi * target_cov))
    
    mod_scalar.dist1 <- 1 / sqrt(det(2 * pi * mod_cov1))
    mod_scalar.dist2 <- 1 / sqrt(det(2 * pi * mod_cov2))
    
    f <- function(x) {
      x_t <- x
      
      mod_gauss_f1 <- log(mod_scalar.dist1) + (-.5*(t(x_t - mod_mean1) %*% mod_inv_cov1 %*% (x_t - mod_mean1)))
      mod_gauss_f2 <- log(mod_scalar.dist2) + (-.5*(t(x_t - mod_mean2) %*% mod_inv_cov2 %*% (x_t - mod_mean2)))
      
      if(mod_gauss_f1 < mod_gauss_f2) {
        gauss_f1 <- scalar.dist1 * exp(-.5*(t(x - target_mean) %*% solve(target_cov) %*% (x - target_mean)) + scale)
      }else{
        0
      }
    }
    
    box.int <- box.width * eigen(target_cov)$values
    box.int <- box.width * c(1,1,1)
    
    results <- cuhre(3, 1, f, 
                     lower = int.region - box.int,
                     upper = int.region + box.int,
                     flags = list(verbose = 1), rel.tol = 10e-3)
    
    return(results)
  }
  
  cluster <- create_cluster(cores = 16)
  cluster_assign_value(cluster, 'integrate_distribution', integrate_distribution)
  cluster_assign_value(cluster, 'find_root', find_root)
  cluster_assign_value(cluster, 'find_roots_optim', find_roots_optim)
  cluster_assign_value(cluster, 'find_roots_suboptim', find_roots_suboptim)
  cluster_library(cluster, 'bbmle')
  cluster_library(cluster, 'purrr')
  cluster_library(cluster, 'R2Cuba')
  cluster_library(cluster, 'bbmle')
  cluster_library(cluster, 'tidyr')
  cluster_library(cluster, 'dplyr')
  cluster_library(cluster, 'broom')
  cluster_library(cluster, 'mvtnorm')
  cluster_library(cluster, 'Rmpfr')

  set_default_cluster(cluster)
  
  scaled.model <- model.wide %>% 
    arrange(TARGET, BIN, eccentricity) %>%
    rowwise() %>%
    mutate(n_obs = length(data_response_vec_0[,1]),
           SUBJECT     = "model",
           roots       = list(find_root(data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1)),
           roots_mle   = list(find_roots_optim(roots[1:3], data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1)),
           coord       = list(roots_mle[1:3]),
           scale_gauss = list(.5*(coord - data_mean_vec_0) %*% solve(data_cov_mat_0) %*% (coord - data_mean_vec_0)),
           mod_roots       = list(find_root(mod_mean0, mod_mean1, mod_cov0, mod_cov1)),
           mod_roots_mle   = list(find_roots_optim(mod_roots[1:3], mod_mean0, mod_mean1, mod_cov0, mod_cov1)),
           mod_roots_mle_2 = list(find_roots_suboptim(mod_roots[1:3], mod_mean0, mod_mean1, data_mean_vec_0, mod_cov0, mod_cov1, data_cov_mat_0)),
           mod_roots_mle_3 = list(find_roots_suboptim(mod_roots[1:3], mod_mean0, mod_mean1, data_mean_vec_1, mod_cov0, mod_cov1, data_cov_mat_1)),
           mod_coord       = list(mod_roots_mle[1:3]),
           mod_coord_2     = list(mod_roots_mle_2[1:3]),
           mod_coord_3     = list(mod_roots_mle_3[1:3]),
           mod_scale_gauss = list(.5*(mod_coord - mod_mean0) %*% solve(mod_cov0) %*% (mod_coord - mod_mean0)),
           mod_scale_gauss_2 = list(.5*(mod_coord_2 - data_mean_vec_0) %*% solve(data_cov_mat_0) %*% (mod_coord_2 - data_mean_vec_0)),
           mod_scale_gauss_3 = list(.5*(mod_coord_3 - data_mean_vec_1) %*% solve(data_cov_mat_1) %*% (mod_coord_3 - data_mean_vec_1))) %>% 
    dplyr::select(TARGET, BIN, mod_scale_gauss, scale_gauss, coord, mod_coord, roots, roots_mle, mod_scale_gauss_2, mod_coord_2, mod_scale_gauss_3, mod_coord_3, mod_roots_mle_2, eccentricity, data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1, mod_mean0, mod_mean1, mod_cov0, mod_cov1, SUBJECT)
  
  scaled.integrate <- scaled.model %>%
    filter(BIN == 3) %>%
    group_by(TARGET, BIN, eccentricity) %>%
    as_tibble() %>%
    nest() %>%
    partition(TARGET, BIN, eccentricity, cluster = cluster) %>%
    mutate(integrate_min = 
             map(data, function(x) 
             {
               mod_scale_gauss <- x$mod_scale_gauss[[1]]
               mod_mean0       <- x$mod_mean0[[1]]
               mod_mean1       <- x$mod_mean1[[1]]
               mod_cov0        <- x$mod_cov0[[1]]
               mod_cov1        <- x$mod_cov1[[1]]
               mod_coord       <- x$mod_coord_2[[1]]
                 
               scale_0     <- x$mod_scale_gauss_2[[1]]
               scale_1     <- x$mod_scale_gauss_3[[1]]
               mean0       <- x$data_mean_vec_0[[1]]
               mean1       <- x$data_mean_vec_1[[1]]
               cov0        <- x$data_cov_mat_0[[1]]
               cov1        <- x$data_cov_mat_1[[1]]
               coord_0     <- x$mod_coord_2[[1]]
               coord_1     <- x$mod_coord_3[[1]]
               
               false_alarm <- integrate_distribution(scale = scale_0, mod_mean0, mod_mean1, mod_cov0, mod_cov1, mean0, cov0, coord_0, box.width = 5)
               miss        <- integrate_distribution(scale = scale_1, mod_mean1,mod_mean0, mod_cov1, mod_cov0, mean1, cov1, coord_1, box.width = 5)
               
               data.frame(false_alarm = false_alarm$value, miss = miss$value)
             }
             )) %>%
    collect() %>%
    unnest(integrate_min, data)

  integrate.result <- scaled.integrate %>% 
    rowwise() %>%
    mutate(error = mpfr(miss, precBits = 120)/exp(mpfr(mod_scale_gauss_3[[1]], precBits = 120)) + 
             mpfr(false_alarm, precBits = 120)/exp(mpfr(mod_scale_gauss_2[[1]], precBits = 120)))
  
  integrate.result$dprime <- -2*qnorm(as.numeric(log(integrate.result$error)), log.p = T)

  integrate.result %>% dplyr::select(BIN, TARGET, SUBJECT, eccentricity, dprime) %>% arrange(TARGET, eccentricity)
  
  integrate.result <- scaled.integrate %>% 
    rowwise() %>% 
    mutate(scale_used = min(c(mod_scale_gauss_2[[1]], mod_scale_gauss_3[[1]]))) %>%
    mutate(dprime = abs(2*qnorm(log(integral[[1]]) - scale_used[[1]], 0, 1, log.p = T))) %>%
    arrange(TARGET, BIN, eccentricity)
  
  return(integrate.result)
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

