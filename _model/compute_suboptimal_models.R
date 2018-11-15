#' Computes the error rate for a suboptimal model.
#' Suboptimal models have no covariance between cues and have information removed for one or two dimensions. 
#'
#' @param valid_dim 
#'
#' @return
#' @export
#'
#' @examples
compute_suboptimal <- function(valid_dim) {
  library(DEoptim)
  
  source('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_model/find_root.R')
  source('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_model/find_roots_suboptim.R')
  
  # Error Rate Model Final Code
  load('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_data/model_wide.rdata')
  
  decode_dim <- valid_dim
  
  model.wide.diag <- model.wide %>%
    rowwise() %>% 
    filter(eccentricity > 1) %>%
    mutate(mod_mean0 = list(data_mean_vec_0 * decode_dim), mod_mean1 = list(data_mean_vec_1 * decode_dim), mod_cov0 = list(diag(data_cov_mat_0[diag(3)==T] * decode_dim + !decode_dim)), mod_cov1 = list(diag(data_cov_mat_1[diag(3)==T] * decode_dim + !decode_dim)))
  
  model.wide.diag <- model.wide.diag %>% 
    group_by(TARGET, BIN, eccentricity) %>%
    nest()
  
  roots.solve <- mclapply(model.wide.diag$data, FUN = function(x) {
    mod_mean0 <- x$mod_mean0[[1]]
    mod_mean1 <- x$mod_mean1[[1]]
    mod_cov0 <- x$mod_cov0[[1]]
    mod_cov1 <- x$mod_cov1[[1]]
    
    data_mean_vec_0 <- x$data_mean_vec_0[[1]]
    data_mean_vec_1 <- x$data_mean_vec_1[[1]]
    data_cov_mat_0 <- x$data_cov_mat_0[[1]]
    data_cov_mat_1 <- x$data_cov_mat_1[[1]]
    
    mod_roots       = find_root_ordered(mod_mean0, mod_mean1, mod_cov0, mod_cov1)
    mod_roots_mle   = find_roots_optim_ordered(mod_roots[1:3], mod_mean0, mod_mean1, mod_cov0, mod_cov1)
    mod_roots_mle_2 = find_roots_suboptim_ordered(mod_roots[1:3], mod_mean0, mod_mean1, data_mean_vec_0, mod_cov0, mod_cov1, data_cov_mat_0)
    mod_roots_mle_3 = find_roots_suboptim_ordered(mod_roots[1:3], mod_mean0, mod_mean1, data_mean_vec_1, mod_cov0, mod_cov1, data_cov_mat_1)
    mod_coord       = mod_roots_mle[1:3]
    mod_coord_2     = mod_roots_mle_2[1:3]
    mod_coord_3     = mod_roots_mle_3[1:3]
    mod_scale_gauss = .5*(mod_coord - mod_mean0) %*% solve(mod_cov0) %*% (mod_coord - mod_mean0)
    mod_scale_gauss_2 = .5*(mod_coord_2 - data_mean_vec_0) %*% solve(data_cov_mat_0) %*% (mod_coord_2 - data_mean_vec_0)
    mod_scale_gauss_3 = .5*(mod_coord_3 - data_mean_vec_1) %*% solve(data_cov_mat_1) %*% (mod_coord_3 - data_mean_vec_1)
        
    return.frame <- list(mod_roots = (mod_roots), mod_roots_mle = (mod_roots_mle), mod_roots_mle_2 = (mod_roots_mle_2), 
                               mod_roots_mle_3 = (mod_roots_mle_3), mod_coord = (mod_coord), mod_coord_2 = (mod_coord_2), 
                               mod_coord_3 = (mod_coord_3), mod_scale_gauss = (mod_scale_gauss), 
                               mod_scale_gauss_2 = (mod_scale_gauss_2), mod_scale_gauss_3 = (mod_scale_gauss_3))
    

    return(return.frame)
  }, mc.cores = 16)
  
  model.wide.diag$roots.solve <- roots.solve
  
  n_row <- nrow(model.wide.diag)
  
  error <- mclapply(1:n_row, FUN = function(x) {
    row <- model.wide.diag[x, ]
    
    data <- row$data[[1]]
    roots.solve <-  row$roots.solve[[1]]
    
    
    
    data_mean_vec_0 <- data$data_mean_vec_0[[1]]
    data_mean_vec_1 <- data$data_mean_vec_1[[1]]
    data_cov_mat_0 <- data$data_cov_mat_0[[1]]
    data_cov_mat_1 <- data$data_cov_mat_1[[1]]
    
    mod_mean0 <- data$mod_mean0[[1]]
    mod_mean1 <- data$mod_mean1[[1]]
    
    mod_cov0 <- data$mod_cov0[[1]]
    mod_cov1 <- data$mod_cov1[[1]]
    
    mod_coord_2 <- roots.solve$mod_coord_2
    mod_coord_3 <- roots.solve$mod_coord_3
    
    sz <- 5
    res <- .05
    
    error.sub.fa   = compute_error(data_mean_vec_0, data_cov_mat_0, mod_mean0, mod_mean1, mod_cov0, mod_cov1, mod_coord_2, sz, res)
    error.sub.miss = compute_error(data_mean_vec_1, data_cov_mat_1, mod_mean1, mod_mean0, mod_cov1, mod_cov0, mod_coord_3, sz, res)
    
    list(error.sub.fa = error.sub.fa, error.sub.miss = error.sub.miss)
  }, mc.cores = 16)
  
  model.wide.diag$error <- error
  
  dprime <- map(model.wide.diag$error, function(x) {
    error <- 1/2 * x$error.sub.fa + 1/2 * x$error.sub.miss; error.log <- as.numeric(log(error)); result <- -2 * qnorm(error.log, log.p = T)
  })
  
  model.wide.diag$dprime <- dprime
  
  model.wide.diag <- model.wide.diag %>% unnest(dprime)
  
  model.wide.diag$sub_type <- list(valid_dim)
  
  model.wide.return <- model.wide.diag %>% select(TARGET, BIN, eccentricity, dprime, sub_type) %>% arrange(TARGET, BIN, eccentricity)
  return(model.wide.return)
}
