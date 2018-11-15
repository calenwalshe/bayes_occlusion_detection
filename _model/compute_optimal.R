#' Title
#'
#' @return
#' @export
#'
#' @examples
compute_optimal_roots <- function() {
  
  source(
    '~/Dropbox/Calen/Work/occluding/detection_model_analysis/_model/compute_error.R'
  )
  source('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_model/find_root.R')
  source(
    '~/Dropbox/Calen/Work/occluding/detection_model_analysis/_model/find_roots_suboptim.R'
  )
  
  library(dplyr)
  library(mvtnorm)
  library(DEoptim)
  library(tidyr)
  library(purrr)
  library(Rmpfr)
  
  # Error Rate Model Final Code
  load(
    '~/Dropbox/Calen/Work/occluding/detection_model_analysis/_data/model_wide.rdata'
  )
  
  # Compute roots
  
  # Nest data
  model.optimal <- model.wide %>%
    arrange(TARGET, BIN, eccentricity) %>%
    filter(eccentricity > 1) %>%
    group_by(BIN) %>%
    nest()
  
  # Number of rows in dataframe
  n_row <- nrow(model.optimal)
  
  # Use mclapply to compute faster
  roots <- mclapply(
    1:n_row,
    FUN = function(x) {
      data <- model.optimal[x,]$data[[1]]
      
      # start dplyr chain
      model.return <- data %>%
        arrange(TARGET, eccentricity) %>%
        filter(eccentricity > 1) %>%
        rowwise() %>%
        mutate(
          n_obs = length(data_response_vec_0[, 1]),
          SUBJECT     = "model",
          # compute roots
          roots       = list(
            find_root(
              data_mean_vec_0,
              data_mean_vec_1,
              data_cov_mat_0,
              data_cov_mat_1
            )
          ),
          roots_mle   = list(
            find_roots_optim(
              roots[1:3],
              data_mean_vec_0,
              data_mean_vec_1,
              data_cov_mat_0,
              data_cov_mat_1
            )
          ),
          coord       = list(roots_mle[1:3]),
          scale_gauss = list(
            .5 * (coord - data_mean_vec_0) %*% solve(data_cov_mat_0) %*% (coord - data_mean_vec_0)
          )
        ) %>%
        # select variables
        dplyr::select(
          TARGET,
          scale_gauss,
          coord,
          roots,
          roots_mle,
          eccentricity,
          data_mean_vec_0,
          data_mean_vec_1,
          data_cov_mat_0,
          data_cov_mat_1,
          SUBJECT
        )
      # end dplyr chain
      return(model.return)
    },
    mc.cores = 16
  )
  
  model.optimal$roots <- roots
  
  model.optimal <- model.optimal %>% 
    unnest(roots)  
  
  return(model.optimal)
}

#' Error rate for optimal model
#'
#' @param target_mean
#' @param target_cov
#' @param mod_mean1
#' @param mod_mean2
#' @param mod_cov1
#' @param mod_cov2
#' @param start.coord
#' @param box.sd
#' @param box.step
#' @param lower
#'
#' @return
#' @export

compute_optimal_model <- function(optimal.roots, width = 5) {
  
  model.optimal <- optimal.roots
  models.list <- mclapply(
    1:15,
    FUN = function(x) {
      model.optimal.1 <- model.optimal %>%
        filter(BIN == x)
      
      model.optimal.2 <- model.optimal.1 %>%
        ungroup() %>%
        group_by(TARGET, eccentricity) %>%
        nest() 
      
      model.optimal.3 <- model.optimal.2 %>%
        mutate(error = map(data, function(x) {
          data_mean_vec_0 <- x$data_mean_vec_0[[1]]
          data_mean_vec_1 <- x$data_mean_vec_1[[1]]
          data_cov_mat_0 <- x$data_cov_mat_0[[1]]
          data_cov_mat_1 <- x$data_cov_mat_1[[1]]
          coord <- x$coord[[1]]
          
          sz <- width 
          res <- .05
          
          # Error rates. Separately for two categories.
          error.opt.fa = compute_error(
            data_mean_vec_0,
            data_cov_mat_0,
            data_mean_vec_0,
            data_mean_vec_1,
            data_cov_mat_0,
            data_cov_mat_1,
            coord,
            sz,
            res
          )
          error.opt.miss = compute_error(
            data_mean_vec_1,
            data_cov_mat_1,
            data_mean_vec_1,
            data_mean_vec_0,
            data_cov_mat_1,
            data_cov_mat_0,
            coord,
            sz,
            res
          )
          
          list(error.opt.fa = error.opt.fa, error.opt.miss = error.opt.miss)
        }))
      
      model.optimal.2$BIN <- x
      
      model.return <- model.optimal.2
      
      return(model.return)
    } ,
    mc.cores = 16
  )
  
  model.error <- do.call(rbind, models.list)
  
  model.error <- model.error %>%
    mutate(opt.fa = map(error, 1), opt.miss = map(error, 2))
  
  # Convert to equivalent dprime
  model.error$dprime <-
    map2(model.error$opt.fa, model.error$opt.miss, function(x, y) {
      error <-
        1 / 2 * x + 1 / 2 * y
      error.log <-
        as.numeric(log(error))
      result <- -2 * qnorm(error.log, log.p = T)
    })
  
  model.error <- model.error %>%
    unnest(dprime)
  
  model.error$sub_type <- factor("optimal")
  model.error$observer <- factor("optimal")
  
  return(model.error)
}
