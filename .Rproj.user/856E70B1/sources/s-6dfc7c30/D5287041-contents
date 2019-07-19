## Compute Suboptimal Models
## Return dataframe with psychometric parameters

library(dplyr)
library(Rmpfr)
library(purrrlyr)
library(purrr)
library(tidyr)
library(parallel)

source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/polar_error.R')
source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/polar_roots.R')
source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/regular_placement.R')
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model_wide.rdata")

n.rows       <- nrow(model.wide)
n.resolution <- 2^15
n.cores      <- 16
bPlot        <- 0

model.wide.nocov <- model.wide
model.wide.nocov$mod_cov0 <- map(model.wide$data_cov_mat_0, function(x) x * diag(3))
model.wide.nocov$mod_cov1 <- map(model.wide$data_cov_mat_1, function(x) x * diag(3))
model.wide.nocov$mod_mean0 <- model.wide.nocov$data_mean_vec_0
model.wide.nocov$mod_mean1 <- model.wide.nocov$data_mean_vec_1
model.wide.nocov$model.type <- "nocov"

model.wide.stdabs <- model.wide
model.wide.stdabs$mod_cov0 <- map2(model.wide$data_cov_mat_0,model.wide$data_cov_mat_1, function(x,y) (x * diag(3) + y * diag(3))/2)
model.wide.stdabs$mod_cov1 <- model.wide.stdabs$mod_cov0
model.wide.stdabs$mod_mean0 <- model.wide.nocov$data_mean_vec_0
model.wide.stdabs$mod_mean1 <- model.wide.nocov$data_mean_vec_1
model.wide.stdabs$model.type <- "stdabs"

model.wide.stdave <- model.wide
model.wide.stdave$mod_cov0 <- map(model.wide$data_cov_mat_0, function(x) x * diag(3))
model.wide.stdave$mod_cov1 <- map(model.wide$data_cov_mat_0, function(x) x * diag(3))
model.wide.stdave$mod_mean0 <- model.wide.nocov$data_mean_vec_0
model.wide.stdave$mod_mean1 <- model.wide.nocov$data_mean_vec_1
model.wide.stdave$model.type <- "stdave"

model.wide.optimal <- model.wide
model.wide.optimal$mod_cov0 <- model.wide$data_cov_mat_0
model.wide.optimal$mod_cov1 <- model.wide$data_cov_mat_1
model.wide.optimal$mod_mean0 <- model.wide$data_mean_vec_0
model.wide.optimal$mod_mean1 <- model.wide$data_mean_vec_1
model.wide.optimal$model.type <- "optimal"

valid.dim <- list(c(1,0,0), c(0,1,0), c(0,0,1))

model.wide$model.type <- list("optimal")
model.wide.nocov$model.type <- list("nocov")

combination <- lapply(valid.dim, FUN = function(x) {
  decode_dim      <- x
  model.wide.diag <- model.wide %>%
    rowwise() %>% 
    mutate(mod_mean0 = list(data_mean_vec_0 * decode_dim), mod_mean1 = list(data_mean_vec_1 * decode_dim), mod_cov0 = list(diag(data_cov_mat_0[diag(3)==T] * decode_dim + !decode_dim)), mod_cov1 = list(diag(data_cov_mat_1[diag(3)==T] * decode_dim + !decode_dim))) %>%
    mutate(model.type = list(decode_dim))
})


models.list <- c(list(model.wide.optimal, model.wide.nocov), combination)

n.resolution <- 2^17
vals <- regular.placement(1, n.resolution)

model.measure <- lapply(models.list, FUN = function(model) {
  n.rows      <- nrow(model)
  model.wide <- model
  
  errors <- mclapply(1:n.rows, FUN = function(x) {
    print(x)
    target_cov_0  <- model.wide$data_cov_mat_0[[x]]
    target_mean_0 <- model.wide$data_mean_vec_0[[x]]
    
    target_cov_1  <- model.wide$data_cov_mat_1[[x]]
    target_mean_1 <- model.wide$data_mean_vec_1[[x]]
    
    mod_cov1 <- model.wide$mod_cov0[[x]]
    mod_cov2 <- model.wide$mod_cov1[[x]]
    
    mod_mean1 <- model.wide$mod_mean0[[x]]
    mod_mean2 <- model.wide$mod_mean1[[x]]
    
    return.val.1 <- polar.error(mod_mean1, mod_mean2, mod_cov1, mod_cov2, target_mean_0, target_cov_0, vals = vals, pr_a = .5)
    return.val.2 <- polar.error(mod_mean2, mod_mean1, mod_cov2, mod_cov1, target_mean_1, target_cov_1, vals = vals, pr_a = .5)
    
    responses <- list(falsealarm = return.val.1, miss = return.val.2)
  }, mc.cores = n.cores)
  
  model.wide$responses  <- errors
  
  return(model.wide)
})

model.all   <- do.call(rbind, model.measure)

# Compute dprime from percent correct
dprime.vals <- map(model.all$responses, function(x) {
  miss        <- x$miss
  falsealarm  <- x$falsealarm
  dprime      <- -2*qnorm(as.numeric(log((x$miss + x$falsealarm)/2)), log.p = T)
})

model.all$dprime  <- map(dprime.vals, 1)

model.all         <- model.all %>% 
  unnest(dprime)

model.all$sub_type <- model.all$model.type
model.all$sub_type <- map(model.all$sub_type, function(x) paste(x, collapse = ' '))
model.all$observer <- "model"

model.all         <- model.all %>% unnest(sub_type)
model.all$SUBJECT <- model.all$sub_type

model.error <- model.all

model.error$observer <- model.error$sub_type

save(file = '~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model_error.rdata', model.error)

# Plot and check optimization
if (bPlot == 1) {
  model.psychometric <- get_model_psychometric(model.error, .07)
  model.psychometric$SUBJECT <- model.psychometric$sub_type
  model.psychometric$se <- 0
  
  load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/human.psychometrics.rdata")
  
  human.psychometrics$se <- 0
  
  plot_psychometric(human.psychometrics)
  plot.model.psychometric(model.psychometric)
  
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R')
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/plot_theme.R')
  model.threshold <- get_threshold(model.psychometric)
  model.threshold$SUBJECT <- model.threshold$sub_type
  human.threshold  <- get_threshold(human.psychometrics) %>% mutate(sub_type = "human")
  
  plot.vals <- plot_publication_thresholds(human.thresholds = human.threshold, model.thresholds = model.threshold, statIn = "Lvals")
  
  fig.1 <- plot.vals + scale_linetype_manual("", labels = c("rcw", "sps", "yhb","Optimal", "No Covariance", "Mahalanobis"), values = c(1,1,1,2,2,2)) +
    scale_colour_manual("", labels = c( "rcw", "sps", "yhb","Optimal", "No Covariance", "Mahalanobis"), values=colour.vals)
  plot(fig.1)  
}



