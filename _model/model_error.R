#' Error rate for
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
f <- function(target_mean, target_cov, mod_mean1, mod_mean2, mod_cov1, mod_cov2, start.coord, box.sd, box.step, lower) {
  width.frac <- box.step
  arr.sz     <- box.sd
  
  L <- solve(t(chol(target_cov)))
  
  target_cov_delt <- L %*% target_cov %*% t(L)
  mod_cov0_delt <- L %*% mod_cov1 %*% t(L)
  mod_cov1_delt <- L %*% mod_cov2 %*% t(L)
  
  x_new <- t(L %*% start.coord)
  mod_mean0_new <- t(L %*% mod_mean1)
  mod_mean1_new <- t(L %*% mod_mean2)
  target_mean_new <- t(L %*% target_mean)
  
  f.grad <- function(x) {
    val <- dmvnorm(c(x[1],x[2],x[3]), mod_mean0_new, mod_cov0_delt, log = T) - dmvnorm(c(x[1],x[2],x[3]), mod_mean1_new, mod_cov1_delt, log = T)  
  }
  
  grad.val <- grad(f.grad, x_new)
  
  dim_cum.old <- which.max(abs(x_new - target_mean_new))
  
  dim_cum     <- which.max(abs(grad.val))
  dim_marg <- setdiff(c(1,2,3), dim_cum)
  
  lower <- ifelse(sign(grad.val[dim_cum]) > 0, T, F)
  
  
  W.1 <- seq(x_new[dim_cum] - arr.sz, x_new[dim_cum] + arr.sz, width.frac)
  W.2 <- expand.grid(Y = seq(x_new[dim_marg[1]] - arr.sz, x_new[dim_marg[1]] + arr.sz, width.frac), Z = seq(x_new[dim_marg[2]] - arr.sz, x_new[dim_marg[2]] + arr.sz, width.frac))
  
  W   <- lapply(W.1, FUN = function(x) tibble(X = x, data = list(data.frame(W.2)))) %>% do.call(rbind, .)
  
  discrim.eval   <- map2(W$X, W$data, function(x, y) {
    coord <- as.matrix(cbind(x, y))
    g <- function(x) {
      coord.new <- matrix(0, nrow(x), 3)
      coord.new[,dim_cum] <- x[,1]
      coord.new[,dim_marg[1]] <- x[,2]
      coord.new[,dim_marg[2]] <- x[,3]
      val <- dmvnorm(coord.new, mod_mean0_new, mod_cov0_delt, log = T) - dmvnorm(coord.new, mod_mean1_new, mod_cov1_delt, log = T)      
    }
    val <- g(coord)
    cbind(y, val)
  }
  )
  
  W$data <- discrim.eval
  
  bound.coords <- map2(W$data[1:length(W$data)-1], W$data[2:length(W$data)], function(x, y) {
    x[x[,3] * y[,3] < 0,]
  })
  
  W[1:nrow(W)-1,]$data <- bound.coords
  W <- W[1:nrow(W)-1,]
  
  all.bound <- W %>%
    rowwise() %>%
    filter(!is.null(W)) %>%
    unnest(data)
  
  all.bound[,c(dim_cum, dim_marg[1], dim_marg[2])] <- all.bound[,1:3] # replace order
  all.bound <- all.bound %>% mutate(X = X - target_mean_new[1], Y = Y - target_mean_new[2], Z = Z - target_mean_new[3])
  
  coord_data <- mpfr(as.matrix(all.bound[, 1:3]), precBits = 120)
  error <- sum(width.frac^2 * (pnorm(coord_data[,dim_cum], 0, 1, lower.tail = lower) * dnorm(coord_data[,dim_marg[1]], 0, 1) * dnorm(coord_data[,dim_marg[2]], 0, 1)))
}

# Error Rate Model Final Code
load('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_data/model_wide.rdata')

dims <- list(c(1,0,0), c(0,1,0), c(0,0,1), c(1,1,0), c(1,0,1), c(0,1,1), c(1,1,1))

all.models <- lapply(dims, FUN = function(x) {
  decode_dim <- x
  
  model.wide.diag <- model.wide %>%
    rowwise() %>% 
    filter(eccentricity > 1) %>%
    mutate(mod_mean0 = list(data_mean_vec_0 * decode_dim), mod_mean1 = list(data_mean_vec_1 * decode_dim), mod_cov0 = list(diag(data_cov_mat_0[diag(3)==T] * decode_dim + !decode_dim)), mod_cov1 = list(diag(data_cov_mat_1[diag(3)==T] * decode_dim + !decode_dim)))
  
  model.wide.diag <- model.wide.diag %>% 
    arrange(TARGET, BIN, eccentricity) %>%
    rowwise() %>%
    mutate(n_obs = length(data_response_vec_0[,1]),
           SUBJECT     = "model",
           roots       = list(find_root(data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1)),
           roots_mle   = list(find_roots_optim(roots[1:3], data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1)),
           coord       = list(roots_mle[1:3]),
           scale_gauss = list(.5*(coord - data_mean_vec_0) %*% solve(data_cov_mat_0) %*% (coord - data_mean_vec_0)),
           mod_roots       = list(find_root_ordered(mod_mean0, mod_mean1, mod_cov0, mod_cov1)),
           mod_roots_mle   = list(find_roots_optim_ordered(mod_roots[1:3], mod_mean0, mod_mean1, mod_cov0, mod_cov1)),
           mod_roots_mle_2 = list(find_roots_suboptim_ordered(mod_roots[1:3], mod_mean0, mod_mean1, data_mean_vec_0, mod_cov0, mod_cov1, data_cov_mat_0)),
           mod_roots_mle_3 = list(find_roots_suboptim_ordered(mod_roots[1:3], mod_mean0, mod_mean1, data_mean_vec_1, mod_cov0, mod_cov1, data_cov_mat_1)),
           mod_coord       = list(mod_roots_mle[1:3]),
           mod_coord_2     = list(mod_roots_mle_2[1:3]),
           mod_coord_3     = list(mod_roots_mle_3[1:3]),
           mod_scale_gauss = list(.5*(mod_coord - mod_mean0) %*% solve(mod_cov0) %*% (mod_coord - mod_mean0)),
           mod_scale_gauss_2 = list(.5*(mod_coord_2 - data_mean_vec_0) %*% solve(data_cov_mat_0) %*% (mod_coord_2 - data_mean_vec_0)),
           mod_scale_gauss_3 = list(.5*(mod_coord_3 - data_mean_vec_1) %*% solve(data_cov_mat_1) %*% (mod_coord_3 - data_mean_vec_1))) %>% 
    dplyr::select(TARGET, BIN, mod_scale_gauss, scale_gauss, coord, mod_coord, roots, roots_mle, mod_scale_gauss_2, mod_coord_2, mod_scale_gauss_3, mod_coord_3, mod_roots_mle_2, eccentricity, data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1, mod_mean0, mod_mean1, mod_cov0, mod_cov1, SUBJECT)
  
  
  models.list <- mclapply(1:15, FUN = function(x) {
    model.wide.temp <- model.wide.diag %>% 
      filter(BIN == x)
    
    model.wide.temp <- model.wide.temp %>%
      ungroup() %>%
      group_by(TARGET, eccentricity) %>%
      nest() %>%
      mutate(error = map(data, function(x) {
        data_mean_vec_0 <- x$data_mean_vec_0[[1]]
        data_mean_vec_1 <- x$data_mean_vec_1[[1]]
        data_cov_mat_0 <- x$data_cov_mat_0[[1]]
        data_cov_mat_1 <- x$data_cov_mat_1[[1]]
        
        mod_mean0 <- x$mod_mean0[[1]]
        mod_mean1 <- x$mod_mean1[[1]]
        
        mod_cov0 <- x$mod_cov0[[1]]
        mod_cov1 <- x$mod_cov1[[1]]
        
        coord <- x$coord[[1]]
        mod_coord_2 <- x$mod_coord_2[[1]]
        mod_coord_3 <- x$mod_coord_3[[1]]
        
        sz <- 5
        res <- .05
        error.opt.fa = f(data_mean_vec_0, data_cov_mat_0, data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1, coord, sz, res)
        error.opt.miss = f(data_mean_vec_1, data_cov_mat_1, data_mean_vec_1, data_mean_vec_0, data_cov_mat_1, data_cov_mat_0, coord, sz, res)
        error.sub.fa   = f(data_mean_vec_0, data_cov_mat_0, mod_mean0, mod_mean1, mod_cov0, mod_cov1, mod_coord_2, sz, res)
        error.sub.miss = f(data_mean_vec_1, data_cov_mat_1, mod_mean1, mod_mean0, mod_cov1, mod_cov0, mod_coord_3, sz, res)
        
        list(error.opt.fa = error.opt.fa, error.opt.miss = error.opt.miss, error.sub.fa = error.sub.fa, error.sub.miss = error.sub.miss)
      }))
    
    model.wide.temp$BIN <- x
    
    model.return <- model.wide.temp %>% 
      select(-data)
    
    return(model.wide.temp)
  },mc.cores = 8)
  
  model.error <- do.call(rbind, models.list)
  
  model.error <- model.error %>% mutate(opt.fa = map(error, 1), opt.miss = map(error,2), sub.fa = map(error,3), sub.miss = map(error,4))
  
  model.error$dprime.opt <- map2(model.error$opt.fa, model.error$opt.miss, function(x,y) {error <- 1/2 * x + 1/2 * y; error.log <- as.numeric(log(error)); result <- -2 * qnorm(error.log, log.p = T)})
  model.error$dprime.sub <- map2(model.error$sub.fa, model.error$sub.miss, function(x,y) {error <- 1/2 * x + 1/2 * y; error.log <- as.numeric(log(error)); result <- -2 * qnorm(error.log, log.p = T)})
  
  model.error <- model.error %>% 
    gather(type, dprime, dprime.opt, dprime.sub) %>%
    unnest(dprime)
  
  model.error$sub_type <- list(decode_dim)
  
  return(model.error)
  
  #lapply(c("vertical", "horizontal", "bowtie", "spot"), FUN = function(x) {
  #  fig <- model.error %>% filter(TARGET == x) %>% ggplot(., aes(x = eccentricity, y = dprime, colour = type)) + geom_point() + geom_line() + theme(aspect.ratio = 1) + facet_wrap(~BIN, nrow = 4, scale = "free_y") + ggtitle(x)
  #  ggsave(filename = paste0(x, '-', sz, ".pdf"), path = "~/Dropbox/Calen/Dropbox/", plot = fig, width = 20, height = 20, units = "in")
  #})
})

marginal.models$sub_type <- map(marginal.models$sub_type, function(x) {(paste0(x, collapse = ''))})
marginal.models          <- unnest(marginal.models, sub_type) 
marginal.models$sub_type <- factor(marginal.models$sub_type)

lapply(c("vertical", "horizontal", "bowtie", "spot"), FUN = function(x) {
  fig <- marginal.models %>% filter(TARGET == x) %>% ggplot(., aes(x = eccentricity, y = dprime, colour = type, linetype = sub_type)) + geom_point() + geom_line() + theme(aspect.ratio = 1) + facet_wrap(~BIN, nrow = 4, scale = "free_y") + ggtitle(x)
  ggsave(filename = paste0(x, '-', sz, ".pdf"), path = "~/Dropbox/Calen/Dropbox/", plot = fig, width = 20, height = 20, units = "in")
})

