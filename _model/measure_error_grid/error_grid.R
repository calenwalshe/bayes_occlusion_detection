g <- function() {
  library(numDeriv)
  library(purrrlyr)
  library(ggplot2)
  sz <- 5
  res <- .05
  
  #tt <- scaled.model %>% ungroup() %>% filter(TARGET == "vertical", BIN == 1, eccentricity %in% unique(eccentricity)[c(1,2,3,4)]) %>%
  model.list <- mclapply(seq(5,6,.02), FUN = function(x) {
    sz <- x
    print(x)
    scaled.model %>% ungroup() %>% filter(TARGET %in% c("vertical", "horizontal", "bowtie", "spot"), eccentricity %in% unique(eccentricity)[c(1,2,3,4)]) %>%
      group_by(TARGET, BIN, eccentricity) %>%
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
        
        error.opt.fa = f(data_mean_vec_0, data_cov_mat_0, data_mean_vec_0, data_mean_vec_1, data_cov_mat_0, data_cov_mat_1, coord, sz, res, lower = F)
        error.opt.miss = f(data_mean_vec_1, data_cov_mat_1, data_mean_vec_1, data_mean_vec_0, data_cov_mat_1, data_cov_mat_0, coord, sz, res, lower = T)
        error.sub.fa   = f(data_mean_vec_0, data_cov_mat_0, mod_mean0, mod_mean1, mod_cov0, mod_cov1, mod_coord_2, sz, res, lower = F)
        error.sub.miss = f(data_mean_vec_1, data_cov_mat_1, mod_mean1, mod_mean0, mod_cov1, mod_cov0, mod_coord_3, sz, res, lower = F)
        
        list(error.opt.fa = error.opt.fa, error.opt.miss = error.opt.miss, error.sub.fa = error.sub.fa, error.sub.miss = error.sub.miss)
      })) %>%
      mutate(sz = x)
    
    }, mc.cores = 16)
  
  tt          <- do.call(rbind, model.list)
  tt <- tt %>% mutate(opt.fa = map(error, 1), opt.miss = map(error,2), sub.fa = map(error,3), sub.miss = map(error,4))
  
  tt$dprime.opt <- map2(tt$opt.fa, tt$opt.miss, function(x,y) {error <- 1/2 * x + 1/2 * y; error.log <- as.numeric(log(error)); result <- -2 * qnorm(error.log, log.p = T)})
  tt$dprime.sub <- map2(tt$sub.fa, tt$sub.miss, function(x,y) {error <- 1/2 * x + 1/2 * y; error.log <- as.numeric(log(error)); result <- -2 * qnorm(error.log, log.p = T)})
  
  tt.collapse <- tt %>% gather(type, dprime, dprime.opt, dprime.sub) %>% unnest(dprime)

  lapply(c("vertical", "horizontal", "bowtie", "spot"), FUN = function(x) {
    fig <- tt.collapse %>% filter(TARGET == x) %>% ggplot(., aes(x = eccentricity, y = dprime, colour = type, shape = as.factor(sz))) + geom_point() + geom_line() + theme(aspect.ratio = 1) + facet_wrap(~BIN, nrow = 4, scale = "free_y") + ggtitle(x)
    ggsave(filename = paste0(x, '-', sz, ".pdf"), path = "~/Dropbox/Calen/Dropbox/", plot = fig, width = 20, height = 20, units = "in")
  })
  
  t_result <- tt %>% unnest(dprime.opt, dprime.sub)
  
  return(t_result)
}


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
  
  #all.bound <- by_row(all.bound, function(x) diag(3) %*% grad(f.grad, as.numeric(x[,1:3])), .collate = "cols")
  all.bound <- all.bound %>% mutate(X = X - target_mean_new[1], Y = Y - target_mean_new[2], Z = Z - target_mean_new[3])
  
  # all.bound.directions <- all.bound %>% rowwise() %>% 
  #   mutate(direction = sign(c(.out1, .out2, .out3)[dim_cum])) %>% 
  #   group_by(direction) %>%
  #   nest()
  
  # results <- by_row(all.bound.directions, function(x) {
  #   lower   <- ifelse(x$direction > 0, T, F)
  #   dim_marg <- setdiff(c(1,2,3), dim_cum)
  #   coord_data <- as.matrix(x$data[[1]][,1:3])
  #   sum(width.frac^2 * (pnorm(coord_data[,dim_cum], 0, 1, lower.tail = lower) * dnorm(coord_data[,dim_marg[1]], 0, 1) * dnorm(coord_data[,dim_marg[2]], 0, 1)))
  # })
  # 
  # error <- sum(unlist(results$.out))

  coord_data <- mpfr(as.matrix(all.bound[, 1:3]), precBits = 120)
  error <- sum(width.frac^2 * (pnorm(coord_data[,dim_cum], 0, 1, lower.tail = lower) * dnorm(coord_data[,dim_marg[1]], 0, 1) * dnorm(coord_data[,dim_marg[2]], 0, 1)))
  
  #all.bound.mat <- mpfr(as.matrix(all.bound[,1:3]), precBits = 120)
  #all.bound.mat.lower <- mpfr(as.matrix(all.bound[all.bound[, 4 + dim_cum] < 0,1:3]), precBits = 120)
  #all.bound.mat.upper <- mpfr(as.matrix(all.bound[all.bound[, 4 + dim_cum] >= 0,1:3]), precBits = 120)
  
  #val.sum.lower <- 0
  #val.sum.upper <- 0
  
  #browser()
  #if(!isempty(all.bound.mat.lower)){
  #  val.sum.lower <- by_row(all.bound.directions, function(x) )
  #}
  
  #if(!isempty(all.bound.mat.upper)) {
  #  val.sum.upper <- sum(width.frac^2 * (pnorm(all.bound.mat.upper[,dim_cum], 0, 1, lower.tail = T) * dnorm(all.bound.mat.upper[,dim_marg[1]], 0, 1) * dnorm(all.bound.mat.upper[,dim_marg[2]], 0, 1)))
  #}
  
  #error <- sum(val.sum.lower) + sum(val.sum.upper)
  #sum(width.frac^2 * (pnorm(all.bound.mat[,dim_cum], 0, 1, lower.tail = lower) * dnorm(all.bound.mat[,dim_marg[1]], 0, 1) * dnorm(all.bound.mat[,dim_marg[2]], 0, 1)))

}





