#' Compute an error rate for two gaussian distributions.
#' The discriminant bound can be suboptimal. This is computed mean mod != target
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
#'
#' @examples
compute_error <- function(target_mean, target_cov, mod_mean1, mod_mean2, mod_cov1, mod_cov2, start.coord, box.sd, box.step, lower) {
  library(numDeriv)
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
  
  coord_data <- mpfr(as.matrix(all.bound[, 1:3]), precBits = 60)
  #coord_data <- as.matrix(all.bound[, 1:3])
  error <- sum(width.frac^2 * (Rmpfr::pnorm(coord_data[,dim_cum], 0, 1, lower.tail = lower) * Rmpfr::dnorm(coord_data[,dim_marg[1]], 0, 1) * Rmpfr::dnorm(coord_data[,dim_marg[2]], 0, 1)))
  #error <- sum(width.frac^2 * (stats::pnorm(coord_data[,dim_cum], 0, 1, lower.tail = lower) * stats::dnorm(coord_data[,dim_marg[1]], 0, 1) * stats::dnorm(coord_data[,dim_marg[2]], 0, 1)))
}

#' Compute an error rate for two gaussian distributions.
#' The discriminant bound can be suboptimal. This is computed mean mod != target
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
#'
#' @examples
compute_error.2 <- function(target_mean, target_cov, mod_mean1, mod_mean2, mod_cov1, mod_cov2, start.coord, box.sd, box.step, lower) {
  library(numDeriv)
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
  
  coord_data <- mpfr(as.matrix(all.bound[, 1:3]), precBits = 60)
  #coord_data <- as.matrix(all.bound[, 1:3])
  error <- sum(width.frac^2 * (Rmpfr::pnorm(coord_data[,dim_cum], 0, 1, lower.tail = lower) * Rmpfr::dnorm(coord_data[,dim_marg[1]], 0, 1) * Rmpfr::dnorm(coord_data[,dim_marg[2]], 0, 1)))
  #error <- sum(width.frac^2 * (stats::pnorm(coord_data[,dim_cum], 0, 1, lower.tail = lower) * stats::dnorm(coord_data[,dim_marg[1]], 0, 1) * stats::dnorm(coord_data[,dim_marg[2]], 0, 1)))
}

