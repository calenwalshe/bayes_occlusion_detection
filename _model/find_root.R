# Find roots of the descriminant function.
# Description: Step candidate points for y, z and look for the root at x.
find_root <- function(mean1, mean2, cov1, cov2) {
  inv_cov1 <- solve(cov1)
  inv_cov2 <- solve(cov2)
  
  mean_cov1 <- inv_cov1 %*% mean1
  mean_cov2 <- inv_cov2 %*% mean2
  
  cov_diff <- -.5 * (inv_cov1 - inv_cov2)
  mean_cov_diff <- mean_cov1 - mean_cov2
  
  a_1 <- cov_diff[1,1]
  a_2 <- cov_diff[2,2]
  a_3 <- cov_diff[3,3]
  
  b_1 <- 2 * cov_diff[1,2]
  b_2 <- 2 * cov_diff[1,3]
  b_3 <- 2 * cov_diff[2,3]
  
  c_1 <- mean_cov_diff[1]
  c_2 <- mean_cov_diff[2]
  c_3 <- mean_cov_diff[3]
  
  d <- -.5 * (t(mean1) %*% inv_cov1 %*% (mean1)) - .5 * log(det(cov1)) + .5 * (t(mean2) %*% inv_cov2 %*% (mean2)) + .5 * log(det(cov2))
  
  y_min <- min(c(mean1[2] - sqrt(cov1[2,2]), mean2[2] - sqrt(cov2[2,2])))
  y_max <- max(c(mean1[2] + sqrt(cov1[2,2]), mean2[2] + sqrt(cov2[2,2])))
  
  z_min <- min(c(mean1[3] - sqrt(cov1[3,3]), mean2[3] - sqrt(cov2[3,3])))
  z_max <- max(c(mean1[3] + sqrt(cov1[3,3]), mean2[3] + sqrt(cov2[3,3])))    
  
  g <- function(y,z) {
    term1 <- (-b_1 * y - b_2 * z - c_1)
    term2 <- sqrt(
      (b_1 * y + b_2 * z + c_1)^2 - 4 * a_1 * (a_2 * y^2 + a_3 * z^2 + b_3 * z*y + c_2 * y + c_3 * z + d))
    
    val1 <- (term1 + term2) / (2*a_1)
    val2 <- (term1 - term2) / (2*a_1)
    return(matrix(c(val1, val2, y, y, z, z), nrow = 2))
  }
  
  search.points <- expand.grid(y = seq(y_min, y_max, length.out = 80), z = seq(z_min, z_max, length.out = 80))
  
  roots <- do.call(rbind, do.call(rbind, apply(search.points, 1, FUN = function(x) list(g(x[1], x[2])))))
  
  discriminant.dist <- 
    apply(roots, 1, FUN = function(x) 
      dmvnorm(x, mean1, cov1, log = T))
  
  ind <- first(which(discriminant.dist == max(discriminant.dist, na.rm = T)))
  
  max_coord <- roots[ind, ]
  max_dist <- discriminant.dist[ind]
  
  
  return(c(max_coord, max_dist))
}  

# Find roots through an optimization method. Use roots found through step based approach.
# Description: Stepping is time consuming. We use an optimization method to search through a
# larger region of feature space to detect a high likelihood root.
find_roots_optim <- function(start_val = c(0,0,0), mean1, mean2, cov1, cov2) {
  inv_cov1 <- solve(cov1)
  inv_cov2 <- solve(cov2)
  
  mean_cov1 <- inv_cov1 %*% mean1
  mean_cov2 <- inv_cov2 %*% mean2
  
  cov_diff <- -.5 * (inv_cov1 - inv_cov2)
  mean_cov_diff <- mean_cov1 - mean_cov2
  
  a_1 <- cov_diff[1,1]
  a_2 <- cov_diff[2,2]
  a_3 <- cov_diff[3,3]
  
  b_1 <- 2 * cov_diff[1,2]
  b_2 <- 2 * cov_diff[1,3]
  b_3 <- 2 * cov_diff[2,3]
  
  c_1 <- mean_cov_diff[1]
  c_2 <- mean_cov_diff[2]
  c_3 <- mean_cov_diff[3]
  
  d <- -.5 * (t(mean1) %*% inv_cov1 %*% (mean1)) - .5 * log(det(cov1)) + .5 * (t(mean2) %*% inv_cov2 %*% (mean2)) + .5 * log(det(cov2))
  
  y_min <- min(c(start_val[2] - 3*sqrt(cov1[2,2]), start_val[2] - 3*sqrt(cov2[2,2])))
  y_max <- max(c(start_val[2] + 3*sqrt(cov1[2,2]), start_val[2] + 3*sqrt(cov2[2,2])))
  
  z_min <- min(c(start_val[3] - 3*sqrt(cov1[3,3]), start_val[3] - 3*sqrt(cov2[3,3])))
  z_max <- max(c(start_val[3] + 3*sqrt(cov1[3,3]), start_val[3] + 3*sqrt(cov2[3,3])))    
  
  g <- function(x) {
    y <- x[1]
    z <- x[2]
    term1 <- (-b_1 * y - b_2 * z - c_1)
    term2 <- sqrt(
      (b_1 * y + b_2 * z + c_1)^2 - 4 * a_1 * (a_2 * y^2 + a_3 * z^2 + b_3 * z*y + c_2 * y + c_3 * z + d))
    
    val1 <- (term1 + term2) / (2*a_1)
    val2 <- (term1 - term2) / (2*a_1)

    d1 <- dmvnorm(c(val1, y, z), mean1, cov1, log = T)
    d2 <- dmvnorm(c(val2, y, z), mean1, cov1, log = T)

    
    m_dist_min <- -max(d1, d2, na.rm =T)
    
    if(m_dist_min == Inf) {
      m_dist_min <- 10^5
    }
    
    return(m_dist_min)
  }
  
  find_single_root <- function(y,z) {
    term1 <- (-b_1 * y - b_2 * z - c_1)
    term2 <- sqrt(
      (b_1 * y + b_2 * z + c_1)^2 - 4 * a_1 * (a_2 * y^2 + a_3 * z^2 + b_3 * z*y + c_2 * y + c_3 * z + d))
    
    val1 <- (term1 + term2) / (2*a_1)
    val2 <- (term1 - term2) / (2*a_1)
    return(matrix(c(val1, val2, y, y, z, z), nrow = 2))
  }
  
  #result <- optim(start_val[2:3], g, lower = c(y_min, z_min), upper = c(y_max, z_max), method = "L-BFGS-B")
  
  result <- DEoptim(g, lower = c(y_min, z_min), upper = c(y_max, z_max))
  root.val <- find_single_root(y = result$optim$bestmem[1], z = result$optim$bestmem[2])
  
  d1 <- dmvnorm(root.val[1,], mean1, cov1, log = T)
  d2 <- dmvnorm(root.val[2,], mean1, cov1, log = T)
  
  if(d2 < d1 | d2 == Inf) {
    return(c(root.val[1,], d1))
  } else {
    return(c(root.val[2,], d2))
  }
}
