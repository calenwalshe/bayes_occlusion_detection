  # Find roots through an optimization method. Use roots found through step based approach.
  # Description: Stepping is time consuming. We use an optimization method to search through a
  # larger region of feature space to detect a high likelihood root.
  find_roots_suboptim <- function(start_val = c(0,0,0), mean1, mean2, target_mean, cov1, cov2, target_cov) {
    ## Transform Order ##
    dim_order <- order(abs(mean1), decreasing = T)
    new_cov <- lapply(list(cov1, cov2, target_cov), FUN = function(x) {
      upper_tri <- x[upper.tri(x)==T]
      new_mat <- diag(3)
      new_mat[1,1] <- diag(x)[dim_order[1]]
      new_mat[2,2] <- diag(x)[dim_order[2]]
      new_mat[3,3] <- diag(x)[dim_order[3]]
      
      new_mat[1,2] <- x[dim_order[1], dim_order[2]]
      new_mat[1,3] <- x[dim_order[1], dim_order[3]]
      new_mat[2,3] <- x[dim_order[2], dim_order[3]]
      
      new_mat[lower.tri(new_mat)==1] <- new_mat[upper.tri(new_mat)==1]
      
      return(new_mat)
    })
    
    #cov1 <- new_cov[[1]]
    #cov2 <- new_cov[[2]]
    #target_cov <- new_cov[[3]]
    
    #mean1 <- mean1[dim_order]
    #mean2 <- mean2[dim_order]
    #target_mean <- target_mean[dim_order]
    #start_val <- start_val[dim_order]
    ##
    
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
      term2 <- suppressWarnings(
        sqrt(
        (b_1 * y + b_2 * z + c_1)^2 - 4 * a_1 * (a_2 * y^2 + a_3 * z^2 + b_3 * z*y + c_2 * y + c_3 * z + d))
      )
      
      val1 <- (term1 + term2) / (2*a_1)
      val2 <- (term1 - term2) / (2*a_1)
      if(any(is.nan(c(val1, val2)))) {
        return(Inf)
      }else{
        d1 <- dmvnorm(c(val1, y, z), target_mean, target_cov, log = T)
        d2 <- dmvnorm(c(val2, y, z), target_mean, target_cov, log = T)
        
        return(-max(d1, d2))
      }
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
    
    result <- DEoptim(g, lower = c(y_min, z_min), upper = c(y_max, z_max), DEoptim.control(reltol = 10e-4))
    root.val <- find_single_root(y = result$optim$bestmem[1], z = result$optim$bestmem[2])
    
    d1 <- dmvnorm(root.val[1,], mean1, cov1, log = T)
    d2 <- dmvnorm(root.val[2,], mean1, cov1, log = T)
    
    if(d2 < d1 | d2 == Inf) {
      return(c(root.val[1,], d1))
      #return(c(root.val[1,dim_order], d1))
    } else {
      return(c(root.val[2,], d2))
      #return(c(root.val[2,dim_order], d2))
    }
  }

  # Find roots through an optimization method. Use roots found through step based approach.
  # Description: Stepping is time consuming. We use an optimization method to search through a
  # larger region of feature space to detect a high likelihood root.
  find_roots_suboptim_ordered <- function(start_val = c(0,0,0), mean1, mean2, target_mean, cov1, cov2, target_cov) {
    ## Transform Order ##
    dim_order <- order(abs(mean1), decreasing = T)
    new_cov <- lapply(list(cov1, cov2, target_cov), FUN = function(x) {
      upper_tri <- x[upper.tri(x)==T]
      new_mat <- diag(3)
      new_mat[1,1] <- diag(x)[dim_order[1]]
      new_mat[2,2] <- diag(x)[dim_order[2]]
      new_mat[3,3] <- diag(x)[dim_order[3]]
      
      new_mat[1,2] <- x[dim_order[1], dim_order[2]]
      new_mat[1,3] <- x[dim_order[1], dim_order[3]]
      new_mat[2,3] <- x[dim_order[2], dim_order[3]]
      
      new_mat[lower.tri(new_mat)==1] <- new_mat[upper.tri(new_mat)==1]
      
      return(new_mat)
    })
    
    cov1 <- new_cov[[1]]
    cov2 <- new_cov[[2]]
    target_cov <- new_cov[[3]]
    
    mean1 <- mean1[dim_order]
    mean2 <- mean2[dim_order]
    target_mean <- target_mean[dim_order]
    start_val <- start_val[dim_order]
    ##
    
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
      term2 <- suppressWarnings(
        sqrt(
          (b_1 * y + b_2 * z + c_1)^2 - 4 * a_1 * (a_2 * y^2 + a_3 * z^2 + b_3 * z*y + c_2 * y + c_3 * z + d))
      )
      
      val1 <- (term1 + term2) / (2*a_1)
      val2 <- (term1 - term2) / (2*a_1)
      if(any(is.nan(c(val1, val2)))) {
        return(Inf)
      }else{
        d1 <- dmvnorm(c(val1, y, z), target_mean, target_cov, log = T)
        d2 <- dmvnorm(c(val2, y, z), target_mean, target_cov, log = T)
        
        return(-max(d1, d2))
      }
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
    
    result <- DEoptim(g, lower = c(y_min, z_min), upper = c(y_max, z_max), DEoptim.control(reltol = 10e-4))
    root.val <- find_single_root(y = result$optim$bestmem[1], z = result$optim$bestmem[2])
    
    d1 <- dmvnorm(root.val[1,], mean1, cov1, log = T)
    d2 <- dmvnorm(root.val[2,], mean1, cov1, log = T)
    
    reorder.1 <- c(0,0,0)
    reorder.2 <- c(0,0,0)
    reorder.1[dim_order] <- root.val[1,]
    reorder.2[dim_order] <- root.val[2,]
    
    if(d2 < d1 | d2 == Inf) {
      #return(c(root.val[1,], d1))
      return(c(reorder.1, d1))
    } else {
      #return(c(root.val[2,], d2))
      return(c(reorder.2, d2))
    }
  }
  
  