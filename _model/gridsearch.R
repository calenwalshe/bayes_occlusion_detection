#' Get best parameters via grid search for the full covariance Gaussian model.
#'
get_model_grid <- function(template.response, human.psychometrics, 
    level_vec, error_fcn = "mle") {
    library(parallel)
    ecc_dat <- template.response %>% dplyr::select(eccentricity) %>% 
        distinct() %>% 
      mutate(eccentricity)
    
    level_vec <- sort(level_vec)
    
    human.psychometrics <- human.psychometrics[, setdiff(names(human.psychometrics), 
        c("L", "C", "S", "statType", "statValue"))] %>% 
      distinct()
    
    human_pc <- merge(ecc_dat, human.psychometrics) %>% 
      group_by(eccentricity, 
        TARGET, BIN, SUBJECT) %>% 
      rowwise() %>% 
      mutate(human_dprime = d0 * e0^b/(e0^b + eccentricity^b), human_pc = pnorm(human_dprime/2), 
        human_pc = ifelse(human_pc == 1, 1 - 1/300, human_pc)) %>% 
        arrange(TARGET, BIN, eccentricity) %>% 
      ungroup()
    
    grid_f <- function(x) {
        if (error_fcn == "mle") {
            error <- get_NLL(x, template.response, human_pc)
        } else {
            error <- get_E(x, template.response, human_pc)
        }
      
        grid_search <- expand.grid(seq(0.001, 1, length.out = 3), 
            seq(0.001, 1, length.out = 10), seq(0.001, 1, 
                length.out = 10), seq(0, 0, length.out = 2))
        
        xy.list <- as.list(as.data.frame(t(grid_search)))
        parameters <- unlist(mclapply(xy.list, FUN = function(x) return(tryCatch(error(x[1], 
            x[2], x[3], x[4]), error = function(e) NULL)), mc.cores = 16))
        
        grid_search$error <- parameters
        
        b_params <- grid_search %>% 
          mutate(eccentricity = x) %>% 
            filter(!is.na(error)) %>% 
          filter(error == min(error)) %>% 
            head(., 1) %>% 
          print
    }
    
    full_cov_best_params <- lapply(level_vec, FUN = grid_f)
    
    full_df <- (do.call(rbind, full_cov_best_params))
    
    names(full_df) <- c("edge", "lum", "pattern", "bias", "error", "eccentricity")
    
    return(full_df)
}


