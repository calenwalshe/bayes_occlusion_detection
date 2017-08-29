#' Compute eccentricity psychometric functions for the model responses
#'
get_model_psychometric <- function(model_dprime) {
    # Compute negative log likelihood.
    get_NLL <- function(current_data) {
        NLL <- function(d0, b, e0) {
            current_data <- current_data %>% mutate(model_logLk_correct = pnorm(0.5 * 
                d0 * (e0^b/(e0^b + eccentricity^b)), log = T, 
                low = T), model_logLk_incorrect = pnorm(0.5 * 
                d0 * (e0^b/(e0^b + eccentricity^b)), log = T, 
                low = F))
            
            NLL_val <- current_data %>% summarize(NLL = sum(percent_correct * 
                -model_logLk_correct) + sum((1 - percent_correct) * 
                -model_logLk_incorrect))
            return(NLL_val)
        }
    }
    
    model_dprime <- model_dprime %>% arrange(BIN) %>% select(-statType, 
        -statValue, -L, -C, -S) %>% distinct()
    
    model_dprime <- model_dprime %>% select(-eccentricity, -dprime, 
        -percent_correct) %>% distinct() %>% mutate(eccentricity = 0, 
        dprime = 4.5, percent_correct = pnorm(dprime/2)) %>% 
        rbind(., model_dprime)
    
    model_dprime <- model_dprime %>% select(-eccentricity, -dprime, 
        -percent_correct) %>% distinct() %>% mutate(eccentricity = 30, 
        dprime = 0, percent_correct = pnorm(dprime/2)) %>% rbind(., 
        model_dprime)
    
    model_dprime <- model_dprime %>% arrange(BIN, TARGET, eccentricity) %>% 
        mutate(function_name = "combined")
    
    model_dprime <- data.frame(model_dprime)
    param_df <- data.frame(BIN = NULL, TARGET = NULL, function_name = NULL, 
        e0 = NULL, b = NULL, d0 = NULL, gamma_bias = NULL)
    
    if (!any(model_dprime$function_name == "combined")) {
        starts <- list(BIN = c(1, 4, 4, 8), TARGET = c("bowtie", 
            "spot", "vertical", "bowtie"), function_name = c("pattern_only", 
            "edge_cos", "edge_cos", "pattern_only"), b = c(1, 
            1, 1, 1), e0 = c(1, 1, 1, 1))
        start_params <- get_start_params(model_dprime, starts)
        
    } else {
        
        starts <- list(BIN = c(10, 10, 9, 9, 8, 15, 15, 8, 14, 
            5, 9, 6, 10, 9, 6, 5, 5, 5, 2, 3, 4, 4, 4, 4, 7), 
            TARGET = c("bowtie", "vertical", "horizontal", "vertical", 
                "bowtie", "bowtie", "vertical", "vertical", "bowtie", 
                "spot", "bowtie", "spot", "spot", "spot", "vertical", 
                "bowtie", "horizontal", "vertical", "bowtie", 
                "bowtie", "vertical", "horizontal", "spot", "bowtie", 
                "bowtie"), function_name = c("combined", "combined", 
                "combined", "combined", "combined", "combined", 
                "combined", "combined", "combined", "combined", 
                "combined", "combined", "combined", "combined", 
                "combined", "combined", "combined", "combined", 
                "combined", "combined", "combined", "combined", 
                "combined", "combined", "combined"), b = c(1, 
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                1, 1, 10, 10, 1, 1, 1, 1, 1), e0 = c(1, 1, 1, 
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                8, 5, 3, 1, 1, 1, 1, 1))
        start_params <- get_start_params(model_dprime, starts)
    }
    
    for (i in unique(model_dprime$BIN)) {
        for (j in unique(model_dprime$TARGET)) {
            for (k in unique(model_dprime$function_name)) {
                
                current_data <- model_dprime %>% filter(BIN == 
                  i, TARGET == j, function_name == k)
                
                NLL <- get_NLL(current_data)
                print(i)
                print(j)
                print(k)
                
                d0 <- max(current_data$dprime)
                
                
                start_b <- start_params[start_params$BIN == i & 
                  start_params$TARGET == j & start_params$function_name == 
                  k, ]$b
                
                start_e0 <- start_params[start_params$BIN == 
                  i & start_params$TARGET == j & start_params$function_name == 
                  k, ]$e0
                
                maxLik <- mle(NLL, start = list(b = start_b, 
                  e0 = start_e0), method = "L-BFGS-B", lower = c(0.001, 
                  0.001), upper = c(40, 40), fixed = list(d0 = d0), 
                  control = list(trace = 2))
                
                param_df <- rbind(param_df, data.frame(BIN = i, 
                  TARGET = j, function_name = k, d0 = coef(maxLik)[1], 
                  b = coef(maxLik)[2], e0 = coef(maxLik)[3], 
                  row.names = NULL, log_likelihood = logLik(maxLik), 
                  gamma_bias = 0))
            }
        }
    }
    
    param_df$BIN <- as.factor(param_df$BIN)
    param_df$SUBJECT <- "model"
    return(param_df)
}

get_start_params <- function(model_dprime, custom_starts) {
    
    start_params <- model_dprime %>% select(TARGET, BIN, function_name) %>% 
        distinct()
    
    start_params$b <- 3
    start_params$e0 <- 8
    
    if (!missing(custom_starts)) {
        for (i in seq(1, length(custom_starts$TARGET))) {
            BIN <- custom_starts$BIN[i]
            TARGET <- custom_starts$TARGET[i]
            function_name <- custom_starts$function_name[i]
            
            if (length(start_params[start_params$BIN == BIN & 
                start_params$TARGET == TARGET & start_params$function_name == 
                function_name, ]$b) > 0) {
                
                start_params[start_params$BIN == BIN & start_params$TARGET == 
                  TARGET & start_params$function_name == function_name, 
                  ]$b <- custom_starts$b[i]
                start_params[start_params$BIN == BIN & start_params$TARGET == 
                  TARGET & start_params$function_name == function_name, 
                  ]$e0 <- custom_starts$e0[i]
            }
        }
    }
    
    
    
    return(start_params)
    
}


#' Return model parameters fit with maximum likelihood
get_model_mle_params <- function(model_responses) {
    library(bbmle)
    library(parallel)
    
    f <- function(d0, e0, b, g) {
        likelihood <- sum(log(pnorm(0.5 * d0 * (e0^b)/(e0^b + 
            hit_vec^b) - g))) + sum(log(pnorm(-0.5 * d0 * (e0^b)/(e0^b + 
            fa_vec^b) - g))) + sum(log(pnorm(-0.5 * d0 * (e0^b)/(e0^b + 
            miss_vec^b) + g))) + sum(log(pnorm(0.5 * d0 * (e0^b)/(e0^b + 
            cr_vec^b) + g)))
        
        nll <- -likelihood
        
        if (is.infinite(nll)) {
            return(10^6)
        } else {
            return(nll)
        }
        
    }
    
    get_params <- function(model_response_condition) {
        hit_vec <- model_response_condition[model_response_condition$HIT == 
            1, "ECCENTRICITY"]
        fa_vec <- model_response_condition[model_response_condition$FALSEALARM == 
            1, "ECCENTRICITY"]
        cr_vec <- model_response_condition[model_response_condition$CORRECTREJECTION == 
            1, "ECCENTRICITY"]
        miss_vec <- model_response_condition[model_response_condition$MISS == 
            1, "ECCENTRICITY"]
        
        environment(f) <- environment()
        return(f)
        
    }
    
    model_response_list <- model_responses %>% group_by(TARGET, 
        BIN, SUBJECT) %>% do(vals = data.frame(.)) %>% select(vals) %>% 
        as.list()
    model_response_list <- model_response_list[[1]]
    
    fcn_vec <- lapply(model_response_list, FUN = function(x) get_params(x))
    
    start_params <- expand.grid(d0 = 4.5, e0 = seq(1, 20, 0.5), 
        b = seq(0, 5, 0.5), g = seq(-10, 10, 2))
    start_params <- data.frame(d0 = 4.5, e0 = 10, b = 2.5, g = 0)
    n_search <- nrow(start_params)
    
    
    
    # Grid search for best starting parameters
    grid_params <- lapply(fcn_vec, FUN = function(x) cbind(start_params, 
        nll = unlist(mclapply(1:n_search, FUN = function(y) x(d0 = start_params[y, 
            1], e0 = start_params[y, 2], b = start_params[y, 
            3], g = start_params[y, 4]), mc.cores = 16))) %>% 
        filter(nll == min(nll)))
    
    start_params <- do.call(rbind, grid_params)
    subject_df <- lapply(model_response_list, FUN = function(x) x[1, 
        c("TARGET", "BIN", "SUBJECT")]) %>% do.call(rbind, .)
    
    # Maximum likelihood for parameters, all parameters free to
    # vary.
    p.1 <- mclapply(1:length(fcn_vec), FUN = function(x) {
        y <- start_params[x, ]
        mle2(fcn_vec[[x]], start = list(e0 = y$e0, b = y$b, g = y$g), 
            lower = c(0, 0, -20), fixed = list(d0 = 4.5)) %>% 
            coef
    }) %>% do.call(rbind, .) %>% data.frame %>% cbind(subject_df, 
        .)  # 
    
    
    fitted.params <- p.1
    return(fitted.params)
}