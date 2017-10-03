#' Return the percent correct value for human observers at all levels of the experiment.
get_human_pc <- function(template.response, human.psychometrics) {
  ecc_dat <- template_response %>% dplyr::select(eccentricity) %>% 
    distinct()
  
  human_psychometrics <- human_psychometrics[, setdiff(names(human_psychometrics), c("L", "C", "S", "statType", "statValue"))] %>%
    distinct()
  
  human_pc <- merge(ecc_dat, human_psychometrics) %>% 
    group_by(eccentricity, TARGET, BIN, SUBJECT) %>%
    rowwise() %>%
    mutate(human_dprime = d0 * e0^b/(e0^b + eccentricity^b), human_pc = pnorm(human_dprime/2), human_pc = ifelse(human_pc == 1, 1 - 1/300, human_pc)) %>%
    arrange(TARGET, BIN, eccentricity)
}

#' Search a global optimum.
get_fitted_search <- function(template_response, human_psychometrics, start_params) {
  library(bbmle)
  library(DEoptim)
  
  human_pc <- get_human_pc(template_response, human.psychometrics)
  
  f <- function(a, b, c0, d, ECC, template_response) {
      error <- get_NLL(ECC, template_response, human_pc)
      
      f_diff_error <- function(x) {
        error(x[1], x[2], x[3], 0)
      }
      
      params <- DEoptim(f_diff_error, lower = c(0, 0, 0), upper = c(2, 2, 2), DEoptim.control(NP = 80, F = 0.8, CR = 0.9, itermax = 70))
      
      #params <- cma_es(c(a,b,c0), f_cma_error, lower = c(0,0,0), upper = c(10000,10000,100000))
      
      return(c(params$optim$bestmem[1], params$optim$bestmem[2], params$optim$bestmem[3], 0, params$optim$bestval, ECC))
  }
  
  parameters_n <- nrow(start_params)  
  
  mle_parameters <- mclapply(1:parameters_n, FUN = function(x) 
    f(start_params[x,][[1]],
      start_params[x, ][[2]],
      start_params[x, ][[3]], 
      start_params[x, ][[4]], 
      start_params[x, ][[6]], 
      template_response), mc.cores = 16)
  
  mle_parameters_df <- do.call(rbind, mle_parameters)
  mle_parameters_df <- data.frame(mle_parameters_df)
  
  names(mle_parameters_df) <- c("edge", "luminance", "pattern", 
                                "bias", "error", "eccentricity")
  
  return(mle_parameters_df)
}

#' Provide a parameter set and get back the responses predicted from the model
get_model_responses <- function(parameters, template_response, 
                                human_psychometrics) {
  ecc_dat <- template_response %>%
    dplyr::select(eccentricity) %>% 
    distinct()

  human_psychometrics <- human_psychometrics[, setdiff(names(human_psychometrics), 
                                                       c("L", "C", "S", "statType", "statValue"))] %>%
    distinct()
  
  human_pc_exp <- merge(ecc_dat, human_psychometrics) %>%
    group_by(eccentricity, TARGET, BIN) %>%
    rowwise() %>% 
    mutate(human_dprime = d0 * e0^b/(e0^b + eccentricity^b), human_pc = pnorm(human_dprime/2), 
           human_pc = ifelse(human_pc == 1, 1 - 1/300, human_pc)) %>% 
    arrange(TARGET, BIN, eccentricity)
  
  get_response <- function(a, b, c0, d, ECC) {
    
    get_response_f <- get_model_response_fcn(ECC = ECC, template_response)
    
    return(get_response_f(a, b, c0, d))
  }
  
  param_len <- nrow(parameters)
  
  model_responses <- mclapply(1:param_len, FUN = function(x) 
    get_response(parameters[x,][[1]], parameters[x, ][[2]], parameters[x, ][[3]], parameters[x, ][[4]], parameters[x, ][[6]]),
    mc.cores = 16)
  
  model_responses <- do.call(rbind, model_responses)
  
  bin_data <- get_experiment_bin_values() %>% merge(., template_response %>% 
                                                      select(eccentricity) %>% distinct())
  
  model_responses <- model_responses %>% mutate(function_name = "combined") %>% 
    arrange(TARGET, BIN, ECCENTRICITY)
  return(model_responses)
}

#' Sets up model performance function for a given eccentricity. Human responses are also supplied and fit values are returned.
get_model_response_fcn <- function(ECC, template_response) {
  load(paste0("~/Dropbox/Calen/Work/data_storage/eccentricity_project/model_distribution_", floor(ECC), ".rdata"))

  environment(get_model_response_weights) <- environment()
  return(get_model_response_weights)
}

#' Returns summary performance on a trial by trial basis.
get_model_response_weights <- function(a, b, c0, d) {
  scaled_model_distribution <- model.data
  
  iter <- seq(1, (nrow(scaled_model_distribution) - 1), 2)
  resp_func <- function(i) {
    BIN_id <- as.numeric(scaled_model_distribution[[i, "BIN"]])
    TARGET_id <- as.character(scaled_model_distribution[[i, 
                                                         "TARGET"]])
    
    responses.present <- scaled_model_distribution[[i,"response_vec"]]
    responses.absent  <- scaled_model_distribution[[i + 1,"response_vec"]]
    
    responses         <- rbind(responses.present, responses.absent)
    
    correct_responses <- c(rep(0, nrow(responses)/2), rep(1, nrow(responses)/2))
    
    mean_vec_abs <- scaled_model_distribution[[i, "mean_mat"]]
    mean_vec_pres <- scaled_model_distribution[[i + 1, "mean_mat"]]
    
    covmat_abs <- scaled_model_distribution[[i, "cov_mat"]]
    covmat_pres <- scaled_model_distribution[[i + 1, "cov_mat"]]
    
    # Additive Noise
    num.response <- nrow(responses)
    response.noise <- responses
    response.noise[, 1] <- response.noise[, 1] + rnorm(num.response, 0, a)
    response.noise[, 2] <- response.noise[, 2] + rnorm(num.response, 0, b)
    response.noise[, 3] <- response.noise[, 3] + rnorm(num.response, 0, c0)
    
    covmat_pres <- cov(response.noise[correct_responses == 1, ])
    covmat_abs <- cov(response.noise[correct_responses == 0, ])
    
    # Additive Noise
    response.noise <- responses
    response.noise[, 1] <- response.noise[, 1] + rnorm(num.response, 0, a)
    response.noise[, 2] <- response.noise[, 2] + rnorm(num.response, 0, b)
    response.noise[, 3] <- response.noise[, 3] + rnorm(num.response, 0, c0)
    
    response.mat <- response.noise
    
    discriminant <- mnormt::dmnorm(x = response.mat, mean = mean_vec_pres, 
                                   varcov = covmat_pres, log = TRUE) - mnormt::dmnorm(x = response.mat, 
                                                                                      mean = mean_vec_abs, varcov = covmat_abs, log = TRUE)
    m.response <- as.numeric(discriminant >= d)
    
    hits <- as.numeric(m.response == 1 & correct_responses == 
                         1)
    miss <- as.numeric(m.response == 0 & correct_responses == 
                         1)
    fa <- as.numeric(m.response == 1 & correct_responses == 
                       0)
    cr <- as.numeric(m.response == 0 & correct_responses == 
                       0)
    correct <- as.numeric(m.response == correct_responses)
    
    model_results <- NULL
    model_results$HIT <- hits
    model_results$FALSEALARM <- fa
    model_results$MISS <- miss
    model_results$CORRECTREJECTION <- cr
    model_results$CORRECT <- correct
    model_results$BIN <- BIN_id
    model_results$TARGET <- TARGET_id
    model_results$ECCENTRICITY <- ECC
    model_results$SUBJECT <- "model"
    model_results <- data.frame(model_results)
    
    model_results <- model_results %>% select(TARGET, BIN, 
                                              SUBJECT, ECCENTRICITY, HIT, FALSEALARM, MISS, CORRECTREJECTION, 
                                              CORRECT)
    return(model_results)
  }
  
  responses <- lapply(iter, FUN = function(x) resp_func(x))
  
  responses <- do.call(rbind, responses)
  
  return(responses)
}

#' Get a dataframe of model responses
get_model_results <- function(ECC, template_response, human_pc_exp) {
  response_fcn <- get_model_response_fcn(ECC, template_response)
  
  human_pc <- human_pc_exp %>% filter(eccentricity == ECC) %>% 
    select(BIN, TARGET, eccentricity, SUBJECT, human_pc)
  
  f <- function(a, b, c0, d) {
    responses <- response_fcn(a, b, c0, d)
    
    r.1 <- responses %>% group_by(BIN, TARGET, ECCENTRICITY) %>% 
      summarize(HIT = sum(HIT)/(n()/2), FALSEALARM = sum(FALSEALARM)/(n()/2), 
                pc = mean(CORRECT)) %>% mutate(HIT = ifelse(HIT == 0, 1/600, HIT),
                                               FALSEALARM = ifelse(FALSEALARM == 0,
                                                                   1/600, FALSEALARM),
                                               HIT = ifelse(HIT == 1, HIT - 1/600, HIT),
                                               FALSEALARM = ifelse(FALSEALARM == 1, FALSEALARM - 1/600,
                                                                   FALSEALARM), 
                                               pc = ifelse(pc == 1, pc - 1/1000, pc), pc = ifelse(pc == 0, 1/1000, pc)) %>%
      mutate(dprime = qnorm(HIT) - qnorm(FALSEALARM)) %>% 
      merge(human_pc, .) %>% 
      mutate(NLL = -(human_pc *log(pc) + (1 - human_pc) * log(1 - pc)), E = abs(human_pc - pc), ESQ = (human_pc - pc)^2)
    return(r.1)
  }
}

#' Computes a negative log likelihood.
get_NLL <- function(ECC, template_response, human_pc_exp) {
  f <- get_model_results(ECC, template_response, human_pc_exp)
  NLL <- function(a, b, c0, d) {
    if (a <= 0 || b <= 0 || c0 <= 0) {
      return(1e+06)
    }
    
    model_response <- f(a, b, c0, d)
    NLL_val <- sum(model_response$NLL)
    return(NLL_val)
  }
}

#' Return an error function based on the squared distance between the model and human percent correct.
get_E <- function(ECC, template_response, human_pc_exp) {
  f <- get_model_results(ECC, template_response, human_pc_exp)
  E <- function(a, b, c, d) {
    if (a == 0 || b == 0 || c == 0) {
      return(1e+06)
    }
    model_response <- f(a, b, c, d)
    e_val <- sum(model_response$ESQ)
    return(e_val)
  }
}

#' Computes a negative log likelihood with only a single fixed value at a given eccentricity.
get_NLL_fixed <- function(ECC, template_response, human_pc_exp) {
  f <- get_model_results(ECC, template_response, human_pc_exp)
  NLL <- function(a, d) {
    
    a <- a
    b <- a
    c0 <- a
    d <- d
    
    model_response <- f(a, b, c0, d)
    NLL_val <- sum(model_response$NLL)
    return(NLL_val)
  }
}
