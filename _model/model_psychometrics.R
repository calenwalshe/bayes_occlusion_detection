#' Compute eccentricity psychometric functions for the model responses
#'
get_model_psychometric <- function(model.dprime, scale.dprime = 1) {
  library(bbmle)
  library(dplyr)
  library(purrr)
  
  model.dprime <- model.dprime %>% 
    mutate(dprime = dprime * scale.dprime) %>%
    mutate(percent_correct = pnorm(dprime/2)) %>%
    select(BIN, TARGET, eccentricity, dprime, SUBJECT) %>%
    unique()
  
  fitted <- model.dprime %>% 
    as_tibble() %>%
    group_by(BIN, TARGET, SUBJECT) %>%
    arrange(eccentricity) %>%
    mutate(d0 = max(dprime)) %>%
    group_by(BIN, TARGET, SUBJECT, d0) %>%
    nest() %>%
    mutate(models = map(data, function(x) { 
      model <- mle2(dprime ~ dnorm(max(dprime) * e0^b/(e0^b + eccentricity^b), sd = 1),
                     start = list(b = 2, e0 = 10),
                     data = x)
    })) %>%
    select(-data)
  
  fitted.params <- fitted %>%
    mutate(params = map(models, coef), 
           b = map(params, c(1)),
           e0 = map(params, c(2)),
           gamma = 0) %>%
    select(-models, -params) %>%
    unnest()
  
  return(fitted.params)
}

fit.scaled.thresholds <- function(human.psychometrics, model.dprime) {
  
  human.thresholds <- get_threshold(human.psychometrics) %>%
    select(TARGET, SUBJECT, BIN, threshold)
  
  f <- function(a) {
    pf <- get_model_psychometric(model.dprime, a)
    
    model <- pf %>% arrange(TARGET, BIN, SUBJECT) %>% get_threshold(.) %>%
      select(TARGET, SUBJECT, BIN, threshold) %>% full_join(human.thresholds, ., by = c("TARGET", "BIN"))
    
    obj <- sum((model$threshold.x - model$threshold.y)^2)
  }
  
  m1 <- optimize(f, lower = c(a = .001), upper = c(a = .25))
  
  return(m1$minimum)
}

