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
           e0 = map(params, c(2))) %>%
    select(-models, -params) %>%
    unnest()
  
  return(fitted.params)
  
}

fit.scaled.thresholds <- function(human.psychometrics, model.dprime) {
  
  human.thresholds <- get_threshold(human.psychometrics)
  
  f <- function(a) {
    pf <- get_model_psychometric(model.dprime, a)
    
    model <- pf %>% arrange(TARGET, BIN, SUBJECT) %>% get_threshold(.)
    
    human <- human.thresholds %>% arrange(TARGET, BIN, SUBJECT)
    
    NLL <- sum(-dnorm(model$threshold, mean = human$threshold, sd = 1, log = TRUE))
    
  }
  
  m1 <- mle2(f, start = list(a = .5), data = model.dprime)
  
  return(coef(m1)[1])
  
  
}

get.model.human.rho <- function(model.thresholds, human.thresolds) {
  bin.values <- get_experiment_bin_values()
  
  m <- model.thresholds %>%
    arrange(TARGET, BIN) %>%
    select(TARGET, BIN, threshold) %>%
    merge(., bin.values) %>%
    arrange(TARGET, statType, BIN)
  
  h <- human.thresholds %>%
    data.frame %>%
    group_by(TARGET, BIN) %>%
    summarize(threshold = mean(threshold)) %>%
    select(TARGET, BIN, threshold) %>%
    merge(., bin.values) %>%
    arrange(TARGET, statType, BIN)
  
  cor.val <- rbind(m %>% mutate(type = "model"), h %>% mutate(type = "human")) %>%
    spread(key = type, value = threshold) %>%
    group_by(TARGET, statType) %>%
    mutate(human = scale(human), model = scale(model)) %>%
    group_by(statType) %>%    
    summarize(cor(human, model))
  
  fig.dat <- rbind(m %>% mutate(type = "model"), h %>% mutate(type = "human")) %>%
    spread(key = type, value = threshold) %>%
    group_by(TARGET, statType) %>%
    mutate(human = scale(human), model = scale(model))
  
  plot.fig <- ggplot(data = fig.dat, aes(x = human, y = model, colour = TARGET)) + geom_point() + facet_grid(~statType)

  save(file = '~/Dropbox/Calen/Dropbox/threshold.scatter.pdf')
}
