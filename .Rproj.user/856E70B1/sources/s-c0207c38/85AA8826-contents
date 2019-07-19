fit.target <- function(subject = "rcw") {

  library(dplyr)
  library(purrr)
  library(purrrlyr)
  
  n.cores <- detectCores(logical = FALSE) - 1
  
  human.responses <- get_human_responses()
  human.responses <- human.responses %>% select(-L,-C, -S, -statType, -statValue) %>% distinct()

  formatted.responses <- human.responses %>%
    filter(SUBJECT == subject) %>%
    group_by(TARGET, SUBJECT, BIN) %>%
    nest() %>%
    mutate(data = map(data, function(data) {format.response(data)})) %>%
    arrange(TARGET, BIN)
  
    fit.separate <- formatted.responses %>%
    mutate(full.model = map(data, function(data) {
      f        <- curry::partial(f.NLL, list(d0 = 4.5, gamma = 0, data = data))
      class(f) <- "function" #stupid bbmle doesn't like the scaffold.

      mle2(f, start = list(e0 = 5, b = 4), lower = c(e0 = .5, b = .5), upper = c(e0 = 23, b = 10), 
           method = "L-BFGS-B") # let's party bbmle
  })) %>%
  mutate(d0 = 4.5, gamma = 0) %>%
  as_tibble()
  
  formatted.responses <- as_tibble(cbind(formatted.responses, data.frame(do.call(rbind, (map(fit.separate$full.model, coef)))))) %>%
    mutate(d0 = 4.5, gamma = 0)
  
  # Compute the negative log likelihood for individual level e0 nested within target b0 and overall d0
  b.NLL <- function(b, bLikelihood = TRUE) {
    
    d0    <- 4.5
    b     <- b
    #b_v  <- params[2]
    #b_h  <- params[3]
    #b_b  <- params[4]
    #b_s  <- params[5]
    
    formatted.responses$d0 <- NULL
    formatted.responses$b  <- NULL
    
    formatted.responses$d0 <- d0
    formatted.responses$b  <- b
    
    # Use currying to fix the beta for the beta level fit.
    f.eval.list <- by_row(formatted.responses, function(data) {
      f <- list(f = curry::partial(f.NLL, list(d0 = data$d0, b = data$b, gamma = 0, data = data$data[[1]])), e0 = data$e0)},
      .to = 'opt.list')
    
    min.vals <- mclapply(f.eval.list$opt.list, FUN = function(x) {
      f        <- x$f
      e0       <- x$e0
      class(f) <- "function"
      mle2(f, start = list(e0 = e0), lower = c(e0 = .5), upper = c(e0 = 23), method = "L-BFGS-B",
           control = list(trace = 0))
    }, mc.cores=n.cores)
  
    if (bLikelihood == TRUE) {
      NLL <- Reduce('+', map(min.vals, function(x) {x@min}))
      return(NLL)
    } else
      
      fitted.parameters         <- map(min.vals, coef)
    
      formatted.responses$e0    <- unlist(map(fitted.parameters, 1))
      formatted.responses$gamma <- unlist(map
                                          (fitted.parameters, 2))
      
      return(formatted.responses)
  }
  
  # Format responses for the objective function
  format.response <- function(response.df) {
    HIT               <- response.df$ECCENTRICITY[response.df$HIT == 1]
    FA                <- response.df$ECCENTRICITY[response.df$FALSEALARM == 1]
    CR                <- response.df$ECCENTRICITY[response.df$CORRECTREJECTION == 1]
    MISS              <- response.df$ECCENTRICITY[response.df$MISS == 1]
    
    return(list(HIT = HIT, FA = FA, CR = CR, MISS = MISS))
  }
  
  # Objective function
  f.NLL <- function(d0, e0, b, gamma, data) {
    HIT  <- data$HIT
    CR   <- data$CR
    FA   <- data$FA
    MISS <- data$MISS
    
    f <-
      Vectorize(curry::partial(f.dprime.eccentricity, list(
        d0 = d0,
        e0 = e0,
        b = b,
        gamma = gamma
      )))
    
    nll.cr = sum(pnorm(1 / 2 * f(CR), log = T))
    nll.hit = sum(pnorm(1 / 2 * f(HIT), log = T))
    nll.miss = sum(pnorm(-1 / 2 * f(MISS), log = T))
    nll.fa = sum(pnorm(-1 / 2 * f(FA), log = T))

    return(-sum(nll.cr, nll.hit, nll.miss, nll.fa))
  }
  
  # Dprime as a function of eccentricity.
  f.dprime.eccentricity <- function(x, d0, e0, b, gamma = 0) {
    dprime <- d0 * e0 ^ b / (e0 ^ b + x ^ b) - gamma
  }

  optim.b_g <- mle2(b.NLL, start = list(b = mean(formatted.responses$b)), lower = c(b = .5), upper = c(b = 10), method = "L-BFGS-B", control = list(trace = 4))
  
  fitted.psychometrics <- b.NLL(coef(optim.b_g), FALSE) %>%
    mutate(gamma = 0) %>%
    get_threshold(.)
  
  return(fitted.psychometrics)
}

simple.fits <- function() {
  conditions <- c("rcw", "sps")
  
  fits <- lapply(conditions, FUN = function(x) {fit.target(x[1])})
  
  fitted.frame <- do.call(rbind, fits$.out)  
}

get_bootstrap <- function(fit.psychometric, n_boot = 100) {
  
  human.responses <- get_human_responses()
  
  human.responses <- human.responses %>% # import responses
    select(-L,-C, -S, -statType, -statValue) %>%
    distinct() %>%
    group_by(TARGET, BIN, SUBJECT) %>%
    nest() 
  
  boot.dat <- fit.psychometric %>% # setup data frame for boostrapping samples.
    select(-data) %>%
    merge(., human.responses) %>%
    as_tibble()
  
  
  boot.fits <- boot.dat %>% # run the bootstrap procedure
    group_by(TARGET, BIN, SUBJECT) %>%
    nest() %>%
    mutate(boot.fit = map(data, function(data) { # mutate and map does the bootstrap on each row separately
      
      e0               <- data$e0 # setup parameters for bootstrap. we are only bootstrapping the e0 parameter
      b                <- data$b
      psychometric.dat <- data$data[[1]]
      d0               <- data$d0
      
      fits <- mclapply(1:n_boot, FUN = function(x) { # using lapply to do the resampling.
        
        formatted.dat <- psychometric.dat %>% 
          sample_frac(1, replace = T) %>% format.response(.) # resampling
        
        f.obj        <- curry::partial(f.NLL, list(d0 = 4.5, b = b, gamma = 0, data = formatted.dat)) # bind the resampled data to the objective fcn.
        class(f.obj) <- "function"
        
        fits <- mle2(f.obj, start = list(e0 = e0), method = "L-BFGS-B", lower = c(e0 = .5), upper = c(e0 = 23),
           control = list(trace = 0)) # estimate the e0 for the resampled data
        
        e0 <- coef(fits)
        
        threshold <- ((d0 * e0 ^ b) / 1 - e0 ^ b) ^ (1 / b)
      }, mc.cores = 15)
      
      unlist(fits)
    }))
  
  boot.fits <- boot.fits %>%
    mutate(se = map(boot.fits$boot.fit, sd)) %>%
    unnest(se)
}
