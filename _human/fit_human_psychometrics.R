# Get fitted psychometrics for all conditions separately.
get_fitted_psychometric_all <- function(human.responses) {
  library(bbmle)
  library(dplyr)
  library(multidplyr)
  library(purrr)
  library(curry)
  
  
  human.responses <- human.responses %>%
    select(
      SUBJECT,
      TARGET,
      BIN,
      SESSION,
      LEVEL,
      TRIAL,
      ECCENTRICITY,
      HIT,
      MISS,
      FALSEALARM,
      CORRECTREJECTION
    ) %>%
    distinct()
  
  # Setup paralell
  cluster <- create_cluster(cores = 16)
  set_default_cluster(cluster)
  cluster_assign_value(cluster, 'f.dprime.eccentricity', f.dprime.eccentricity)
  cluster_assign_value(cluster, 'f.NLL', f.NLL)
  cluster_library(cluster, 'purrr')
  cluster_library(cluster, 'curry')
  cluster_library(cluster, 'DEoptim')
  cluster_library(cluster, 'dplyr')
  
  # Setup wide data
  human.responses.wide <- human.responses %>%
    arrange(SUBJECT, TARGET, BIN, SESSION, LEVEL, TRIAL, ECCENTRICITY) %>%
    group_by(SUBJECT, TARGET, BIN, HIT, MISS, FALSEALARM, CORRECTREJECTION) %>%
    summarize(ECCENTRICITY = list(ECCENTRICITY)) %>%
    filter(!SUBJECT %in% c("jsa", "yhb")) %>%
    mutate(RESPONSE = ifelse(HIT == 1, "HIT", ifelse(
      MISS == 1,
      "MISS",
      ifelse(
        CORRECTREJECTION == 1,
        "CORRECTREJECTION",
        ifelse(FALSEALARM == 1, "FALSEALARM", "MISS")
      )
    ))) %>%
    ungroup() %>%
    select(SUBJECT, TARGET, BIN, ECCENTRICITY, RESPONSE)
  
  human.responses.nested <- human.responses.wide %>%
    group_by(SUBJECT, TARGET, BIN) %>%
    nest()
  
  # Fit all conditions separately. 
  fit.all <- human.responses.nested %>%
    partition(TARGET, cluster = cluster) %>%
    mutate(params = map(data, function(data) {
      d0    <- 5
      gamma <- 0
      
      HIT <- unlist(data[data$RESPONSE == "HIT", 1]$ECCENTRICITY)
      FA <-
        unlist(data[data$RESPONSE == "FALSEALARM", 1]$ECCENTRICITY)
      MISS <- unlist(data[data$RESPONSE == "MISS", 1]$ECCENTRICITY)
      CR <-
        unlist(data[data$RESPONSE == "CORRECTREJECTION", 1]$ECCENTRICITY)
      
      data.1 <- NULL
      data.1$HIT  <- HIT
      data.1$FA   <- FA
      data.1$CR   <- CR
      data.1$MISS <- MISS
      
      f <-
        partial(f.NLL, list(data = data.1)) # Partially apply the data frame to the objective function
      
      g <- function(x) {
        f(d0, x[1], x[2], x[3]) # Wrap the objective function in the format used by DEoptim.
      }
      
      starts <-
        DEoptim(
          g,
          lower = c(e0 = 0, b = 0, gamma = 0),
          upper = c(e0 = 20, b = 20, gamma = 0),
          control = list(
            trace = 6,
            reltol = .00001,
            steptol = 10
          )
        )
      
      return(starts$optim)
    })) %>%
    collect() %>%
    mutate(
      e0 = map(params, c("bestmem", "e0")),
      b = map(params, c("bestmem", "b")),
      gamma = map(params, c("bestmem", "gamma"))
    ) %>%
    unnest(e0, b, gamma)
  
  fit.all$d0 <- 5
  
  return(fit.all)
}

# Get nested optimization
get_fitted_nested <- function(fitted.psychometric) {
  library(curry)
  library(dplyr)
  library(parallel)
  
  get_fitted.b <- function(b, e) {
    f.list <- lapply(data, function(sub.data) {
      f.eval <- function(x) {
        HIT <- unlist(sub.data[sub.data$RESPONSE == "HIT", 1]$ECCENTRICITY)
        FA <-
          unlist(sub.data[sub.data$RESPONSE == "FALSEALARM", 1]$ECCENTRICITY)
        MISS <-
          unlist(sub.data[sub.data$RESPONSE == "MISS", 1]$ECCENTRICITY)
        CR <-
          unlist(sub.data[sub.data$RESPONSE == "CORRECTREJECTION", 1]$ECCENTRICITY)
        
        data <- list(
          HIT = HIT,
          FA = FA,
          MISS = MISS,
          CR = CR
        )
        
        f <- curry::partial(f.NLL, list(
          b = b,
          d0 = 5,
          gamma = 0, 
          data = data
        ))
        
        f(x[1])
        
      }
    })
    
    min.e0 <- min(fitted.psychometric$e0)
    max.e0 <- max(fitted.psychometric$e0)
    
    min.gamma <- min(fitted.psychometric$gamma)
    max.gamma <- max(fitted.psychometric$gamma)
    
    initial.pop <-
      as.matrix((data.frame(
        e0 = fitted.psychometric$e0)) %>% 
          sample_n(10))
    
    results.l2 <- mclapply(f.list, function(f.apply) {
      optim.result <- DEoptim(
        f.apply,
        lower = c(e0 = min.e0),
        upper = c(e0 = max.e0),
        control = list(
          reltol = .001,
          steptol = 5,
          trace = 0,
          NP = 10,
          initialpop = initial.pop
        )
      )$optim
      
      #list(optim.result = optim.result)
      return(optim.result)
    }, mc.cores = 16)
    
    NLL <-
      do.call(sum, map(results.l2, c("bestval")))
    
    l2.store <- list(b = b, optim = results.l2, NLL = NLL)
    
    storage <<- append(storage, list(result = l2.store))
    return(NLL)
  }  
  
  data <- fitted.psychometric$data
  
  storage <- list() # create storage container that is used to store results from nested optimization
  
  min.b <- min(fitted.psychometric$b)
  max.b <- max(fitted.psychometric$b)
  
  initial.b <- as.matrix(data.frame(b = fitted.psychometric$b) %>% sample_n(10)) # initial parameter for beta 

  results.l1 <-
    DEoptim(
      get_fitted.b,
      lower = c(b = min.b),
      upper = c(b = max.b),
      e = environment(),
      DEoptim.control(
        reltol = .001,
        steptol = 5,
        trace = 1,
        NP = 10,
        initialpop = initial.b
      )
    )

  NLL <- unlist(map(storage, "NLL"))
  best.model <- storage[[which(NLL == min(NLL))]]

  fit.nested    <- fitted.psychometric %>% 
    ungroup()
  
  fit.nested$d0 <- 5

  fit.nested$b <- best.model$b
  fit.nested[, c("e0")] <-
    do.call(rbind, map(best.model$optim, c("bestmem")))
  fit.nested$gamma <- 0
  
  return(fit.nested)
}

get_bootstrap_e0 <- function(fit.psychometric, n_boot = 100) {
  b     <- mean(fit.psychometric$b)
  gamma <- mean(fit.psychometric$gamma)
  d0    <- mean(fit.psychometric$d0)
  data  <- fit.psychometric$data
  
  data <- fit.psychometric$data
  #
  data.1 <- map(data, function(data) {
    d.1 <- unnest(data)
  })
  
  fit.psychometric$data <- data.1
  
  min.e0 <- min(fit.psychometric$e0)
  max.e0 <- max(fit.psychometric$e0)
  initial.e0 <- as.matrix(sample(fit.psychometric$e0, 10, replace = F))
  
  
  boot.samples <- mclapply(1:n_boot, function(x) by_row(fit.psychometric, function(row) {
    data <- row$data[[1]]
    
    data <- sample_frac(data, 1, replace = T)
    
    HIT  <- unlist(data[data$RESPONSE == "HIT", 2]$ECCENTRICITY)
    FA   <-  unlist(data[data$RESPONSE == "FALSEALARM", 2]$ECCENTRICITY)
    MISS <- unlist(data[data$RESPONSE == "MISS", 2]$ECCENTRICITY)
    CR   <- unlist(data[data$RESPONSE == "CORRECTREJECTION", 2]$ECCENTRICITY)
    
    data.1 <- NULL
    data.1$HIT  <- sample(HIT, replace = T)
    data.1$FA   <- sample(FA, replace = T)
    data.1$CR   <- sample(CR, replace = T)
    data.1$MISS <- sample(MISS, replace = T)
    
    f <- curry::partial(f.NLL, list(d0 = 5, b = b, data = data.1, gamma = 0))
    
    results.l1 <-
      DEoptim(
        f,
        lower = c(e0 = min.e0),
        upper = c(e0 = max.e0),
        DEoptim.control(
          reltol = .01,
          steptol = 5,
          trace = 1,
          NP = 10,
          initialpop = initial.e0
        )
      )
    
    return(results.l1$optim$bestmem)
  }, .to = "boot_e0"))
   
}

# Refit psychometric function fit with a varying beta with a beta now fixed by target.
get_fit_fixed_b <- function(human.psychometrics) {
  grouped.b <- human.psychometrics %>%
    group_by(SUBJECT) %>%
    mutate(b = mean(b))
  
  fitted.e0  <- grouped.b %>%
    group_by(SUBJECT, TARGET, BIN) %>%
    nest(-SUBJECT, -TARGET) %>%
    mutate(e0 = map(data, function(data) {
      sub.data <- data$data[[1]]
      
      HIT  <- unlist(sub.data[sub.data$RESPONSE == "HIT", 1]$ECCENTRICITY)
      FA   <-
        unlist(sub.data[sub.data$RESPONSE == "FALSEALARM", 1]$ECCENTRICITY)
      MISS <-
        unlist(sub.data[sub.data$RESPONSE == "MISS", 1]$ECCENTRICITY)
      CR   <-
        unlist(sub.data[sub.data$RESPONSE == "CORRECTREJECTION", 1]$ECCENTRICITY)
      
      responses <- list(
        HIT = HIT,
        FA = FA,
        MISS = MISS,
        CR = CR
      )
      
      b     <- data$b
      d0    <- data$d0
      gamma <- data$gamma

      f <- curry::partial(f.NLL, list(d0 = d0, b = b, gamma = gamma, data = responses))
      #f(1)
      optim.eval <- optimize(f, c(0,20))$min
    })) %>%
    unnest(e0) %>%
    select(SUBJECT, TARGET, BIN, e0)
    
  fitted.psychometric <- grouped.b %>%
    select(-e0) %>%
    merge(fitted.e0, ., by = c("SUBJECT", "TARGET", "BIN")) %>%
    as_tibble() %>%
    select(-data, -params) %>%
    arrange(SUBJECT, TARGET, BIN)
}

# Objective function
f.NLL <- function(d0, e0, b, gamma, data) {
  HIT <- data$HIT
  CR <- data$CR
  FA <- data$FA
  MISS <- data$MISS
  
  f <-
    Vectorize(partial(f.dprime.eccentricity, list(
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

