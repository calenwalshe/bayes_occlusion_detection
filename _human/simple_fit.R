simple.fit <- function(target = "vertical") {

  human.responses <- get_human_responses()
  human.responses <- human.responses %>% select(-L,-C, -S, -statType, -statValue) %>% distinct() %>%
    filter(SUBJECT != "yhb")

  formatted.responses <- human.responses %>%
    group_by(TARGET, BIN) %>%
    nest() %>%
    mutate(data = map(data, function(data) {format.response(data)})) %>%
    arrange(TARGET, BIN) %>%
    filter(TARGET == target)
  
    fit.separate <- formatted.responses %>%
    mutate(full.model = map(data, function(data) {
      f        <- curry::partial(f.NLL, list(d0 = 4.5, gamma = 0, data = data))
      class(f) <- "function" #stupid bbmle doesn't like the scaffold.

      mle2(f, start = list(e0 = 5, b = 4), 
           lower = c(e0 = .5, b = 0),
           upper = c(e0 = 40, b = 40),
           method = "L-BFGS-B") # let's party bbmle
  })) %>%
  mutate(d0 = 4.5, gamma = 0) %>%
  as_tibble()
  
  formatted.responses <- as_tibble(cbind(formatted.responses, data.frame(do.call(rbind, (map(fit.separate$full.model, coef))))))
  #cluster_copy(cluster, formatted.responses)
  
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
      f     <- x$f
      e0    <- x$e0
      class(f) <- "function"
      mle2(f, start = list(e0 = e0), method = "L-BFGS-B", lower = c(e0 = .5), upper = c(e0 = 23),
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
    FA                <- response.df$ECCENTRICITY[response.df$MISS == 1]
    CR                <- response.df$ECCENTRICITY[response.df$CORRECTREJECTION == 1]
    MISS              <- response.df$ECCENTRICITY[response.df$MISS == 1]
    
    return(list(HIT = HIT, FA = FA, CR = CR, MISS = MISS))
  }
  
  # Objective function
  f.NLL <- function(d0, e0, b, gamma, data) {
    HIT <- data$HIT
    CR <- data$CR
    FA <- data$FA
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
  
  optim.b_g <- mle2(b.NLL, start = list(b = mean(formatted.responses$b)), method="L-BFGS-B", lower = c(b = .5), upper = c(b = 10), control = list(trace = 4))
  
  fitted.psychometrics <- b.NLL(coef(optim.b_g), FALSE) %>% mutate(gamma = 0) %>% get_threshold(.)
  return(fitted.psychometrics)
}

fits <- lapply(list('vertical','horizontal', 'bowtie', 'spot'), FUN = function(x) simple.fit(x))

fits <- do.call(rbind, fits)

save(file = '~/Dropbox/Calen/Dropbox/fits.rdata', fits)