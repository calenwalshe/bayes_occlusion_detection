#' Produce a data frame from the raw data that contains summary statistics.
get_human_detect <- function(human_data) {
  if (missing(human_data)) {
    error("Missing human data")
  }
  
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R')
  experiment_bin_values <- get_experiment_bin_values()
  
  human_detect <- human_data %>%
    group_by(SUBJECT, ECCENTRICITY, BIN, TARGET) %>%
    dplyr::summarize(
      hit = sum(HIT == 1) / (sum(HIT == 1) + sum(MISS == 1)),
      miss = sum(MISS == 1) / (sum(MISS == 1) + sum(HIT == 1)),
      falsealarm = sum(FALSEALARM == 1) / (sum(FALSEALARM == 1) + sum(CORRECTREJECTION == 1)),
      correctrejection = sum(CORRECTREJECTION == 1) / (sum(CORRECTREJECTION == 1) + sum(FALSEALARM == 1)),
      percent_correct = mean(CORRECT)
    ) %>% mutate(
      hit_adj = ifelse(hit == 0, 1 / 1200, ifelse(hit == 1, 1199 / 1200, hit)),
      falsealarm_adj = ifelse(
        falsealarm == 0,
        1 / 1200,
        ifelse(falsealarm == 1, 1199 / 1200, falsealarm)
      ),
      percent_correct_adj = ifelse(percent_correct == 1, 2399 / 2400, percent_correct),
      percent_correct_adj = ifelse(percent_correct == 0, 1 / 2400, percent_correct)
    ) %>%
    mutate(dprime = qnorm(hit_adj) - qnorm(falsealarm_adj)) %>%
    mutate(bias = -1 / 2 * (qnorm(hit) + qnorm(falsealarm))) %>%
    merge(., experiment_bin_values, by = c("BIN", "TARGET")) %>%
    rename(eccentricity = ECCENTRICITY) %>%
    arrange(SUBJECT, BIN, TARGET, eccentricity, percent_correct) %>%
    data.frame()
  
  human_detect        <- human_detect %>%
    arrange(SUBJECT, statType, TARGET, BIN, eccentricity)
  
  return(human_detect)
}

# Compute summary statistics for the model from the raw responses
get_model_detect <- function(model_data) {
  experiment_bin_values <- get_experiment_bin_values()
  
  model_detect <- model_data %>%
    group_by(SUBJECT, ECCENTRICITY, BIN, TARGET) %>%
    dplyr::summarize(
      hit = sum(HIT == 1) / (sum(HIT == 1) + sum(MISS == 1)),
      miss = sum(MISS == 1) / (sum(MISS == 1) + sum(HIT == 1)),
      falsealarm = sum(FALSEALARM == 1) / (sum(FALSEALARM == 1) + sum(CORRECTREJECTION == 1)),
      correctrejection = sum(CORRECTREJECTION == 1) / (sum(CORRECTREJECTION == 1) + sum(FALSEALARM == 1)),
      percent_correct = mean(CORRECT)
    ) %>% mutate(
      hit_adj = ifelse(hit == 0, 1 / 1200, ifelse(hit == 1, 1199 / 1200, hit)),
      falsealarm_adj = ifelse(
        falsealarm == 0,
        1 / 1200,
        ifelse(falsealarm == 1, 1199 / 1200, falsealarm)
      ),
      percent_correct_adj = ifelse(percent_correct == 1, 2399 / 2400, percent_correct),
      percent_correct_adj = ifelse(percent_correct == 0, 1 / 2400, percent_correct)
    ) %>%
    mutate(dprime = (qnorm(hit_adj) - qnorm(falsealarm_adj))) %>%
    mutate(bias = dprime / 2 - qnorm(hit_adj)) %>%
    merge(., experiment_bin_values, by = c("BIN", "TARGET")) %>%
    rename(eccentricity = ECCENTRICITY) %>%
    arrange(SUBJECT, BIN, TARGET, eccentricity, percent_correct) %>%
    data.frame()
  
  model_detect$TARGET <-
    factor(
      model_detect$TARGET,
      levels = c("vertical", "horizontal", "bowtie", "spot"),
      ordered = T
    )
  
  model_detect        <- model_detect %>%
    arrange(SUBJECT, statType, TARGET, BIN, eccentricity)
  
  return(model_detect)
}

#' Add thresholds to a dataframe that contains psychometric parameters
get_threshold <- function(psychometric_parameters) {
  approxExtrap <- Hmisc::approxExtrap
  thresholds <- psychometric_parameters %>%
    rowwise() %>%
    mutate(threshold = ((d0 *
                           e0 ^ b) / 1 - e0 ^ b) ^ (1 / b),
           threshold_pc = (e0^b * ((2 * qnorm(.7))/d0)^-1 - e0^b)^(1/b)
           )
  return(thresholds)
}

#' Compute eccentricity threshold from linearly interpolated model responses
#'
#' @param psychometric 
#' @param eccentricity 
#'
#' @return
#' @export
#'
#' @examples
get_interpolated_threshold <- function(model.psychometric, human.psychometrics) {
  library(tidyr)
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_human/human_psychometrics_tools.R')
  
  #model.error.1 <- model.psychometric %>%
  #  mutate(dprime = dprime * scale) %>%
  #  #select(BIN, TARGET, eccentricity, observer, dprime) %>%
  #  group_by(BIN, TARGET, observer) %>% 
  #  filter(observer == "optimal") %>%
  #  nest() 
  
  model.psychometric$threshold <- map(model.psychometric$data, function(x) {
    ln.mod <- approxExtrap(x$dprime, x$eccentricity, 1)$y
  }) %>% unlist()
  
  model.psychometric <- model.psychometric %>% select(-data)
  
  return(model.psychometric)
}

#' Add thresholds to a dataframe that contains psychometric parameters
get_threshold_extrap <- function(psychometric_parameters) {
  approxExtrap <- Hmisc::approxExtrap
  
  thresholds <- by_row(psychometric_parameters, function(x) {
    dprime.dat <- x$data[[1]]
    
    bExtrap <- !((max(dprime.dat$dprime) > 1) & min((dprime.dat$dprime) < 1))

    if (bExtrap) {
      threshold <- approxExtrap(dprime.dat$dprime, dprime.dat$eccentricity, 1)$y
    } else {
      threshold <- ((x$d0 * x$e0 ^ x$b) / 1 - x$e0 ^ x$b) ^ (1 / x$b)
    }
    
    return(threshold)
  }, .to = "threshold") %>% unnest(threshold)
  
  return(thresholds)
}

#' Bootstrap thresholds from fitted psychometric function
#'
#' @return
#' @export
#'
#' @examples
bootstrap.psychometrics <- function(human.psychometrics, n_samples) {
  library(bbmle)
  library(broom)
  library(parallel)
  library(R.matlab)
  library(tidyr)
  library(bbmle)
  
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_human/export_responses.R')
  load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/human.psychometrics.rdata")
  
  raw.data <- export.responses()
  
  raw.data$observer <- raw.data$SUBJECT
  raw.data$SUBJECT  <- NULL
  
  f <- function() {
    data.grouped <- raw.data %>%
      rename(eccentricity = ECCENTRICITY) %>%
      group_by(TARGET, observer, BIN) %>%
      sample_frac(1, replace = T) %>%
      filter(TRIAL != 1) %>%
      nest()

    data.grouped <- left_join(data.grouped, human.psychometrics %>% select(TARGET, observer, BIN, b, d0), by = c("TARGET", "observer", "BIN"))
    
    data.grouped <- data.grouped %>% 
      group_by(TARGET, observer, BIN) %>% 
      nest()
    
    many.models.fixed.b <- map(data.grouped$data, function(x) {
      data          <- x$data[[1]]
      eccentricity  <- data$eccentricity
      HIT           <- data$HIT
      CORRECTREJECTION <- data$CORRECTREJECTION
      
      b.fix  <- x$b
      d0.fix <- x$d0
      
      f.mle <- function(x) {
        e0 <- x[1]
        sample <- (1 / 2 * (-1 + (HIT |
                                    CORRECTREJECTION) * 2) * d0.fix * e0 ^ b.fix / (e0 ^ b.fix + eccentricity ^ b.fix))
        
        -sum(pnorm(sample, mean = 0, sd = 1, log.p = T))
      }
      
      optim(par = c(e0 = 5), f.mle, method = "Brent", lower = 0, upper = 15)
    })
    
    data.grouped$parameters = map(many.models.fixed.b, function(x) data.frame(as.list(c(e0 = x$par))))
    
    data.grouped <- data.grouped %>% unnest()
    
    data.grouped$data <- map(data.grouped$data, function(x) x %>% group_by(eccentricity) %>% summarize(pc = mean(CORRECT), dprime = qnorm(mean(HIT)) - qnorm(mean(FALSEALARM))))
    
    data.return <- data.grouped %>%
      get_threshold(.)
    
    return(data.return)
  }
  
  all.boots <- mclapply(1:n_samples, FUN = function(x) f(), mc.cores = 8)
  
  tt <- map(all.boots, "threshold") %>% do.call(rbind,.)
  
  ttt <- apply(tt, 2, sd)
  
  human.psychometrics$sd <- ttt
  
  return(human.psychometrics)
}

#' Visualize a single psychometric function
plot_psychometric <- function(psychometrics, out_path = "~/Dropbox/Calen/Dropbox/") {
  library(dplyr)
  library(ggplot2)
  library(purrrlyr)
  library(purrr)
  
  psychometrics <- psychometrics %>% arrange(observer, BIN, TARGET)
  
  n_row <- nrow(psychometrics)
  
  plot.fig.data <- lapply(1:n_row, FUN = function(x) {
    data <- psychometrics[x, ]
    # model fit
    eccentricity <- seq(0, 25, .1)
    d0 <- data$d0
    e0 <- data$e0
    b  <- data$b
    gamma  <- data$gamma
    dprime <- d0 * e0^b/(e0^b + eccentricity^b)
    pc     <- (pnorm(dprime/2 - gamma) + (1 - pnorm(-dprime/2 - gamma)))/2
    
    model <- data.frame(TARGET = data$TARGET, BIN = data$BIN, observer = data$observer, eccentricity = eccentricity, pc = pc, dprime = dprime, data = factor("psychometric"))
    
    
    # raw data
    raw.data <- psychometrics[[x, "data"]]
    eccentricity <- raw.data$eccentricity
    dprime       <- raw.data$dprime
    pc           <- raw.data$pc
    
    response.dat <- data.frame(TARGET = data$TARGET, BIN = data$BIN, observer = data$observer, eccentricity = eccentricity, pc = pc, dprime = dprime, data = factor("response"))
    
    all.dat <- rbind(model, response.dat)
    
  })
  
  plot.data <- do.call(rbind, plot.fig.data) %>% as_tibble()
  
  fig <- ggplot() + 
    geom_point(data = plot.data %>% filter(data == "response"), 
               aes(x = eccentricity, y = pc, colour = observer)) +
    geom_line(data = plot.data %>% filter(data == "psychometric"), 
              aes(x = eccentricity, y = pc, colour = observer)) +
    facet_grid(TARGET ~ BIN) +
    theme(aspect.ratio = 1)
  
  ggsave(filename = '~/Dropbox/Calen/Dropbox/psychometrics.pdf', plot = fig, device = 'pdf', width = 40, height = 30, units = "in")
  
}

#' Plot thresholds measured from the human and ideal performance.
plot_publication_thresholds <-
  function(human.thresholds,
           model.thresholds,
           statIn = "Lvals", name.add = "1") {
    
    library(ggthemes)
    
    human.thresholds <- human.thresholds %>% ungroup() %>%
      mutate(TARGET = factor(TARGET, levels = c("vertical", "horizontal", "bowtie", "spot"))) %>% group_by(TARGET, BIN) %>% summarize(threshold = mean(threshold)) %>% mutate(observer = "ave", se = 0)
    
    model.thresholds <- model.thresholds %>%
      mutate(TARGET = factor(TARGET, levels = c("vertical", "horizontal", "bowtie", "spot")))
    
    human.threshold.1 <- human.thresholds %>%
      group_by(TARGET, BIN, observer) %>%
      dplyr::summarize(se = mean(se), threshold = mean(threshold)) %>%
      select(TARGET, BIN, threshold, observer, se) %>%
      as_tibble()
    
    model.threshold.1 <- model.thresholds %>%
      select(TARGET, BIN, threshold, observer,se) %>%
      as_tibble()
    
    threshold_values <- full_join(model.threshold.1, human.threshold.1)
    
    bin_values <- get_experiment_bin_values()
    
    d.1 <- merge(threshold_values, bin_values) %>%
      group_by(TARGET, BIN, statValue, statType, observer) %>%
      dplyr::summarize(se = mean(se), threshold = mean(threshold)) %>%
      filter(statType == statIn) %>%
      arrange(TARGET, BIN, observer)
    
    d.1$linetype <- as.factor(ifelse(d.1$observer %in% c("rcw", "sps", "yhb", "ave"), "1", "2"))
    
    d.1 <- d.1 %>% filter(observer %in% c("ave", "rcw", "sps", "yhb", "1 0 0", "0 1 0", "0 0 1", "optimal"))
    
    label.set <- data.frame(names = c("ave", "rcw", "sps", "yhb", "optimal", "mahalanobis", "nocov", "1 0 0", "0 1 0", "0 0 1"), pretty.names = c("Average Human", "rcw", "sps", "yhb", "Optimal", "Mahalanobis", "No Covariance", "Edge", "Luminance", "Pattern"))
    
    my.labels <- (label.set[label.set$names %in% d.1$observer,])
    
    d.1$observer <- factor(d.1$observer, levels = my.labels[,1], labels = my.labels[,2])
    
    t.1.plot <-
      ggplot(data = d.1, aes(x = statValue, y = threshold,
                             colour = observer, linetype = linetype)) +
      geom_point(size = 3) +
      geom_line(size = 1.25) +
      facet_wrap(~ TARGET) +
      theme(aspect.ratio = 1) +
      coord_cartesian(ylim = c(0, 60)) +
      ylab("Eccentricity Threshold (\U00B0)")
    
    if (statIn == "Lvals") {
      t.1.plot <- t.1.plot + expand_limits(x = c(0,1), y = c(0, 23)) + xlab("Background Luminance (%)")
    } else if (statIn == "Cvals") {
      t.1.plot <- t.1.plot + expand_limits(x = c(0, 1), y = c(0, 15)) + xlab("Background Contrast (RMS)")
    } else if (statIn == "Svals") {
      t.1.plot <- t.1.plot + expand_limits(x = c(.4,1), y = c(0,18)) + xlab("Background Similarity")
    }
    plot(t.1.plot)
    ggsave(t.1.plot, file = paste0('~/Dropbox/Calen/Dropbox/', statIn, '_', name.add, "_thresholds.pdf"), scale= 3.5)
    
    return(t.1.plot)
  }
