#' Import human data from raw text.
get_human_responses <-
  function(path = '~/Dropbox/Calen/Work/occluding/detection_model/_data/exported/human_data.txt') {
    library(dplyr)
    library(tidyr)
    
    bin_values <- get_experiment_bin_values()
    
    human_data <- read.table(path, sep = "\t", header = T)
    
    human_data <- human_data %>% filter(!SUBJECT %in% c("jsa", "yhb"), TRIAL != 1) %>% 
      rename(BIN = BINS)
    
    #human_data <- human_data %>% group_by(SUBJECT, BIN, TARGET) %>%
    #  mutate(n_ecc = length(unique((ECCENTRICITY)))) %>% filter(!(SESSION ==
    #                                                                2 & n_ecc > 5))
    
    human_data$CORRECT <- ifelse(human_data$HIT == 1 | human_data$CORRECTREJECTION == 1, 1, 0)
    
    human_data$BIN <- factor(human_data$BIN)
    
    human_data <- merge(human_data, bin_values) %>% 
      arrange(SUBJECT, TARGET, BIN, statType, SESSION, ECCENTRICITY) %>%
      as_tibble()
    
    human_data %>% rename(SUBJECT = observer)
    return(human_data)
  }

#' Produce a data frame from the raw data that contains summary statistics.
get_human_detect <- function(human_data) {
  if (missing(human_data)) {
    error("Missing human data")
  }
  
  source('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_model/import_model.R')
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
    mutate(bias = dprime / 2 - qnorm(hit_adj)) %>%
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
  thresholds <- psychometric_parameters %>%
    rowwise() %>%
    mutate(threshold = ((d0 *
                           e0 ^ b) / 1 - e0 ^ b) ^ (1 / b),
           threshold_pc = (e0^b * ((2 * qnorm(.7))/d0)^-1 - e0^b)^(1/b))
  
  return(thresholds)
}

#' Plot thresholds measured from the human and ideal performance.
plot_publication_thresholds <-
  function(human.thresholds,
           model.thresholds,
           statIn = "Lvals") {
    
    library(ggthemes)
    human.thresholds <- human.thresholds %>%
      mutate(TARGET = factor(TARGET, levels = c("vertical", "horizontal", "bowtie", "spot")))
    
    model.thresholds <- model.thresholds %>%
      mutate(TARGET = factor(TARGET, levels = c("vertical", "horizontal", "bowtie", "spot")))
    
    human.threshold.1 <- human.thresholds %>%
      group_by(TARGET, BIN, SUBJECT) %>%
      dplyr::summarize(se = se, threshold = mean(threshold)) %>%
      select(TARGET, BIN, threshold, SUBJECT, se) %>%
      as_tibble()
    
    model.threshold.1 <- model.thresholds %>%
      select(TARGET, BIN, threshold, SUBJECT) %>%
      as_tibble()
    
    threshold_values <- full_join(model.threshold.1, human.threshold.1)
    
    threshold_values <- threshold_values %>% mutate(linetype = ifelse(SUBJECT %in% c("dprime.opt", "dprime.sub") ,1,2))
    
    bin_values <- get_experiment_bin_values()
    
    d.1 <- merge(threshold_values, bin_values) %>%
      group_by(TARGET, BIN, statValue, statType, SUBJECT) %>%
      dplyr::summarize(se = se, threshold = mean(threshold), linetype = linetype) %>%
      filter(statType == statIn) %>%
      arrange(TARGET, BIN, SUBJECT)
    
    d.1$linetype <- as.factor(d.1$linetype)
    t.1.plot <-
      ggplot(data = d.1, aes(x = statValue, y = threshold,
                             colour = SUBJECT, group = SUBJECT)) +
      geom_point(size = 2) +
      geom_line(size = 1.25) +
      facet_wrap(~ TARGET) +
      theme_set(theme_bw(base_size = 35))  +# pre-set the bw theme.
      theme(aspect.ratio = 1, legend.key.height=unit(3,"line")) +
      scale_color_brewer(name = "Subject", palette = "Dark2") +
      ylab("Eccentricity Threshold (º)")
    
    if (statIn == "Lvals") {
      t.1.plot <- t.1.plot + expand_limits(x = c(0,1), y = c(0, 23)) + xlab("Background Luminance (%)") +
        geom_errorbar(aes(ymin=threshold-se, ymax=threshold+se), width = 4)
    } else if (statIn == "Cvals") {
      t.1.plot <- t.1.plot + expand_limits(x = c(0, 1), y = c(0, 15)) + xlab("Background Contrast (RMS)") +
        geom_errorbar(aes(ymin=threshold-se, ymax=threshold+se), width = .1)
    } else if (statIn == "Svals") {
      t.1.plot <- t.1.plot + expand_limits(x = c(.4,1), y = c(0,18)) + xlab("Background Similarity") +
        geom_errorbar(aes(ymin=threshold-se, ymax=threshold+se), width = .025)
    }
    plot(t.1.plot)
    ggsave(t.1.plot, file = paste0('~/Dropbox/Calen/Work/occluding/detection_model_analysis/presentations/', statIn, "_thresholds.pdf"), scale= 1.35)
  }
