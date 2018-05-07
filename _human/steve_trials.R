steve.data <- function() {
  
  steve.data <- steve.detect %>%
    filter(SUBJECT == "sps") %>%
    select(-L,-C,-S,-statType, -statValue) %>%
    arrange(BIN, TARGET, SUBJECT, SESSION) %>%
    group_by(BIN, TARGET, SUBJECT) %>%
    filter(length(unique(eccentricity)) > 5) %>%
    nest()
  
  by_row(steve.data, ..f = function(row) {
    condition <- row$data[[1]]
    fig <- ggplot(condition, aes(x = eccentricity, y = dprime, colour = as.factor(SESSION))) + geom_point() + ggtitle(paste0(row$TARGET, row$BIN)) + theme(legend.title=element_blank())
    
    ggsave(filename = paste0('~/Dropbox/Calen/Dropbox/', ' BIN ', head(row$BIN,1), ' TARGET ', head(row$TARGET,1), '.pdf'), plot = fig, device = 'pdf')
  }
  )
  
  get_steve_detect <- function(human_data) {
    if (missing(human_data)) {
      error("Missing human data")
    }
    experiment_bin_values <- get_experiment_bin_values()
    
    human_detect <- human_data %>%
      group_by(SUBJECT, SESSION, ECCENTRICITY, BIN, TARGET) %>%
      summarize(
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
    
    human_detect$TARGET <-
      factor(
        human_detect$TARGET,
        levels = c("vertical", "horizontal", "bowtie", "spot"),
        ordered = T
      )
    
    human_detect        <- human_detect %>%
      arrange(SUBJECT, statType, TARGET, BIN, eccentricity)
    
    return(human_detect)
  }
  #' Import human data from raw text.
  get_steve_responses <-
    function(path = '~/Dropbox/Calen/Work/occluding/detection_model/_data/exported/human_data.txt') {
      library(dplyr)
      library(tidyr)
      
      bin_values <- get_experiment_bin_values()
      
      human_data <- read.table(path, sep = "\t", header = T)
      
      human_data <- human_data %>% filter(!SUBJECT %in% c("jsa",
                                                          "yhb"), TRIAL != 1) %>%
        rename(BIN = BINS)
      
      human_data <- human_data %>% group_by(SUBJECT, BIN, TARGET) %>%
        mutate(n_ecc = length(unique((ECCENTRICITY))))
      
      human_data$CORRECT <-
        ifelse(human_data$HIT == 1 | human_data$CORRECTREJECTION ==
                 1, 1, 0)
      
      human_data$BIN <- factor(human_data$BIN)
      
      human_data <-
        merge(human_data, bin_values) %>% arrange(SUBJECT,
                                                  TARGET, BIN, statType, SESSION, ECCENTRICITY)
      
      return(human_data)
    }
}