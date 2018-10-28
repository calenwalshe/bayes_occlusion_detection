# SVM for Target Detection
library("e1071")
library(dplyr)
library(tidyr)
library(purrr)
library(purrrlyr)
library(ggplot2)


load('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_data/template_response.rdata')
load('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_data/model_wide.rdata')

template_response <- template_response %>%
  select(-one_of("L","C","S","statType", "statValue")) %>%
  unique()

gauss.params <- model_wide %>% select(eccentricity, TARGET, BIN, data_sd_vec_0, data_sd_vec_1)

response_sd <- gather(gauss.params, "TPRESENT", "response_sd", 4:5)
response_sd <- response_sd %>% mutate(TPRESENT = ifelse(TPRESENT == "data_sd_vec_0", "absent", "present"))


template.response.grouped <- template_response %>%
  spread(function_name, TRESP) %>%
  group_by(PYRAMIDLVL, BIN, TARGET, eccentricity, TPRESENT) %>%
  nest() 

template.response.sd.join <- template.response.grouped %>%
  left_join(., response_sd, by = c("eccentricity", "TARGET", "BIN", "TPRESENT")) %>%
  as_tibble()

template.response.sd.join$TPRESENT <- as.factor(template.response.sd.join$TPRESENT)

add_noise <- 0

noise.data <- template.response.sd.join %>% 
  by_row(., function(row) {
    data <- row$data[[1]]
    response_sd <- row$response_sd[[1]]
    
    n_obs <- nrow(data)
    data[, "edge_cos"] <- data[, "edge_cos"] + rnorm(n_obs, 0, sd = response_sd[1] * add_noise)
    data[, "mean_only"] <- data[, "mean_only"] + rnorm(n_obs, 0, sd = response_sd[2] * add_noise)
    data[, "pattern_only"] <- data[, "pattern_only"] + rnorm(n_obs, 0, sd = response_sd[3] * add_noise)
    
    data.frame(PATCHID = data$PATCHID, edge_cos = data$edge_cos, mean_only = data$mean_only, pattern_only = data$pattern_only)
  }, .to = "noise.dat") %>%
  unnest(noise.dat)


svm.results <- noise.data %>%
  group_by(eccentricity, BIN, TARGET, eccentricity) %>%
  nest() %>%
  mutate(model = map(data, function(data) {
    svm.untuned <- svm(TPRESENT ~ mean_only + pattern_only + edge_cos, data = data, probability = TRUE)
    #svm.tuned   <- tune.svm(TPRESENT ~ mean_only + edge_cos + pattern_only, data = data, kernel = "radial", cost=10^(-1:2), gamma = c(.5, 1, 2))
    #svm.final   <- svm(TPRESENT ~ mean_only + pattern_only + edge_cos, data = data, kernel = "radial" , gamma = svm.tuned$best.parameters$gamma, cost = svm.tuned$best.parameters$cost)
  }))


svm.test.performance <- svm.results %>%
  mutate(response = map2(model, data, function(model, data) { 
    features <- data %>% select(pattern_only, edge_cos, mean_only)
    response <- predict(model, features, probability = TRUE)
    labels <- attributes(response)$levels
    probabilities <- attributes(response)$probabilities
    data.frame(PATCHID = data$PATCHID, lik = list(log(probabilities)), model_response = response, correct_response = data$TPRESENT)
  })) %>%
  unnest(response) %>%
  mutate(LLR.svm = lik.present - lik.absent) %>%
  rename(TPRESENT = correct_response)

# Import human responses
human.responses <- export.responses()

# Format human responses
human.responses <- human.responses %>%
  filter(TRIAL != 1) %>%
  rename(eccentricity_human = ECCENTRICITY)

human.responses$TPRESENT <-
  ifelse(human.responses$TPRESENT == 1, "present", "absent")

# Merge dataframes
human.svm <-
  inner_join(human.responses,
             svm.test.performance,
             by = c("BIN", "PATCHID", "TARGET", "TPRESENT"))

bin_values <- get_experiment_bin_values()

# Interpolate model responses at human eccentricities
interpolated.svm <- human.svm %>% 
  group_by(TRIAL, TPRESENT, SUBJECT, BIN, TARGET, PATCHID, eccentricity_human) %>%
  nest() %>%
  mutate(svm.llr.interpolated = map2(eccentricity_human,data, function(x,y) {
    svm.llr.interpolated = approxExtrap(y$eccentricity, y$LLR.svm, xout = x)$y
  }
  )) %>%
  unnest(svm.llr.interpolated, .preserve = data)

         