#' Visualize a single psychometric function
plot_single_psychometric <- function(psychometric_parameters, 
                                     performance_measures, bin, target, out_path = "~/Dropbox/Calen/Dropbox/tmp_images/") {
  sub_params <- psychometric_parameters %>% filter(BIN == bin, 
                                                   TARGET == target) %>% data.frame()
  
  sub_response <- performance_measures %>% filter(BIN == bin, 
                                                  TARGET == target) %>% data.frame()
  
  responses <- lapply(1:nrow(sub_params), FUN = function(x) data.frame(SUBJECT = sub_params[x, 
                                                                                            "SUBJECT"], eccentricity = seq(0, 30, 0.1), pc = pnorm(sub_params[x, 
                                                                                                                                                              "d0"]/2 * (sub_params[x, "e0"]^sub_params[x, "b"])/(sub_params[x, 
                                                                                                                                                                                                                             "e0"]^sub_params[x, "b"] + seq(0, 30, 0.1)^sub_params[x, 
                                                                                                                                                                                                                                                                                   "b"])))) %>% do.call(rbind, .)
  
  responses <- responses %>% mutate(BIN = bin, TARGET = target)
  
  p <- ggplot(responses, aes(x = eccentricity, y = pc, colour = SUBJECT)) + 
    geom_line() + geom_point(data = sub_response, aes(x = eccentricity, 
                                                      y = percent_correct)) + scale_y_continuous(limits = c(0.5, 
                                                                                                            1)) + ggtitle(paste0(target, " ", bin))
  
  ggsave(last_plot(), file = paste0(out_path, bin, "-", target, 
                                    ".pdf"))
}

plot_eccentricity <- function(psychometric_parameters, performance_measures, 
                              measure = "pc", out_path = "~/Dropbox/Calen/Dropbox/") {
  library(ggplot2)
  library(dplyr)
  
  sub_params <- psychometric_parameters
  
  sub_response <- performance_measures
  
  sub_params <- merge(sub_params, data.frame(eccentricity = unique(sub_response$eccentricity)))
  
  pc <- lapply(1:nrow(sub_params), FUN = function(x) pc = pnorm(sub_params[x, 
                                                                           "d0"]/2 * (sub_params[x, "e0"]^sub_params[x, "b"])/(sub_params[x, 
                                                                                                                                          "e0"]^sub_params[x, "b"] + sub_params[x, "eccentricity"]^sub_params[x, 
                                                                                                                                                                                                              "b"]))) %>% do.call(rbind, .)
  
  sub_params$pc <- pc
  
  human_dat <- sub_params %>% group_by(TARGET, BIN, eccentricity) %>% 
    summarize(pc = mean(pc)) %>% ungroup() %>% mutate(dprime = 2 * 
                                                        as.numeric(qnorm(pc)), SUBJECT = "human") %>% data.frame
  
  model_dat <- performance_measures %>% select(TARGET, BIN, 
                                               eccentricity, percent_correct, dprime) %>% rowwise() %>% 
    mutate(pc = percent_correct, dprime = 2 * as.numeric(qnorm(pc)), 
           percent_correct = NULL, SUBJECT = "model") %>% data.frame
  
  dat <- rbind(human_dat, model_dat)
  
  all.fig <- ggplot(dat, aes_string(x = "BIN", y = measure, 
                                    group = "TARGET", colour = "SUBJECT")) + geom_point() + 
    facet_grid(~eccentricity)
  
  ggsave(last_plot(), file = paste0(out_path, "all_fig", ".pdf"), 
         width = 50, height = 50/(1920/1080), units = c("cm"))
}

# Plot thresholds derived from the psychometric function.
plot_thresholds <- function(threshold_values, out_path = "~/Dropbox/Calen/Dropbox/tmp_images/") {
  bin_values <- get_experiment_bin_values()
  
  d.1 <- merge(threshold_values, bin_values)
  
  t.1.plot <- ggplot(data = d.1, aes(x = statValue, y = threshold, 
                                     colour = SUBJECT)) + geom_point() + facet_wrap(statType ~ 
                                                                                      TARGET, scales = "free_x") + geom_line()
  
  ggsave(last_plot(), file = paste0(out_path, "thresholds.pdf"))
}

#' Plot the fitted psychometric functions and the psychometric data together.
plot_human_psychometrics <- function(params, response_data, out_path, 
                                     target = "vertical") {
  
  params <- params %>% filter(TARGET %in% target)
  response_data <- response_data %>% filter(TARGET %in% target)
  
  interpolated_responses <- generate_epf_observations(params, 
                                                      gamma_bias = TRUE)
  
  interpolated_responses <- interpolated_responses %>% filter(psychometric_type %in% 
                                                                c("percent_correct", "no_bias")) %>% mutate(has_bias = ifelse(gamma_bias == 
                                                                                                                                0, FALSE, TRUE)) %>% mutate(percent_correct = PC)
  
  interpolated_responses <- as.data.frame(lapply(interpolated_responses,
                                                 function(x) if (is.factor(x)) 
                                                   factor(x) else x)) %>% filter(has_bias == FALSE)
  
  response_data <- as.data.frame(lapply(response_data, function(x) if (is.factor(x)) 
    factor(x) else x))
  
  interpolated_responses$BIN <- factor((interpolated_responses$BIN), 
                                       levels = seq(1, 15, 1))
  
  psychometrics_percentcorrect <- ggplot(interpolated_responses %>% filter(statType == "Lvals"), 
                                         aes(x = eccentricity, y = percent_correct, group = SUBJECT)) + 
    geom_line() + 
    geom_point(data = response_data %>% filter(statType == "Lvals"), aes(x = eccentricity, y = percent_correct)) +
    facet_wrap(~statValue, scales = "free") +
    scale_x_continuous(limits = c(0, 30)) +
    scale_y_continuous(limits = c(0.25, 1)) +
    theme(axis.line = element_line()) + 
    ggtitle(target)
  
  ggsave("~/Dropbox/Calen/Dropbox/lvals.pdf", height = 40, width = 40, units = "in", device = "pdf")    
  
  psychometrics_percentcorrect <- ggplot(interpolated_responses %>% filter(statType == "Cvals"), 
                                         aes(x = eccentricity, y = percent_correct, group = SUBJECT)) + 
    geom_line() + 
    geom_point(data = response_data %>% filter(statType == "Cvals"), aes(x = eccentricity, y = percent_correct)) +
    facet_wrap(~statValue, scales = "free") +
    scale_x_continuous(limits = c(0, 30)) +
    scale_y_continuous(limits = c(0.25, 1)) +
    theme(axis.line = element_line()) + 
    ggtitle(target)
  
  ggsave("~/Dropbox/Calen/Dropbox/cvals.pdf", height = 40, width = 40, units = "in", device = "pdf")        
  
  psychometrics_percentcorrect <- ggplot(interpolated_responses %>% filter(statType == "Svals"), 
                                         aes(x = eccentricity, y = percent_correct, group = SUBJECT)) + 
    geom_line() + 
    geom_point(data = response_data %>% filter(statType == "Svals"), aes(x = eccentricity, y = percent_correct)) +
    facet_wrap(~statValue, scales = "free") +
    scale_x_continuous(limits = c(0, 30)) +
    scale_y_continuous(limits = c(0.25, 1)) +
    theme(axis.line = element_line()) + 
    ggtitle(target)    
  
  ggsave("~/Dropbox/Calen/Dropbox/svals.pdf", height = 40, width = 40, units = "in", device = "pdf")            
}

#' Returns interpolated observations for the eccentricity psychometric functions stored in params
generate_epf_observations <- function(params, gamma_bias = FALSE) {
  if (missing(params)) {
    error("Parameters missing")
  }
  
  x <- seq(0, 30, 0.01)
  interpolated_data <- data.frame(BIN = NULL, TARGET = NULL, 
                                  SUBJECT = NULL, eccentricity = NULL, PC = NULL, gamma_bias = NULL, 
                                  psychometric_type = NULL)
  
  for (i in unique(params$BIN)) {
    for (j in unique(params$TARGET)) {
      for (k in unique(params$SUBJECT)) {
        sub_params <- params %>% 
          filter(BIN == i, TARGET == j, SUBJECT == k) %>% .[1, ]
        
        d0 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                         eccentricity = x, PC = with(sub_params, pnorm(0.5 * 
                                                                         dmax * (e0^b)/(e0^b + x^b), mean = 0, sd = 1)), 
                         gamma_bias = 0, psychometric_type = "no_bias")
        
        d1 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                         eccentricity = x, PC = with(sub_params, pnorm(0.5 * 
                                                                         dmax * (e0^b)/(e0^b + x^b) - g, 
                                                                       mean = 0, sd = 1)), gamma_bias = sub_params$g, 
                         psychometric_type = "hit")
        
        d2 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                         eccentricity = x, PC = with(sub_params, pnorm(-0.5 * 
                                                                         dmax * (e0^b)/(e0^b + x^b) - g, 
                                                                       mean = 0, sd = 1)), gamma_bias = sub_params$g, 
                         psychometric_type = "falsealarm")
        
        d3 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                         eccentricity = x, PC = with(sub_params, 0.5 * 
                                                       pnorm(0.5 * dmax * (e0^b)/(e0^b + x^b) - 
                                                               sub_params$g) + 0.5 * (1 - pnorm(-0.5 * 
                                                                                                  dmax * (e0^b)/(e0^b + x^b) - sub_params$g))), 
                         gamma_bias = sub_params$g, psychometric_type = "percent_correct")
        
        interpolated_data <- rbind(interpolated_data, d0, d1, d2, d3)
      }
    }
    
  }
  
  experiment_bin_values <- get_experiment_bin_values()
  interpolated_data <- merge(interpolated_data, experiment_bin_values, 
                             by = c("BIN", "TARGET"))
  return(interpolated_data)
}


