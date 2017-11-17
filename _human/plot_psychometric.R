#' Visualize a single psychometric function
plot_single_psychometric <- function(psychometric_parameters, 
                                     performance_measures, bin, target, out_path = "~/Dropbox/Calen/Dropbox/") {
  
  library(dplyr)
  library(ggplot2)
  sub_params <- psychometric_parameters %>%
    filter(BIN == bin, TARGET == target) %>% 
    data.frame()
  
  sub_response <- performance_measures %>%
    filter(BIN == bin, TARGET == target) %>% 
    data.frame()
  
  responses <- lapply(1:nrow(sub_params), 
                      FUN = function(x) 
                        data.frame(SUBJECT = sub_params[x, "SUBJECT"], eccentricity = seq(0, 30, 0.1), 
                                   pc = pnorm(sub_params[x, "d0"]/2 * (sub_params[x, "e0"]^sub_params[x, "b"])/(sub_params[x, "e0"]^sub_params[x, "b"] + seq(0, 30, 0.1)^sub_params[x, "b"])),
                                   dprime = (sub_params[x, "d0"] * (sub_params[x, "e0"]^sub_params[x, "b"])/(sub_params[x, "e0"]^sub_params[x, "b"] + seq(0, 30, 0.1)^sub_params[x, "b"])))) %>% 
    do.call(rbind, .)
  
  responses <- responses %>%
    mutate(BIN = bin, TARGET = target)
  
  p <- ggplot(responses, aes(x = eccentricity, y = dprime, colour = SUBJECT)) + 
    geom_line() +
    geom_point(data = sub_response %>% filter(!SUBJECT == "sps"), aes(x = eccentricity, y = dprime)) +
    ggtitle(paste0(target, " ", bin))
  
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