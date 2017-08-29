#' Import human data
#'
#' @param path
#'
#' @return
#' @export
#'
#' @import dplyr magrittr
#'
#' @examples
get_human_responses <- function(path = '~/Dropbox/Calen/Work/natural_masking_aux/_data/exported/human_data.txt') {
    library(dplyr)
    library(tidyr)
  
    bin_values <- get_experiment_bin_values()
    
    human_data <- read.table(path, sep = "\t", header = T)
    
    human_data <- human_data %>% filter(!SUBJECT %in% c("jsa", 
        "yhb"), TRIAL != 1) %>% 
      rename(BIN = BINS)
    
    human_data <- human_data %>% group_by(SUBJECT, BIN, TARGET) %>% 
        mutate(n_ecc = length(unique((ECCENTRICITY)))) %>% filter(!(SESSION == 
        2 & n_ecc > 5))
    
    human_data$CORRECT <- ifelse(human_data$HIT == 1 | human_data$CORRECTREJECTION == 
        1, 1, 0)
    
    human_data$BIN <- factor(human_data$BIN)
    
    human_data <- merge(human_data, bin_values) %>% arrange(SUBJECT, 
        TARGET, BIN, statType, SESSION, ECCENTRICITY)
    
    return(human_data)
}

#' Convert human_data into a dataframe containing percent correct data
get_human_detect <- function(human_data) {
    if (missing(human_data)) {
        error("Missing human data")
    }
    experiment_bin_values <- get_experiment_bin_values()
    
    human_detect <- human_data %>% group_by(SUBJECT, ECCENTRICITY, 
        BIN, TARGET) %>% summarize(hit = sum(HIT == 1)/(sum(HIT == 
        1) + sum(MISS == 1)), miss = sum(MISS == 1)/(sum(MISS == 
        1) + sum(HIT == 1)), falsealarm = sum(FALSEALARM == 1)/(sum(FALSEALARM == 
        1) + sum(CORRECTREJECTION == 1)), correctrejection = sum(CORRECTREJECTION == 
        1)/(sum(CORRECTREJECTION == 1) + sum(FALSEALARM == 1)), 
        percent_correct = mean(CORRECT)) %>% mutate(hit_adj = ifelse(hit == 
        0, 1/1200, ifelse(hit == 1, 1199/1200, hit)), falsealarm_adj = ifelse(falsealarm == 
        0, 1/1200, ifelse(falsealarm == 1, 1199/1200, falsealarm)), 
        percent_correct_adj = ifelse(percent_correct == 1, 2399/2400, 
            percent_correct), percent_correct_adj = ifelse(percent_correct == 
            0, 1/2400, percent_correct)) %>% mutate(dprime = qnorm(hit_adj) - 
        qnorm(falsealarm_adj)) %>% mutate(bias = dprime/2 - qnorm(hit_adj)) %>% 
        merge(., experiment_bin_values, by = c("BIN", "TARGET")) %>% 
        rename(eccentricity = ECCENTRICITY) %>% arrange(SUBJECT, 
        BIN, TARGET, eccentricity, percent_correct) %>% data.frame()
    
    return(human_detect)
}

#' Estimate psychometric parameters for human detection experiment.
get_human_psychometric_params <- function(human_responses) {
    library(bbmle)
    library(parallel)
  
    # Remove replicated rows
    human_responses <- human_responses[, setdiff(names(human_responses), 
        c("L", "C", "S", "statType", "statValue"))] %>%
      distinct()
    
    # Likelihood function.
    f <- function(d0, e0, b, g) {
        likelihood <- sum(log(pnorm(0.5 * d0 * (e0^b)/(e0^b + 
            hit_vec^b) - g))) + sum(log(pnorm(-0.5 * d0 * (e0^b)/(e0^b + 
            fa_vec^b) - g))) + sum(log(pnorm(-0.5 * d0 * (e0^b)/(e0^b + 
            miss_vec^b) + g))) + sum(log(pnorm(0.5 * d0 * (e0^b)/(e0^b + 
            cr_vec^b) + g)))
        
        
        nll <- -likelihood
        
        if (is.infinite(nll)) {
            return(10000)
        } else {
            return(nll)
        }
        
    }
    
    get_params <- function(human_response_condition) {
        hit_vec <- human_response_condition[human_response_condition$HIT == 
            1, "ECCENTRICITY"]
        fa_vec <- human_response_condition[human_response_condition$FALSEALARM == 
            1, "ECCENTRICITY"]
        cr_vec <- human_response_condition[human_response_condition$CORRECTREJECTION == 
            1, "ECCENTRICITY"]
        miss_vec <- human_response_condition[human_response_condition$MISS == 
            1, "ECCENTRICITY"]
        
        environment(f) <- environment()
        return(f)
        
    }
    
    human_response_list <- human_responses %>% 
      group_by(TARGET, 
        BIN, SUBJECT) %>% do(vals = data.frame(.)) %>% 
      select(vals) %>% 
        as.list()
    human_response_list <- human_response_list[[1]]
    
    fcn_vec <- lapply(human_response_list, FUN = function(x) get_params(x))
    
    start_params <- expand.grid(d0 = 4.5, e0 = seq(1, 20, 0.5), 
        b = seq(0, 5, 0.5), g = seq(-2, 2, 0.2))
    n_search <- nrow(start_params)
    
    # Grid search for best starting parameters
    grid_params <- lapply(fcn_vec, FUN = function(x) cbind(start_params, 
        nll = unlist(mclapply(1:n_search, FUN = function(y) x(d0 = start_params[y, 
            1], e0 = start_params[y, 2], b = start_params[y, 
            3], g = start_params[y, 4]), mc.cores = 16))) %>% 
        filter(nll == min(nll)))
    
    start_params <- do.call(rbind, grid_params)
    subject_df <- lapply(human_response_list, FUN = function(x) x[1, 
        c("TARGET", "BIN", "SUBJECT")]) %>% do.call(rbind, .)
    
    # ---
    # Step 1.    
    # Maximum likelihood for parameters, all parameters free to
    # vary.
    p.1 <- mclapply(1:length(fcn_vec), FUN = function(x) {
        y <- start_params[x, ]
        mle2(fcn_vec[[x]], start = list(e0 = y$e0, b = y$b, g = y$g), 
            fixed = list(d0 = 4.5)) %>% coef
    }, mc.cores = 16) %>% do.call(rbind, .) %>% data.frame %>% 
        cbind(subject_df, .)  # 
    
    # --- 
    # Step 2.
    p.1 <- p.1 %>% group_by(SUBJECT) %>% mutate(b = mean(b))
    
    # Maximum likelihood for parameters, only e0 can vary.  b
    # fixed at the mean.
    p.2 <- mclapply(1:length(fcn_vec), FUN = function(x) {
        y <- p.1[x, ]
        mle2(fcn_vec[[x]], start = list(e0 = y$e0, g = y$g), 
            fixed = list(d0 = 4.5, b = y$b)) %>% coef
    }, mc.cores = 16) %>% do.call(rbind, .) %>%
      data.frame %>% 
        cbind(subject_df, .)
    
    # ---
    # Step 3.
    # All parameters fixed, but beta varies.
    p.3 <- mclapply(1:length(fcn_vec), FUN = function(x) {
        y <- p.2[x, ]
        mle2(fcn_vec[[x]], start = list(b = y$b), fixed = list(d0 = 4.5, 
            e0 = y$e0, g = y$g)) %>% coef
    }, mc.cores = 16) %>% do.call(rbind, .) %>% data.frame %>% 
        cbind(subject_df, .)
    
    # ---
    # Step 4.
    p.3 <- p.3 %>%
      group_by(SUBJECT) %>%
      mutate(b = mean(b))
    
    # All parameters vary, beta fixed at the mean.
    p.4 <- mclapply(1:length(fcn_vec), FUN = function(x) {
        y <- p.3[x, ]
        mle2(fcn_vec[[x]], start = list(e0 = y$e0, g = y$g), 
            fixed = list(d0 = 4.5, b = y$b)) %>% coef
    }, mc.cores = 16) %>% do.call(rbind, .) %>%
      data.frame %>% 
        cbind(subject_df, .)
    
    # --- 
    # Step 5.
    # Return final model.
    fitted.params <- p.4
    return(fitted.params)
}

#' Bootstrap thresholds for human psychometric functions.
bootstrap_thresholds_human <- function(human_responses) {
  # Remove replicated rows
  human_responses <- human_responses[, setdiff(names(human_responses), 
                                               c("L", "C", "S", "statType", "statValue"))] %>%
    distinct()
  
  # Likelihood function.
  f <- function(d0, e0, b, g) {
    likelihood <- sum(log(pnorm(0.5 * d0 * (e0^b)/(e0^b + hit_vec^b) - g))) +
      sum(log(pnorm(-0.5 * d0 * (e0^b)/(e0^b + fa_vec^b) - g))) +
      sum(log(pnorm(-0.5 * d0 * (e0^b)/(e0^b + miss_vec^b) + g))) +
      sum(log(pnorm(0.5 * d0 * (e0^b)/(e0^b + cr_vec^b) + g)))
    
    
    nll <- -likelihood
    
    if (is.infinite(nll)) {
      return(10000)
    } else {
      return(nll)
    }
    
  }  
  
  get_params <- function(human_response_condition) {
    hit_vec <- human_response_condition[human_response_condition$HIT == 
                                          1, "ECCENTRICITY"]
    fa_vec <- human_response_condition[human_response_condition$FALSEALARM == 
                                         1, "ECCENTRICITY"]
    cr_vec <- human_response_condition[human_response_condition$CORRECTREJECTION == 
                                         1, "ECCENTRICITY"]
    miss_vec <- human_response_condition[human_response_condition$MISS == 
                                           1, "ECCENTRICITY"]
    
    environment(f) <- environment()
    return(f)
    
  }
  
  human_response_list <- human_responses %>% 
    group_by(TARGET, 
             BIN, SUBJECT) %>% do(vals = data.frame(.)) %>% 
    select(vals) %>% 
    as.list()
  human_response_list <- human_response_list[[1]]
  
  fcn_vec <- lapply(human_response_list, FUN = function(x) get_params(x))
  
  start_params <- expand.grid(d0 = 4.5, e0 = seq(1, 20, 0.5), 
                              b = seq(0, 5, 0.5), g = seq(-2, 2, 0.2))
  n_search <- nrow(start_params)  
  
  
  # Grid search for best starting parameters
  grid_params <- lapply(fcn_vec, FUN = function(x) cbind(start_params,
                                                         nll = unlist(mclapply(1:n_search, FUN =
                                                                                 function(y) x(d0 = start_params[y, 1], e0 = start_params[y, 2], b = start_params[y, 3], g = start_params[y, 4]), mc.cores = 16))) %>% 
                          filter(nll == min(nll)))
  
  start_params <- do.call(rbind, grid_params)
  subject_df <- lapply(human_response_list, FUN = function(x) x[1, 
                                                                c("TARGET", "BIN", "SUBJECT")]) %>% do.call(rbind, .)  
  
  
  boot_data <- function(human_responses) {
    boot_response <- human_responses %>% 
      group_by(TARGET,BIN,SUBJECT) %>%
      sample_frac(1, replace = TRUE)
    
    human_response_list <- boot_response %>% 
      group_by(TARGET, 
               BIN, SUBJECT) %>% do(vals = data.frame(.)) %>% 
      select(vals) %>% 
      as.list()
    
    human_response_list <- human_response_list[[1]]
    
    fcn_vec <- lapply(human_response_list, FUN = function(x) get_params(x))    
  }
  
  boot_fcn <- function() {
    fcn_vec <- boot_data(human_responses)
    
    # ---
    # Step 1.    
    # Maximum likelihood for parameters, all parameters free to
    # vary.
    p.1 <- mclapply(1:length(fcn_vec), FUN = function(x) {
      y <- start_params[x, ]
      mle2(fcn_vec[[x]], start = list(e0 = y$e0, b = y$b, g = y$g), 
           fixed = list(d0 = 4.5)) %>% coef
    }, mc.cores = 16) %>% do.call(rbind, .) %>% data.frame %>% 
      cbind(subject_df, .)  #     

    
    thresholds <- get_threshold(p.1)
        
    return(thresholds)
  }
  
  boot_samples <- replicate(10, boot_fcn(), simplify = F)
  
  all_boot <- do.call(rbind, boot_samples)
}

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
    measure = "pc", out_path = "~/Dropbox/Calen/Dropbox/tmp_images/") {
    
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

plot_thresholds <- function(threshold_values, out_path = "~/Dropbox/Calen/Dropbox/tmp_images/") {
    bin_values <- get_experiment_bin_values()
    
    d.1 <- merge(threshold_values, bin_values)
    
    t.1.plot <- ggplot(data = d.1, aes(x = statValue, y = threshold, 
        colour = SUBJECT)) + geom_point() + facet_wrap(statType ~ 
        TARGET, scales = "free_x") + geom_line()
    
    ggsave(last_plot(), file = paste0(out_path, "thresholds.pdf"))
}

#' Add thresholds to a dataframe that contains psychometric parameters
#'
#' @param psychometric_parameters
#' @param threshold_value
get_threshold <- function(psychometric_parameters) {
    p.1 <- psychometric_parameters %>% 
      rowwise() %>% 
      mutate(threshold = ((d0 * 
        e0^b)/1 - e0^b)^(1/b))
}

plot_publication_thresholds <- function(human.thresholds, model.thresholds, 
    out_path = "~/Dropbox/Calen/Dropbox/tmp_images/", statIn = "Lvals") {
    library(ggthemes)
    
    human.threshold.1 <- human.thresholds %>% group_by(TARGET, 
        BIN) %>% summarize(threshold = mean(threshold)) %>% mutate(SUBJECT = "human") %>% 
        data.frame
    
    model.threshold.1 <- model.thresholds %>% select(TARGET, 
        BIN, threshold, SUBJECT, -function_name) %>% data.frame
    
    threshold_values <- rbind(model.threshold.1, human.threshold.1) %>% 
        data.frame
    
    bin_values <- get_experiment_bin_values()
    d.1 <- merge(threshold_values, bin_values) %>% group_by(TARGET, 
        BIN, statValue, statType, SUBJECT) %>% summarize(threshold = mean(threshold)) %>% 
        filter(statType == statIn)
    
    
    t.1.plot <- ggplot(data = d.1, aes(x = statValue, y = threshold, 
        colour = SUBJECT)) + geom_point() + geom_line() + facet_grid(~TARGET) + 
        theme_light() + theme(aspect.ratio = 1) + ylab("Eccentricity Threshold")
    
    if (statIn == "Lvals") {
        t.1.plot <- t.1.plot + ylim(0, 25) + expand_limits(x = 0, 
            y = 1) + xlab("Background Luminance")
    } else if (statIn == "Cvals") {
        t.1.plot <- t.1.plot + ylim(0, 25) + expand_limits(x = 0, 
            y = 1) + xlab("Background Contrast")
    } else if (statIn == "Svals") {
        t.1.plot <- t.1.plot + ylim(0, 25) + xlim(0.4, 0.9) + 
            expand_limits(x = 0, y = 0.9) + xlab("Background Similarity")
    }
    
    ggsave(last_plot(), file = paste0(out_path, statIn, "_thresholds.pdf"))
}

