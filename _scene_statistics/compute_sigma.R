# Import the template response data from disk.
# Description: A dataframe with sigma and mean responses for all
# natural scene conditions.
import_stats <- function(file_path = '~/Dropbox/Calen/Work/occluding/occlusion_detect/_scene_statistics/exported_scene_stats/') {
  library(dplyr, quietly = TRUE)
  library(purrrlyr, quietly = TRUE)
  library(purrr, quietly = TRUE)
  library(tidyr, quietly = TRUE)
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R')

  files <- list.files(paste0(file_path), full.names = T, recursive = F)
  files <- files[grepl('txt', files)]
  
  ecc_deg = c(0.000000, 1.381106, 4.621447, 11.509662, 23.099524)
  
  scene_stats <- lapply(files, FUN = function(x) read.table(x, header = T, sep = "\t"))
  
  scene_stats_all <- do.call(rbind, scene_stats) %>%
    rename(response = tSigma)
  
  scene_stats_all$TARGET <- factor(scene_stats_all$TARGET, levels = c(1,2,3,4), labels = c("vertical", "horizontal", "bowtie", "spot")) # create text labels for factors

  bin_values <- get_bin_values() %>%
    select(-TARGET) %>%
    dplyr::rename(TARGET = TARGET_NAME)
  
  response_stats     <- scene_stats_all %>%
    merge(., bin_values) %>%
    arrange(TARGET, Lvals, Cvals, Svals)
  
  response_stats <- by_row(response_stats, function(x) {
    ecc_deg[x$eccentricity]
  }, .to = "ecc_deg") %>%
    unnest(ecc_deg)
  
  return(response_stats)
}

# Generate figures with the standard deviations of the response statistics
plot_response_sigma <- function(response_stats.df) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrrlyr)
  
  
  fig_stats <- response_stats.df %>% 
    group_by(response_type, statistic, eccentricity) %>%
    nest()
  
  by_row(fig_stats, function(x) {
    data <- x$data[[1]]
    
    ggplot(data, aes(x = Svals, y = response, colour = as.factor(Cvals))) +
      geom_point() +
      facet_wrap(Lvals~TARGET, scales = "free_y", ncol = 4) +
      theme(aspect.ratio = 1)
    
    ggsave(filename = paste0('~/Dropbox/Calen/Dropbox/scene_statistics_figures/', x$eccentricity[[1]], '-', x$response_type[[1]], '-', x$statistic[[1]], '.pdf'), width = 30, height = 30)
  })
}

# Fit to Sebastian, 2016
plot_sebastian_fit <- function(response_stats.df) {
  # Load required libraries
  library(bbmle, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(knitr, quietly = TRUE)
  library(kableExtra, quietly = TRUE)
  
  # Filter statistics to contain approx. range used in the Sebastian 2016
  scene.stats.pnas <- response_stats.df %>%
    filter(Lvals < 60, Cvals < .30, Svals < .52) %>%
    filter(!is.nan(response), statistic == "sigma", response_type == "tMatch")
  
  # Fit the separable model to the range ranges in Sebastian 2016.
  stats.nested.pnas <- scene.stats.pnas %>%
    group_by(TARGET) %>%
    nest() %>%
    mutate(models = map(data, function(x) {
      m1 <- mle2(response ~ dnorm(mean = k0 * (I(Lvals) + a)*(I(Cvals) + b)*(I(Svals) + c), sd = 1),
                 start = list(k0 = 0, a = 0,b = 0,c = 0),
                 data = x)
    })) %>%
    mutate(prediction = map(models, predict)) 

  # Compute predictions for the data from the fitted model
  stats.nested.pnas.prediction <- stats.nested.pnas %>%
    unnest(data, prediction) %>%
    gather(type, value, c("response", "prediction"))
  
  # Extract parameters from fitted model.
  stats.nested.pnas.parameters <- stats.nested.pnas %>%
    mutate(parameters = map(models, . %>% coef %>% as.list() %>% as_tibble())) %>%
    unnest(parameters)
  
  # Print the model parameters to disk.
  stats.nested.pnas.parameters %>% 
    select(TARGET, k0, a, b, c) %>%
    kable(format = 'latex') %>%
    kable_as_image(filename = '~/Dropbox/Calen/Dropbox/fovea_parameters_pnas', file_format = 'png')  
  
  # Create figure with scene statistics and predictions.
  pnas.range <- ggplot(stats.nested.pnas.prediction, aes(x = Svals, y = value, colour = as.factor(Cvals), shape = TARGET)) +
    geom_point(data = . %>% filter(type == "response"), size = 5) +
    geom_line(data = . %>% filter(type == "prediction")) +
    facet_wrap(~Lvals, ncol= 4) +
    theme_set(theme_gray(base_size = 30)) +
    list(theme(legend.position = "bottom",
               legend.title = element_text(size=20),
               plot.title = element_text(size = 30),
               axis.title = element_text(size=20),
               axis.text = element_text(size=15),
               legend.text = element_text(size=15))) +
    theme(aspect.ratio = 1)
  
  # Store a figure to disk. 
  ggsave(filename = '~/Dropbox/Calen/Dropbox/pnas.range.pdf', width = 20, height = 20)
}

# Fit a function to the standard standard deviations of the template responses.
fit.template.stats <- function(scene_statistics, eccLvl = c(1)) {
  library(bbmle)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(knitr)
  library(kableExtra)
  library(RColorBrewer)
  library(ggthemes)
  library(DEoptim)
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/plot_theme.R')
  
  pal <- brewer.pal(10, "Set3")
  
  # Filter data frame to contain relevant data
  template.stats <- scene_statistics %>%
    filter(eccentricity %in% eccLvl, response_type == "tMatch", statistic == "sigma", !is.nan(response)) %>%
    group_by(TARGET, eccentricity, ecc_deg) %>%
    nest()
    
  # Fit model
  m.1 <- template.stats %>% 
    mutate(prediction = map(data, function(data) {

      
      f <- function(data) {
        Lvals <- data$Lvals
        Cvals <- data$Cvals
        Svals <- data$Svals
        
        f.return <- function(x) {
          k0 <- x[1]
          L_p <- x[2]
          C_p <- x[3]
          S_p <- x[4]
          S_exp <- x[5]
          
          val <- k0 * ((Lvals) + L_p)*((Cvals) + C_p)*((Svals) + S_p)^S_exp
        }
      }
      
      g   <- f(data) 
      
      f.1 <- function(x) {
        response <- data$response
        response_hat <- g(x)
        
        sum(sqrt((response_hat - response)^2), na.rm = T)
      }
      
      model <- DEoptim(f.1, lower = c(k0 = 0,L_p = 0,C_p = 0,S_p = 0,S_exp = 0), upper = c(k0 = 30,L_p = 5,C_p = 5,S_p = 5,S_exp = 30))
      
      parameters <- model$optim$bestmem
      
      prediction <- list(parameters = parameters, prediction = g(parameters))
      }))
  
  m.1$parameters <- map(m.1$prediction,1)
  m.1$prediction <- map(m.1$prediction,2)
  
  m.2 <- m.1 %>% unnest(data, prediction, .drop = F)
  
  model.response <- m.2 %>% 
    gather("response", "value", prediction, response)

  # Plot Template Response Statistics
  model.response %>%
    filter(L == 5) %>%
    rename(Contrast = Cvals, Similarity = Svals) %>%
    group_by(ecc_deg) %>%
    nest() %>%
    mutate(fig.vis = map2(ecc_deg, data, function(ecc_arg, model_arg) {
      fig.1 <- ggplot(data = model_arg, aes(x = Similarity , y = value, colour = as.factor(Contrast))) +
        geom_point(data = model_arg %>% filter(response == "response"), size = 4) +
        geom_line(data = model_arg %>% filter(response == "prediction"), size = 1.5) +
        facet_wrap(~ TARGET, ncol = 2) +
        theme_set(theme_bw(base_size = 20))  +# pre-set the bw theme.
        theme(aspect.ratio = 1) +
        theme.grant +
        expand_limits(x = c(.45, .85), y = c(0, .015)) +
        ylab(expression("Template Response"~ (sigma))) +
        scale_colour_discrete(name = "Contrast (RMS)")
      
      plot(fig.1)
  
      ggsave(plot = fig.1, filename = paste0('~/Dropbox/Calen/Dropbox/template_sigma_', ecc_arg,'.pdf'), scale = 1.35, useDingbats = FALSE)
    }))
  
  fitted.params <- cbind(m.1[, c(1,2,3)], do.call(rbind, m.1$parameters))
  
  fitted.params %>%
    kable(format = 'latex') %>%
    kable_as_image('~/Dropbox/Calen/Dropbox/template_sigma_params', file_format = 'pdf')
  
  fitted.params.long <- fitted.params %>%
    gather("param_label", "v", 4:8)
  
  return(model.response)
}

# Fit a function to the standard standard deviations of the template responses.
fit.edge.stats <- function(scene_statistics, eccLvl = c(1)) {
  library(bbmle)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(knitr)
  library(kableExtra)
  library(RColorBrewer)
  library(ggthemes)
  library(DEoptim)
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/plot_theme.R')
  
  pal <- brewer.pal(10, "Set3")
  
  # Filter data frame to contain relevant data
  edge.stats <- scene_statistics %>%
    filter(eccentricity %in% eccLvl, response_type == "Eabs", statistic == "sigma", !is.nan(response)) %>%
    group_by(TARGET, L, C, Lvals, Cvals, ecc_deg) %>% 
    summarize(response = mean(response)) %>%
    group_by(TARGET, ecc_deg) %>%
    nest()
  
  # Fit model
  m.1 <- edge.stats %>% 
    mutate(prediction = map(data, function(data) {
      f <- function(data) {
        Lvals <- data$Lvals
        Cvals <- data$Cvals
        
        f.return <- function(x) {
          k0 <- x[1]
          a <- x[2]
          b <- x[3]
          c <- x[4]
          d <- x[5]
          
          val <- k0 * (I(Lvals) + b)^a*(I(Cvals) + d)^c
        }
      }
      
      g   <- f(data) 
      
      f.1 <- function(x) {
        response <- data$response
        response_hat <- g(x)
        
        sum(sqrt((response_hat - response)^2), na.rm = T)
      }
      
      model <- DEoptim(f.1, lower = c(k0 = 0,a = 0,b = 0,c = 0,d = 0), upper = c(k0 = 500,a = 5,b = 5,c = 5,d = 30))
      
      parameters <- model$optim$bestmem
      
      prediction <- list(parameters = parameters, prediction = g(parameters))
    }))
  
  m.1$parameters <- map(m.1$prediction,1)
  m.1$prediction <- map(m.1$prediction,2)
  
  m.2 <- m.1 %>% unnest(data, prediction, .drop = F)
  
  model.response <- m.2 %>% 
    gather("response", "value", prediction, response)
  
  fig <- model.response %>%
    rename(Luminance = Lvals, Contrast = Cvals) %>%
    group_by(ecc_deg) %>%
    nest() %>%
    mutate(fig.vis = map2(ecc_deg, data, function(ecc_arg, model_arg) {
      fig.1 <- ggplot(data = model_arg, aes(x = Contrast , y = value, colour = as.factor(Luminance))) +
        geom_point(data = model_arg %>% filter(response == "response"), size = 4) +
        geom_line(data = model_arg %>% filter(response == "prediction"), size = 1.5) +
        facet_wrap(~ TARGET, ncol = 2) +
        theme_set(theme_bw(base_size = 20))  +# pre-set the bw theme.
        theme(aspect.ratio = 1) +
        theme.grant +
        expand_limits(x = c(.45, .85), y = c(0, .015)) +
        ylab(expression("Edge Response"~ (sigma))) +
        scale_colour_discrete(name = "Luminance (%)")
      
      plot(fig.1)
      
      ggsave(plot = fig.1, filename = paste0('~/Dropbox/Calen/Dropbox/edge_sigma_', ecc_arg,'.pdf'), scale = 1.35, useDingbats = FALSE)
    }))
  
  
  fitted.params <- cbind(m.1[, c(1,2)], do.call(rbind, m.1$parameters))
  
  fitted.params %>%
    kable(format = 'latex') %>%
    kable_as_image(., filename = '~/Dropbox/Calen/Dropbox/edge_sigma_params', file_format = 'pdf')
  
  fitted.params.long <- fitted.params %>%
    gather("param_label", "v", 3:7)
  
  return(model.response)
}

# Luminance Statistics
# Fit a function to the standard standard deviations of the template responses.
fit.L.stats <- function(scene_statistics, eccLvl = c(1)) {
  library(bbmle)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(knitr)
  library(kableExtra)
  library(RColorBrewer)
  
  pal <- brewer.pal(10, "Set3")
  
  # Filter data frame to contain relevant data
  L.stats <- scene_statistics %>%
    filter(eccentricity %in% eccLvl, response_type == "L", statistic == "sigma", !is.nan(response)) %>%
    group_by(TARGET, eccentricity, ecc_deg) %>%
    nest()
  
  # Fit model
  m.1 <- L.stats %>% 
    mutate(prediction = map(data, function(data) {
      
      
      f <- function(data) {
        Lvals <- data$Lvals
        Cvals <- data$Cvals
        Svals <- data$Svals
        
        f.return <- function(x) {
          k0 <- x[1]
          L_p <- x[2]
          C_p <- x[3]
          S_p <- x[4]
          S_exp <- x[5]
          
          val <- k0 * ((Lvals) + L_p)*((Cvals) + C_p)^S_exp*((Svals) + S_p)
        }
      }
      
      g   <- f(data) 
      
      f.1 <- function(x) {
        response <- data$response
        response_hat <- g(x)
        
        sum(sqrt((response_hat - response)^2), na.rm = T)
      }
      
      model <- DEoptim(f.1, lower = c(k0 = 0,L_p = 0,C_p = 0,S_p = 0,S_exp = 0), upper = c(k0 = 5,L_p = 5,C_p = 5,S_p = 5,S_exp = 30))
      
      parameters <- model$optim$bestmem
      
      prediction <- list(parameters = parameters, prediction = g(parameters))
    }))
  
  m.1$parameters <- map(m.1$prediction,1)
  m.1$prediction <- map(m.1$prediction,2)
  
  m.2 <- m.1 %>% unnest(data, prediction, .drop = F)
  
  model.response <- m.2 %>% 
    gather("response", "value", prediction, response)
  
  # Plot Template Response Statistics
  model.response %>%
    filter(C == 5) %>%
    rename(Similarity = Svals, Luminance = Lvals) %>%
    group_by(ecc_deg) %>%
    nest() %>%
    mutate(fig.vis = map2(ecc_deg, data, function(ecc_arg, model_arg) {
      fig.1 <- ggplot(data = model_arg, aes(x = Similarity , y = value, colour = as.factor(Luminance))) +
        geom_point(data = model_arg %>% filter(response == "response")) +
        geom_line(data = model_arg %>% filter(response == "prediction")) +
        facet_wrap(~ TARGET, ncol = 2) +
        theme_set(theme_bw(base_size = 15))  +# pre-set the bw theme.
        theme(aspect.ratio = 1) +
        expand_limits(x = c(.45, .85), y = c(0, .015)) +
        scale_color_manual(name = "Contrast (RMS)",values = pal) +
        ylab("Template Response (sd)")
      
      plot(fig.1)
      
      ggsave(plot = fig.1, filename = paste0('~/Dropbox/Calen/Dropbox/template_sigma_', ecc_arg,'.pdf'), scale = 1.35, useDingbats = FALSE)
    }))
  
  fitted.params <- cbind(m.1[, c(1,2,3)], do.call(rbind, m.1$parameters))
  
  fitted.params %>%
    kable(format = 'latex') %>%
    kable_as_image('~/Dropbox/Calen/Dropbox/template_sigma_params', file_format = 'pdf')
  
  fitted.params.long <- fitted.params %>%
    gather("param_label", "v", 4:8)
  
  return(model.response)
}
