# Plot a figure with dprime for all measured conditions
plot_dprime <- function(dprime.df) {
  library(dplyr)
  library(ggplot2)

  p1 <- ggplot(data = dprime.df, 
               aes(x = eccentricity, y = dprime, colour = SUBJECT)) +
    geom_point() +
    facet_grid(BIN ~ TARGET, scales = "free_y")
  
  ggsave(file = '~/Dropbox/Calen/Dropbox/figure1.pdf', p1, width = 30, height = 30)
}

# Visualize the efficiency of the human observers.
plot_efficiency <- function(model.psychometrics, plotType = 1) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(grid)
  
  summarize <- dplyr::summarise
  
  source('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_model/import_model.R')
  
  
  bin.values <- get_experiment_bin_values()
  
  efficiency.df <- left_join(model.psychometrics, bin.values, by = c("TARGET", "BIN"))
  
  efficiency.df <- efficiency.df %>% mutate(efficiency = 1 / dprime_at_threshold) %>% group_by(TARGET, BIN, statType, statValue) %>% summarize(efficiency = mean(efficiency))
  
  
  
  if(plotType == 1) { # cardinal axis
    
    # Efficiency by statistic    
    fig <- efficiency.df %>%
      group_by(statType) %>%
      nest() %>%
      mutate(fig = map2(statType, data, function(statType, data) {
        efficiency.df <- data
        
        x.lab <- ifelse(statType == "Lvals", "Luminance (%)", ifelse(statType == "Cvals", "Contrast (RMS)", "Similarity"))
        
        fig.1 <- ggplot(efficiency.df, aes(x = statValue, y = efficiency, colour = TARGET)) +
          geom_point(size = 2.5) + 
          geom_line(size = 1.5) +
          theme_set(theme_bw(base_size = 35))  +# pre-set the bw theme.
          scale_color_brewer(name = "Target", palette = "Dark2") +
          theme(aspect.ratio = 1, axis.title.y = element_text(angle = 0, vjust = .5),
                panel.border = element_rect(size = 1), legend.title=element_text(size=10), 
                legend.text=element_text(size=9)) +
          expand_limits(x = c(0), y = c(0, .3)) +
          xlab(x.lab) +
          ylab(expression(sqrt(eta))) +
        ggtitle(first(data$sub_type))
        
        plot(fig.1)
        ggsave(file = paste0("~/Dropbox/Calen/Dropbox/efficiency.fig.", statType, '.', first(data$model), ".statistic.pdf"))
      }))

    
    
   } else if(plotType == 2) {
    efficiency.by.ecc <- efficiency.df %>% group_by(TARGET, eccentricity) %>% summarize(eff.avg = mean(efficiency))
    
    # Efficiency by statistic    
    fig.1 <- ggplot(efficiency.by.ecc, aes(x = eccentricity, y = eff.avg, colour = TARGET)) +
      geom_point(size = 2.5) + 
      geom_line(size = 1.5) +
      theme_set(theme_bw(base_size = 35)) +# pre-set the bw theme.
      geom_point() + 
      geom_line() +
      scale_color_brewer(name = "Target", palette = "Dark2") +
      theme(aspect.ratio = 1, axis.title.y = element_text(angle = 0, vjust = .5), legend.key.height=unit(3,"line"),
            panel.border = element_rect(size = 1)) +
      xlab("Eccentricity (º)") +
      ylab(expression(sqrt(eta)))
  
    plot(fig.1)  
    ggsave(file = "~/Dropbox/Calen/Work/occluding/detection_model_analysis/presentations/vss_2018/efficiency_eccentricity.pdf")
   }
  }

#' Compute the efficiency of the optimal and human responses
get_efficiency <- function(human.psychometric, optimal.observer) {
  library(purrrlyr)
  
  model.dat <- optimal.observer %>%
    select(SUBJECT, BIN, TARGET, eccentricity, dprime)
  
  human.dat <- human.psychometric %>%
    select(SUBJECT, BIN, TARGET, d0, e0, b, gamma) %>%
    mutate(e0 = as.numeric(e0))
  
  optim.human.dat <- merge(model.dat, human.dat, by = c("TARGET", "BIN"))

  efficiency <- by_row(optim.human.dat, function(row) {
    e0  <- row$e0
    b   <- row$b
    d0  <- row$d0
    ecc <- row$eccentricity
    
    dprime.human <- d0 * e0^b/(e0^b + ecc^b)
  }, .to = "dprime.human") %>%
    unnest(dprime.human) %>%
    mutate(efficiency = dprime.human/dprime) %>%
    select(BIN, TARGET, eccentricity, efficiency) %>%
    arrange(BIN, TARGET, eccentricity)

  bin_labels <- get_experiment_bin_values()
  efficiency <- merge(efficiency, bin_labels)
  return(efficiency)
}

get.model.cor.thresholds <- function(model.thresholds, human.thresholds) {
  bin.values <- get_experiment_bin_values()

  m <- model.thresholds %>%
    arrange(TARGET, BIN) %>%
    select(TARGET, BIN, threshold) %>%
    mutate(SUBJECT = "model")

  h <- human.thresholds %>%
    group_by(TARGET, BIN) %>%
    summarize(threshold = mean(threshold)) %>%
    mutate(SUBJECT = "human")
  
  #human.model <- full_join(m, h, by = c("TARGET", "BIN"))

  human.model <- merge(m, h, by = c("TARGET", "BIN"))  
  
  human.model <- merge(bin.values, human.model, by = c("TARGET", "BIN"))

  cor.by.group <- human.model %>%
    group_by(statType) %>%
    summarize(r2 = cor(threshold.x, threshold.y)^2)
  
  fig.dat <- rbind(m %>% mutate(type = "model"), h %>% mutate(type = "human")) %>%
    spread(key = type, value = threshold) %>%
    group_by(TARGET, statType) %>%
    mutate(human = scale(human), model = scale(model)) %>%
    merge(cor.by.group) %>%
    group_by(statType)
  
  #levels(fig.dat$statType) <- c("Luminance", "Contrast", "Similarity")
  plot.fig <- ggplot(data = fig.dat, aes(x = human, y = model, colour = TARGET)) +
    geom_point() +
    facet_grid(~statType) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_set(theme_bw(base_size = 15))  +# pre-set the bw theme.
    theme(aspect.ratio = 1) +
    expand_limits(x = c(-2, 2), y = c(-2, 2)) +
    scale_color_brewer(palette = "Dark2") +
    ylab("Model Threshold") +
    xlab("Human Threshold") +
    geom_text(data=subset(fig.dat, TARGET == 'vertical' & BIN == 3), aes(x = -1.25, y = 1.5, label = paste0('r^2 = ', as.character(round(r2, 3)))), colour = "black")
    
  
  ggsave(file = '~/Dropbox/Calen/Work/occluding/detection_model_analysis/presentations/vss_2018/threshold.scatter.pdf', plot.fig, scale = .6)
  
  plot.fig
}

get.model.cor.dprime <- function(model.dprime, human.psychometrics) {
  bin.values <- get_experiment_bin_values()

  human <- human.psychometrics %>%
    select(TARGET, BIN, e0, d0, b, gamma, SUBJECT)
  
  model <- model.dprime %>%
    mutate(dprime * 1) %>%
    ungroup() %>%
    arrange(TARGET, BIN) %>%
    select(TARGET, eccentricity, BIN, dprime)
  
  human.model <- model %>% 
    merge(., human) %>%
    mutate(model.dprime = dprime, dprime = NULL) %>% 
    mutate(human.dprime = d0 * e0^b/(e0^b + eccentricity^b)) %>%
    as_tibble() %>%
    merge(., bin_values)
  
  cor.val <- human.model %>%
    filter(!(statType == "Lvals" & BIN %in% c(1,2))) %>%
    group_by(statType) %>%
    summarize(cor(human.dprime, model.dprime)^2)
  
  fig.dat <- human.model %>%
    group_by(statType) %>%
    mutate(human.dprime = scale(human.dprime), model.dprime = scale(model.dprime))
  
  plot.fig <- ggplot(data = fig.dat, aes(x = scale(human.dprime), y = scale(model.dprime), colour = TARGET)) + 
    geom_point() + 
    facet_grid(~statType) + 
    geom_smooth(method = "lm", se = FALSE) +
    theme(aspect.ratio = 1)
  
  save(file = '~/Dropbox/Calen/Dropbox/threshold.scatter.pdf', plot.fig)
}

plot.model.psychometric <- function(model.psychometric, model.dprime) {
  library(ggplot2)
  
  model.dprime.nest <- model.dprime %>%
    group_by(TARGET, BIN, SUBJECT) %>%
    nest()
  
  model.data <- merge(model.psychometric, model.dprime.nest) %>% as_tibble()
  
  by_row(model.data, function(row) {
    empirical.dat <- row$data[[1]]
    eccentricity <- seq(0, 23, .1)
    d0 <- row$d0
    b  <- row$b
    e0 <- row$e0
    gamma <- row$gamma
    dprime.obs   <- d0 * e0^b/(e0^b + eccentricity^b)
    
    plot.dat <- data.frame(SUBJECT = row$SUBJECT, eccentricity = eccentricity, dprime = dprime.obs)
    fig      <- ggplot(empirical.dat, aes(x = eccentricity, y = dprime)) + geom_point() + geom_line(data = plot.dat, aes(x = eccentricity, y = dprime))
    
    ggsave(filename = paste0('~/Dropbox/Calen/Dropbox/', row$TARGET,'-', row$BIN, '-', row$SUBJECT, '.pdf'), device = "pdf", plot = fig)
  }, .to = "f.ggplot")
}

plot.bivariate.condition <- function(model.wide) {
  dat.fovea <- model.wide %>% filter(BIN == 3, TARGET == "vertical", eccentricity == 4.621447) %>% select(BIN, TARGET, eccentricity, data_response_vec_0,data_response_vec_1)
  
  long.frame <- by_row(dat.fovea, function(row) {
    BIN <- row$BIN
    TARGET <- row$TARGET
    eccentricity <- row$eccentricity
    
    data.abs <- row$data_response_vec_0 %>% data.frame
    data.pres <- row$data_response_vec_1 %>% data.frame
    
    data.long.abs <- data.abs %>% mutate(BIN = BIN, TARGET = TARGET, eccentricity = as.factor(paste(round(eccentricity,3), 'º')), present = "absent")
    data.long.pres <- data.pres %>% mutate(BIN = BIN, TARGET = TARGET, eccentricity = as.factor(paste(round(eccentricity,3), 'º')), present = "present")
    
    data.long <- rbind(data.long.abs, data.long.pres)
    
  }, .to = 'd.frame')
  
  all.obs <- do.call(rbind, long.frame$d.frame) %>% group_by(TARGET, BIN, eccentricity) %>%
    group_by(TARGET, BIN, eccentricity, present) %>%
    mutate(edge = scale(edge), pattern = scale(pattern), mean = scale(mean))
  
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(25)
  
  fig.edge.pattern <- ggplot(all.obs, aes(x = edge, y = pattern)) +
    stat_bin2d(bins=50) + 
    scale_fill_gradientn(colours=r) +
    theme_bw(base_size = 20) +
    #theme_bw(base_size = 15) +# pre-set the bw theme.
    theme(aspect.ratio = 1) +
    facet_wrap(~present, ncol = 2) +
    xlab("Edge Response") +
    ylab("Pattern Response") +
    coord_cartesian(xlim = c(-5,5), ylim = c(-4, 4))
  
  plot(fig.edge.pattern)
  
  ggsave(filename = '~/Dropbox/Calen/Work/occluding/detection_model_analysis/presentations/vss_2018/covariance.edge.pattern.pdf', plot = fig.edge.pattern, width = 8)
  
  #
  
  fig.edge.luminance <- ggplot(all.obs, aes(x = edge, y = mean)) +
    stat_bin2d(bins=50) + 
    scale_fill_gradientn(colours=r) +
    theme_bw(base_size = 20) +
    #theme_bw(base_size = 15) +# pre-set the bw theme.
    theme(aspect.ratio = 1) +
    facet_wrap(~present, ncol = 2) +
    xlab("Edge Response") +
    ylab("Luminance Response") +
    coord_cartesian(xlim = c(-5,5), ylim = c(-4, 4))
  
  plot(fig.edge.luminance)
  
  ggsave(filename = '~/Dropbox/Calen/Work/occluding/detection_model_analysis/presentations/vss_2018/covariance.edge.luminance.pdf', plot = fig.edge.luminance, width = 8)

  #
  
  fig.pattern.luminance <- ggplot(all.obs, aes(x = pattern, y = mean)) +
    stat_bin2d(bins=50) + 
    scale_fill_gradientn(colours=r) +
    theme_bw(base_size = 20) +
    #theme_bw(base_size = 15) +# pre-set the bw theme.
    theme(aspect.ratio = 1) +
    facet_wrap(~present, ncol = 2) +
    xlab("Pattern Response") +
    ylab("Luminance Response") +
    coord_cartesian(xlim = c(-5,5), ylim = c(-4, 4))
  
  plot(fig.pattern.luminance)
  
  ggsave(filename = '~/Dropbox/Calen/Work/occluding/detection_model_analysis/presentations/vss_2018/covariance.pattern.luminance.pdf', plot = fig.pattern.luminance, width = 8)  
  
}    

