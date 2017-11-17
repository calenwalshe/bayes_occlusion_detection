# Import the template response data from disk.
import_template_stats <- function(file_path = '~/Dropbox/Calen/Dropbox/') {
  
  template.sigma  <- read.table(paste0(file_path, 'tMatch_sigma.txt'), header = T, sep = '\t') %>% 
    mutate(type = "sigma")
  template.mean   <- read.table(paste0(file_path, 'tMatch_mean.txt'), header = T, sep = '\t') %>% 
    mutate(type = "mean")
  
  template.stats     <- rbind(template.sigma, template.mean) %>% 
    spread(type, tSigma) %>%
    mutate(PRESENT = 0)
  
  template.stats$TARGET <- factor(template.stats$TARGET, labels = c("vertical", "horizontal", "bowtie", "spot"))
  
  bin.values        <- get_bin_values() %>%
    mutate(TARGET = TARGET_NAME, TARGET_NAME = NULL)

  template.stats     <- template.stats %>%
    merge(., bin.values)
  
  return(template.stats)
}

# Import the template response data from disk.
import_luminance_stats <- function(file_path = '~/Dropbox/Calen/Dropbox/') {
  
  luminance.sigma  <- read.table(paste0(file_path, 'L_sigma.txt'), header = T, sep = '\t') %>% 
    mutate(type = "sigma")
  luminance.mean   <- read.table(paste0(file_path, 'L_mean.txt'), header = T, sep = '\t') %>% 
    mutate(type = "mean")
  
  luminance.stats     <- rbind(luminance.sigma, luminance.mean) %>% 
    spread(type, tSigma) %>%
    mutate(PRESENT = 0)
  
  luminance.stats$TARGET <- factor(luminance.stats$TARGET, labels = c("vertical", "horizontal", "bowtie", "spot"))
  
  bin.values        <- get_bin_values() %>%
    mutate(TARGET = TARGET_NAME, TARGET_NAME = NULL)
  
  luminance.stats     <- luminance.stats %>%
    merge(., bin.values)
  
  return(luminance.stats)
}

# Import the edge response data from disk.
import_edge_stats <- function(file_path = '~/Dropbox/Calen/Dropbox/') {
  
  edge.sigma  <- read.table(paste0(file_path, 'Eabs_sigma.txt'), header = T, sep = '\t') %>% 
    mutate(type = "sigma")
  edge.mean   <- read.table(paste0(file_path, 'Eabs_mean.txt'), header = T, sep = '\t') %>% 
    mutate(type = "mean")
  
  edge.stats     <- rbind(edge.sigma, edge.mean) %>% 
    spread(type, tSigma) %>%
    mutate(PRESENT = 0)
  
  edge.stats$TARGET <- factor(edge.stats$TARGET, labels = c("vertical", "horizontal", "bowtie", "spot"))
  
  bin.values        <- get_bin_values() %>%
    mutate(TARGET = TARGET_NAME, TARGET_NAME = NULL)
  
  edge.stats     <- edge.stats %>%
    merge(., bin.values)
  
  return(edge.stats)
}

# Compute and save the edge standard deviations.
edge_stat <- function(statType = 'Eabs', target = 'vertical') {
  library(dplyr)
  library(ggplot2)
  
  edge.stats <- import_edge_stats() %>%
    filter(TARGET == target)
  
  edge.all       <- edge.stats %>%
    group_by(L,C,S,PYRAMIDLVL,TARGET) %>%
    arrange(L,C,S,TARGET, PRESENT)
  
  #template.stats$Cvals <- as.factor(template.stats$Cvals)
  edge.all$Svals <- as.factor(edge.all$Svals)
  
  fig <- ggplot(data = edge.all, aes(x = Svals, y = sigma, shape = TARGET, colour = as.factor(Cvals))) +
    geom_point() + 
    facet_wrap(~Lvals, ncol = 3) +
    theme(aspect.ratio = 1)# +
    #scale_color_brewer(palette="Spectral")
  
  ggsave(fig, file = paste0('~/Dropbox/Calen/Dropbox/', target, '_', statType, '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig)
}

# Compute and save the template response standard deviations.
template_stat <- function(statType = 'tSigma', target = "vertical") {
  
  library(dplyr)
  library(ggplot2)
  
  template.stats <- import_template_stats() %>%
    filter(TARGET == target)
  
  template.all       <- template.stats %>%
    group_by(L,C,S,PYRAMIDLVL,TARGET) %>%
    arrange(L,C,S,TARGET, PRESENT)
  
  #template.stats$Cvals <- as.factor(template.stats$Cvals)
  template.all$Svals <- as.factor(template.all$Svals)
  
  fig <- ggplot(data = template.all, aes(x = Cvals, y = sigma, shape = TARGET, colour = Svals)) +
    geom_point() + 
    facet_wrap(~Lvals, ncol = 3, scales = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig, file = paste0('~/Dropbox/Calen/Dropbox/', target, '_', statType, '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig)
}

# Luminance statistics
lum_stat <- function(statType = 'L', target = "vertical") {
  library(dplyr)
  library(ggplot2)
  
  template.stats <- import_luminance_stats() %>%
    filter(TARGET == target)
  
  template.all       <- template.stats %>%
    group_by(L,C,S,PYRAMIDLVL,TARGET) %>%
    arrange(L,C,S,TARGET, PRESENT)
  
  #template.stats$Cvals <- as.factor(template.stats$Cvals)
  template.all$Svals <- as.factor(template.all$Svals)
  
  fig <- ggplot(data = template.all, aes(x = Cvals, y = sigma, shape = TARGET, colour = Svals)) +
    geom_point() + 
    facet_wrap(~Lvals, ncol = 3) +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig, file = paste0('~/Dropbox/Calen/Dropbox/', target, '_', statType, '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig)
}

# Function edge_historgram
# Description: Compute histograms of edge responses in bins. 
edge_histogram <- function(target = 1) {
  library(ggplot2)
  library(dplyr)
  
  lum <- 5
  
  alledgeSigmaAbs <- read.table('~/Dropbox/Calen/Dropbox/Epres_AllResponse_mean.txt', header = T, sep = '\t') %>%
    mutate(PRESENT = as.factor(0))
  
  edgeStatistics     <- alledgeSigmaAbs %>% 
    filter(TARGET == target)
  
  fig_lum <- ggplot(data = edgeStatistics, aes(x = resp, colour = as.factor(lum))) +
    geom_freqpoly() +
    facet_wrap(~C, scale = "free_y") +
    theme(aspect.ratio = 1)
  
  fig_lum <- ggplot(data = edgeStatistics %>% filter(L == lum)
                    , aes(x = resp, colour = as.factor(lum))) +
    geom_histogram(binwidth = .01) +
    facet_wrap(~C, scale = "free_y") + 
    theme(aspect.ratio = 1)  
  
  ggsave(fig_lum, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'edge_hist_', lum, target, '_', '.pdf'), width = 50, height = 50, units = "cm")
}

# Fit a function to the standard standard deviations of the template responses.
fit.template.stats <- function(target) {
  template.stats <- import_template_stats() %>%
    filter(!is.nan(sigma)) %>%
    filter(TARGET == target)

  Lvals <- template.stats$Lvals 
  Cvals <- template.stats$Cvals
  Svals <- template.stats$Svals
  
  sigma <- template.stats$sigma

  m1 <- mle2(sigma ~ dnorm(mean = k0 * (Lvals + a)*(I(Cvals) + b)*(I(Svals)^d + c), sd = 1),
             start = list(k0 = 1, a = 0,b = 0,c = 0, d = 0),
             data = data.frame(Lvals = Lvals, Cvals = Cvals, Svals = Svals, sigma = sigma))
  
  print((var(sigma) - var(residuals(m1)))/var(sigma))
  
  template.stats.predict       <- template.stats
  template.stats.predict$sigma <- as.numeric(predict(m1))
  
  template.stats.predict       <- template.stats.predict %>%
    group_by(L,C,S,PYRAMIDLVL,TARGET) %>%
    arrange(L,C,S,TARGET, PRESENT)
  
  template.stats.predict$Svals <- as.numeric(template.stats.predict$Svals)
  template.stats$Svals <- as.numeric(template.stats$Svals)
  
  # Plot Edge Statistic
  
  fig <- ggplot(data = template.stats, aes(x = Svals , y = (sigma), colour = as.factor(Cvals))) +
    geom_point() + 
    geom_line(data = template.stats.predict, aes(x = Svals, y = (sigma),  colour = as.factor(Cvals))) + 
    facet_wrap(~Lvals, ncol = 3) +
    theme(aspect.ratio = 1)# +
  
  plot(fig)
  
  ggsave(file = paste0('~/Dropbox/Calen/Dropbox/template_', target, '_separable.pdf'), fig, width = 40, height = 40)
}

# Fit a function to the standard deviation of the edge responses.
fit.edge.stats <- function(target) {
  library(bbmle)
  
  edge.stats <- import_edge_stats() %>%
    filter(!is.nan(sigma)) %>%
    filter(TARGET == target) %>%
    group_by(Lvals, Cvals) %>%
    summarize(sigma = mean(sigma))
    
  
  Lvals <- edge.stats$Lvals 
  Cvals <- edge.stats$Cvals
  
  sigma <- edge.stats$sigma
  
  m1 <- mle2(sigma ~ dnorm(mean = (a*Lvals) + (b*Cvals) + c*Lvals*Cvals, sd = 1),
             start = list(a = 1,b = 1,c = 1),
             data = data.frame(Lvals = Lvals, Cvals = Cvals, sigma = sigma))
  
  print((var(sigma) - var(residuals(m1)))/var(sigma))
  
  edge.stats.predict       <- edge.stats
  edge.stats.predict$sigma <- as.numeric(predict(m1))
  
    # Plot Edge Statistic
  
  fig <- ggplot(data = edge.stats, aes(x = Cvals , y = sigma, colour = as.factor(Lvals))) +
    geom_point() + 
    geom_line(data = edge.stats.predict, aes(x = Cvals, y = sigma,  colour = as.factor(Lvals))) +
    theme(aspect.ratio = 1)# +
  
  plot(fig)
  
  ggsave(file = paste0('~/Dropbox/Calen/Dropbox/edge_', target, '_separable.pdf'), fig, width = 40, height = 40)
}

# Fit a function to the standard deviation of the luminance responses.
fit.lum.stats <- function(target) {
  library(bbmle)
  
  lum.stats <- import_luminance_stats() %>%
    filter(!is.nan(sigma)) %>%
    filter(TARGET == target) %>%
    group_by(Lvals, Cvals) %>%
    summarize(sigma = mean(sigma))
  
  
  Lvals <- lum.stats$Lvals 
  Cvals <- lum.stats$Cvals
  
  sigma <- lum.stats$sigma
  
  m1 <- mle2(sigma ~ dnorm(mean = (a*Lvals), sd = 1),
             start = list(a = 1),
             data = data.frame(Lvals = Lvals, Cvals = Cvals, sigma = sigma))
  
  print((var(sigma) - var(residuals(m1)))/var(sigma))
  
  lum.stats.predict       <- lum.stats
  lum.stats.predict$sigma <- as.numeric(predict(m1))
  
  # Plot Edge Statistic
  
  fig <- ggplot(data = lum.stats, aes(x = Cvals , y = sigma, colour = as.factor(Lvals))) +
    geom_point() + 
    geom_line(data = lum.stats.predict, aes(x = Cvals, y = sigma,  colour = as.factor(Lvals))) +
    theme(aspect.ratio = 1)# +
  
  plot(fig)
  
  ggsave(file = paste0('~/Dropbox/Calen/Dropbox/edge_', target, '_separable.pdf'), fig, width = 40, height = 40)
}

