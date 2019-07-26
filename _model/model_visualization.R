plot.model.psychometric <- function(model.psychometrics, path = "~/Dropbox/Calen/Dropbox/model.psychometrics.pdf") {
  library(ggplot2)

  dprime.vals <- map(model.psychometrics$data, function(x) data.frame(eccentricity = x$eccentricity, dprime = x$dprime))
  
  
  model.psychometrics$response <- dprime.vals
  
  map.observer <- model.psychometrics %>% group_by(TARGET, BIN, observer) %>% nest()
  
  psychometric.dat <- map(map.observer$data, function(x) {
    dprime_hat <- x$d0 * x$e0^x$b / (x$e0^x$b + seq(0, 26, .1)^x$b)
    
    data.frame(eccentricity = seq(0, 26, .1), dprime = dprime_hat)
  })
  
  model.psychometrics$psychometric.dat <- psychometric.dat
  
  psychometric.dat.1 <- model.psychometrics %>% select(-response) %>% unnest(psychometric.dat)
  psychometric.dat.2 <- model.psychometrics %>% select(-psychometric.dat) %>% unnest(response)
  
  fig.1 <- ggplot(data = psychometric.dat.1, aes(x = eccentricity, y = dprime, colour = observer)) + 
    geom_line() + 
    geom_point(data = psychometric.dat.2, aes(x = eccentricity, y = dprime, colour = observer)) + 
    facet_wrap(TARGET ~ BIN, ncol = 15, scale = "free") + 
    theme_set(theme_gray(base_size = 18)) +
    theme(aspect.ratio = 1)
  
  ggsave(plot = fig.1, filename = path, device = "pdf", width = 40, height = 30)
}

plot.correlation.observers <- function(observer.dprime) {
  
  obs.1 <- observer.dprime %>% select(BIN, TARGET, observer, eccentricity, dprime) %>% spread(observer, dprime)
  
  obs.1[,4] <- scale(obs.1[,4])
  obs.1[,6] <- scale(obs.1[,6])
  
  fig.10 <- ggplot(obs.1, aes_string(x = "optimal", y = "mahalanobis", colour = "TARGET")) + 
    geom_point(size = point.sz) +
    labs(colour = "") +
    xlab(expression("d'"["norm"]~~"(Mahalanobis)")) +
    ylab(expression("d'"["norm"]~~"(Optimal)")) +
    theme.1 +
    scale_colour_manual(values = colours.targets)
    
  ggsave(fig.10, filename = '~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/cue_correlation.pdf', scale = 2)
}

plot.correlation.optimal.human <- function() {
  library('latex2exp')
  library(dplyr)
  
  summarize = dplyr::summarise
  
  load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/human.psychometrics.rdata")
  load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.psychometrics.scaled.rdata")
  
  human.threshold <- get_threshold(human.psychometrics)
  model.threshold <- get_threshold(model.psychometrics.scaled)
  
  human.threshold <- human.threshold[, c("TARGET", "BIN", "observer", "threshold")] 
  human.threshold <- human.threshold %>% group_by(BIN, TARGET) %>% summarize(threshold = mean(threshold)) %>% mutate(observer = "human")
  
  model.threshold <- model.threshold[, c("TARGET", "BIN", "observer", "threshold")] 
  
  threshold.1 <- rbind(human.threshold %>% data.frame(), model.threshold %>% data.frame())
  
  threshold.2 <- threshold.1 %>% filter(observer %in% c("human", "optimal")) %>% group_by(BIN, TARGET) %>% spread(observer, threshold)
  
  threshold.3 <- threshold.2
  
  threshold.3[,3] <- scale(threshold.2[,3])
  threshold.3[,4] <- scale(threshold.2[,4])
  
  threshold.3$outlier <- as.factor(ifelse((threshold.3[,3] > 2) | (threshold.3[,4] > 2), 1, 0))

  threshold.3$observer <- factor(threshold.3$observer, levels = c("optimal", "nocov", "mahalanobis"), labels = c("Optimal", "Mahalanobis", "No Covariance"))
  
  x.lab <- TeX("$Human Eccentricity (\\degree)$")
  fig.2 <- ggplot(threshold.3, aes_string(x = "human", y = "optimal", colour = "TARGET")) + 
    geom_point(size = point.sz) +
    labs(colour = "") +
    xlab(x.lab) +
    theme.1 +
    xlab("Eccentricity Threshold (z)") +
    ylab("Eccentricity Threshold (z)") +
    theme.1 +
    scale_colour_manual(values = colours.targets) + 
    guides(shape = FALSE) +
    geom_smooth(method = "lm", se = F, size = 2,colour = colour.regression) +
    expand_limits(x = c(-2,4), y = c(-2,4))
  
  plot(fig.2)
  
  ggsave(filename = '~/Dropbox/Calen/Work/occluding/paper/figures/correlations/model_human_correlation.pdf', plot = fig.2, scale = 3, useDingbats = FALSE)
}

plot.quick.thresholds <- function(model.psychometrics) {
  bin.values <- get_experiment_bin_values()
  model.thresholds <- get_threshold(model.psychometrics)
  
  plot.values <- merge(bin.values, model.thresholds) %>% as_tibble()
  theme_set(theme_bw())  
  ggplot(plot.values, aes(x = statValue, y = threshold, colour = observer)) + geom_point() + geom_line() + facet_wrap(statType~TARGET, scales = "free") + theme(aspect.ratio = 1) + coord_cartesian(ylim = c(0, 15))
}

plot.quick.dprime <- function(model.error) {
  source('~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R')
  #load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/old_model_error/model.error.rdata")
  load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.error2d.rdata")
  bin.values <- get_experiment_bin_values()
  mod.marg <- by_row(model.error, function(x) {
    themeans <- abs(x$data_mean_vec_1[[1]] - x$data_mean_vec_0[[1]])
    names(themeans) <- c("edge", "pattern")
    thesigma <- sqrt(x$data_sd_vec_0[[1]]^2)
    
    thedprime <- themeans / thesigma
    
    data.frame(edge = thedprime[1], pattern = thedprime[2])
  }, .to = "output") %>% select(BIN, TARGET, eccentricity, output)
  
  mod.marg.1 <- mod.marg %>% unnest(output)
  
  mod.marg.2 <- merge(mod.marg.1, bin.values)
  
  mod.marg.3 <- mod.marg.2 %>% filter(statType == "Svals")#, round(eccentricity) %in% c(0, 11))
  
  model.error.1 <- merge(model.error, bin.values) %>% filter(observer == "optimal") %>% as_tibble() %>% filter(statType == "Cvals")
  
  theme_set(theme_bw())
  fig.1 <- ggplot(mod.marg.3, aes(x = statValue, y = sqrt(pattern^2 + edge^2), colour = as.factor(eccentricity))) + geom_point() + geom_line() + facet_grid(~TARGET) + theme(aspect.ratio = 1) + xlab("Luminance") + ylab("dprime")
  
  fig.2 <- ggplot(mod.marg.3, aes(x = statValue, y = edge, colour = as.factor(eccentricity))) + geom_point() + geom_line() + facet_grid(~TARGET) + theme(aspect.ratio = 1) + xlab("Luminance") + ylab("dprime")
  
  ggsave(file = "~/Dropbox/Calen/Dropbox/dprime.combined.new.pdf", fig.1, scale = 3)
  ggsave(file = "~/Dropbox/Calen/Dropbox/dprime.mean.new.pdf", fig.2, scale = 3)
  
}