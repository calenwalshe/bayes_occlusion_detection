#' Plot formatted and publication ready thresholds.
#'
#' @return
#' @export
#'
#' @examples
plot_eccentricity_thresholds <- function() {
  library(ggthemes)
  summarize <- dplyr::summarise
  
  # Load and Source
  source(
    "~/Dropbox/Calen/Work/occluding/occlusion_detect/_human/human_psychometrics_tools.R"
  )
  source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R")
  load(
    "~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.psychometrics.scaled.rdata"
  )
  load(
    "~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/human.psychometrics.rdata"
  )
  #
  
  model.psychometrics.scaled <- model.psychometrics.scaled %>% filter(observer %in% c("optimal"))
  
  model.thresholds  <- get_threshold_extrap(model.psychometrics.scaled)
  human.thresholds  <- get_threshold(human.psychometrics)
  
  human.thresholds <- human.thresholds %>% ungroup() %>%
    mutate(TARGET = factor(TARGET, levels = c("vertical", "horizontal", "bowtie", "spot"))) %>%
    group_by(TARGET, BIN) %>% summarize(threshold = mean(threshold), se = sqrt(mean(sd^2))) %>%
    mutate(observer = "ave") %>%
    as_tibble() %>% 
    mutate(bExtrap = FALSE)
  
  model.thresholds <-
    model.thresholds %>% mutate(TARGET = factor(TARGET, levels = c(
      "vertical", "horizontal", "bowtie", "spot"
    )))
  
  
  threshold_values <- full_join(model.thresholds, human.thresholds)
  
  bin_values <- get_experiment_bin_values()
  d.1        <- merge(threshold_values, bin_values) %>%
                  arrange(TARGET, BIN, observer)
  
  d.1$linetype <-
    as.factor(ifelse(d.1$observer %in% c("ave", "rcw", "sps", "yhb"), "1", "2"))
  
  observer <- unique(d.1$observer)
  
  d.1$linetype <- as.factor(ifelse(d.1$observer %in% c("rcw", "sps", "yhb", "ave"), "1", "2"))
  
  d.1 <- d.1 %>% filter(observer %in% c("ave", "rcw", "sps", "yhb", "1 0 0", "0 1 0", "0 0 1", "optimal"))
  
  colour.vals <- RColorBrewer::brewer.pal(10, "Spectral")
  label.set <- data.frame(names = c("ave", "rcw", "sps", "yhb", "optimal", "mahalanobis",
                                    "nocov", "1 0 0", "0 1 0", "0 0 1"), 
                          pretty.names = c("Average Human", "rcw", "sps", "yhb", "Optimal", 
                                           "Mahalanobis", "No Covariance", "Edge", "Luminance", "Pattern"),
                          colours = as.character(rev(colour.vals[1:10])),
                          linetype = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2))
  
  label.set$colours <- as.character(label.set$colours)
  label.set$colours[1] <- "#3288bd"
  label.set$colours[5] <- "#fc8d59"
  
  my.labels <- (label.set[label.set$names %in% d.1$observer,])
  
  d.1$observer <- factor(d.1$observer, levels = my.labels[,1], labels = my.labels[,2])
  
  ### PLOTTING ###
  
  axis.text.sz <- 8
  axis.title.sz <- 12
  tick.sz       <- 1
  axis.ticks.length.sz = .1
  legend.key.size <- 15
  legend.text.size <- 10
  panel.border.sz  <- 2
  strip.text.sz    <- 12
  legend.line.sz   <- 1
  width.fig  <- 7
  height.fig <- 7
  
  pt.sz <- 2
  line.sz <- .5
  
  theme_set(theme_base())
  theme_update(aspect.ratio = 1, 
               axis.title = element_text(size = axis.title.sz),
               plot.background = element_blank(),
               axis.ticks = element_line(size = tick.sz),
               axis.ticks.length = unit(axis.ticks.length.sz, "cm"),
               panel.background = element_rect(fill = "white"),
               #axis.ticks.length=unit(tick.length, "cm"),
               #axis.text.x = element_text(margin=unit(axis.text.shift.x, "cm")), 
               #axis.text.y = element_text(margin=unit(axis.text.shift.y, "cm")),
               axis.text = element_text(size = axis.text.sz),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               legend.title=element_blank(),
               legend.text = element_text(size = legend.text.size),
               legend.key.size = unit(legend.key.size, "pt"),
               panel.border = element_rect(fill=NA, colour = "#19181A", size=panel.border.sz),
               plot.title = element_text(hjust = 0.5),
               #legend.key=element_blank(),
               strip.text = element_text(size = strip.text.sz),
               legend.key.width = unit(legend.line.sz,"cm"))
  
  # Luminance
  breaks.luminance <- round(unique(d.1[d.1$statType == "Lvals",]$statValue/100),2)
  error.width <- (max(breaks.luminance) - min(breaks.luminance))/25
  
  t.1.plot <- d.1 %>% filter(statType == "Lvals") %>% mutate(statValue = statValue / 100) %>%
    ggplot(data = ., aes(x = statValue, y = threshold,
                         colour = interaction(observer, linetype), linetype = interaction(observer,linetype))) +
    geom_point(size = pt.sz, aes(shape = bExtrap)) +
    geom_line(size = line.sz) +
    facet_wrap(~ TARGET) +
    xlab("Background Luminance (Rel)") +
    ylab("Eccentricity Threshold (\U00B0)") +
    geom_errorbar(data = . %>% filter(observer %in% c("Average Human", "rcw", "sps", "esm")), aes(ymin=threshold-se, ymax=threshold+se), width = error.width) +
    scale_linetype_manual("", labels = my.labels[,2], values = my.labels[,4]) +
    scale_colour_manual("", labels = my.labels[,2], values= as.character(my.labels[,3])) +
    scale_x_continuous() +
    scale_y_continuous(breaks = c(10, 20, 30, 40, 50)) +
    expand_limits(y = c(10, 50)) +
    coord_cartesian(ylim = c(0, 55)) +
    guides(shape = FALSE)
  
  plot(t.1.plot)  
  ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/eccentricity_thresholds/_figures/threshold_luminance.pdf", 
         t.1.plot, width = width.fig, height = height.fig, useDingbats=FALSE)
  
  # Contrast
  breaks.contrast <- round(d.1[d.1$statType == "Cvals",]$statValue %>% unique(),2)*100
  error.width <- (max(breaks.contrast) - min(breaks.contrast))/25
  
  t.1.plot <- d.1 %>% filter(statType == "Cvals") %>%
    ggplot(data = ., aes(x = statValue*100, y = threshold,
                         colour = interaction(observer, linetype), linetype = interaction(observer,linetype))) +
    geom_point(size = pt.sz, aes(shape = bExtrap)) +
    geom_line(size = line.sz) +
    facet_wrap(~ TARGET) +
    expand_limits(x = c(0,1), y = c(0, 30)) + 
    xlab("Background Contrast (RMS)") +
    ylab("Eccentricity Threshold (\U00B0)") +
    geom_errorbar(data = . %>% filter(observer %in% c("Average Human", "rcw", "sps", "esm")), aes(ymin=threshold-se, ymax=threshold+se), width = error.width) +
    scale_linetype_manual("", labels = my.labels[,2], values = my.labels[,4]) +
    scale_colour_manual("", labels = my.labels[,2], values= as.character(my.labels[,3])) +
    scale_y_continuous(breaks = seq(0, 32, 8)) +
    guides(shape = FALSE)
  
  plot(t.1.plot)  
  ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/eccentricity_thresholds/_figures/threshold_contrast.pdf", t.1.plot, width = width.fig, height = height.fig, useDingbats=FALSE)
  
  # Similarity
  breaks.similarity <- round(d.1[d.1$statType == "Svals",]$statValue %>% unique(),2)
  error.width <- (max(breaks.similarity) - min(breaks.similarity))/25
  t.1.plot <- d.1 %>% filter(statType == "Svals") %>%
    ggplot(data = ., aes(x = statValue, y = threshold,
                         colour = interaction(observer, linetype), linetype = interaction(observer,linetype))) +
    geom_point(size = pt.sz, aes(shape = bExtrap)) +
    geom_line(size = line.sz) +
    facet_wrap(~ TARGET) +
    expand_limits(x = c(.45,1), y = c(0, 25)) + 
    xlab("Background Similarity") +
    ylab("Eccentricity Threshold (\U00B0)") +
    geom_errorbar(aes(ymin=threshold-se, ymax=threshold+se), width = error.width) +
    scale_linetype_manual("", labels = my.labels[,2], values = my.labels[,4]) +
    scale_colour_manual("", labels = my.labels[,2], values= as.character(my.labels[,3])) +
    scale_x_continuous(breaks = seq(.45, .95, .15)) +
    guides(shape = FALSE)
  
  plot(t.1.plot)  
  ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/eccentricity_thresholds/_figures/threshold_similarity.pdf", t.1.plot, width = width.fig, height = height.fig, useDingbats=FALSE)
}

