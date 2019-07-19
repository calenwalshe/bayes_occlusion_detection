### PLOTTING ###

library(RColorBrewer)

source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R")

load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/efficiency.rdata")

model.efficiency$linetype <-
  as.factor(ifelse(model.efficiency$observer %in% c("ave", "rcw", "sps", "yhb"), "1", "2"))

observer <- unique(model.efficiency$observer)

model.efficiency$linetype <- as.factor(ifelse(model.efficiency$observer %in% c("rcw", "sps", "yhb", "ave"), "1", "2"))

model.efficiency <- model.efficiency %>% filter(observer %in% c("ave", "rcw", "sps", "yhb", "1 0 0", "0 1 0", "0 0 1", "optimal"))

colour.vals <- RColorBrewer::brewer.pal(10, "Spectral")
label.set <- data.frame(names = c("ave", "rcw", "sps", "yhb", "optimal", "mahalanobis",
                                  "nocov", "1 0 0", "0 1 0", "0 0 1"), 
                        pretty.names = c("Average Human", "rcw", "sps", "yhb", "Optimal", 
                                         "Mahalanobis", "No Covariance", "Edge", "Luminance", "Pattern"),
                        colours = as.character(colour.vals[1:10]),
                        linetype = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2))

my.labels <- (label.set[label.set$names %in% model.efficiency$observer,])

model.efficiency$observer <- factor(model.efficiency$observer, levels = my.labels[,1], labels = my.labels[,2])

bin.values <- get_experiment_bin_values()

model.efficiency <- merge(model.efficiency, bin.values)

model.efficiency$efficiency <- log2(model.efficiency$efficiency)

colour.vals <- RColorBrewer::brewer.pal(8, "Spectral")

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

pt.sz <- 1.5
line.sz <- .75

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
breaks.luminance <- round(model.efficiency[model.efficiency$statType == "Lvals",]$statValue %>% unique(),0)
t.1.plot <- model.efficiency %>% filter(statType == "Lvals") %>%
  ggplot(data = ., aes(x = statValue, y = efficiency,
                       colour = observer)) +
  geom_point(size = pt.sz) +
  geom_line(size = line.sz) +
  facet_wrap(~ TARGET) +
  xlab("Background Luminance (%)") +
  ylab("Log Efficiency") +
  #scale_linetype_manual("", labels = c("Optimal", "No Covariance", "Mahalanobis"), values = rep(c("solid", "dashed"), each = 3)) +
  scale_colour_manual(values=as.character(my.labels[,3])) +
  scale_x_continuous(breaks = breaks.luminance) +
  #scale_y_continuous(breaks = seq(0, 1, .1)) +
  expand_limits(y = 0)

plot(t.1.plot)
ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/efficiency/_figures/efficiency_luminance.pdf", t.1.plot, width = width.fig, height = height.fig, useDingbats=FALSE)

# Contrast
breaks.contrast <- round(model.efficiency[model.efficiency$statType == "Cvals",]$statValue %>% unique(),0)
t.1.plot <- model.efficiency %>% filter(statType == "Cvals") %>%
  ggplot(data = ., aes(x = statValue, y = efficiency,
                       colour = observer)) +
  geom_point(size = pt.sz) +
  geom_line(size = line.sz) +
  facet_wrap(~ TARGET) +
  xlab("Background Contrast (%)") +
  ylab("Log Efficiency") +
  #scale_linetype_manual("", labels = c("Optimal", "No Covariance", "Mahalanobis"), values = rep(c("solid", "dashed"), each = 3)) +
  scale_colour_manual(values= as.character(my.labels[,3])) +
  scale_x_continuous(breaks = breaks.contrast) +
  #scale_y_continuous(breaks = c(0, .1, .2, .3)) +
  expand_limits(y = 0)

plot(t.1.plot)  
ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/efficiency/_figures/efficiency_contrast.pdf", t.1.plot, width = width.fig, height = height.fig, useDingbats=FALSE)

# Similarity
breaks.similarity <- round(model.efficiency[model.efficiency$statType == "Svals",]$statValue %>% unique(),2)
t.1.plot <- model.efficiency %>% filter(statType == "Svals") %>%
  ggplot(data = ., aes(x = statValue, y = efficiency,
                       colour = observer)) +
  geom_point(size = pt.sz) +
  geom_line(size = line.sz) +
  facet_wrap(~ TARGET) +
  xlab("Similarity") +
  ylab("Log Efficiency") +
  #scale_linetype_manual("", labels = c("Optimal", "No Covariance", "Mahalanobis"), values = rep(c("solid", "dashed"), each = 3)) +
  scale_colour_manual(values = as.character(my.labels[,3])) +
  scale_x_continuous(breaks = seq(.4, 1, .1)) +
  #scale_y_continuous(breaks = seq(0, .5, .1)) +
  expand_limits(y = 0)

plot(t.1.plot)  
ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/efficiency/_figures/efficiency_similarity.pdf", t.1.plot, width = width.fig, height = height.fig, useDingbats=FALSE)
