#' Publication ready eccentricity psychometric function
#' @return
#' @export
#'
#' @examples
plot_eccentricity_psychometric <- function() {
  library(ggthemes)
  load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/human.psychometrics.rdata")
  load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.psychometrics.scaled.rdata")
  
  model.psychometrics.scaled$data <- map(model.psychometrics.scaled$data, function(x) {
    data <- x %>% select(eccentricity, dprime) %>% mutate(pc = pnorm(dprime/2))
  })
  
  model.psychometrics.scaled.1 <- model.psychometrics.scaled %>% 
    select(TARGET, BIN, e0, data, b, d0, observer) %>% 
    mutate(gamma = 0) %>%
    filter(BIN == 3, observer == "optimal")
  
  human.psychometrics.1        <- human.psychometrics %>% select(TARGET, BIN, e0, data, b, d0, observer, gamma)
  
  combined.psychometrics <- rbind(model.psychometrics.scaled.1, human.psychometrics.1)
  
  psychometric.values <- combined.psychometrics %>% 
    group_by(TARGET, observer, BIN) %>%
    nest()
  
  psychometric.values$pc_hat <- map(psychometric.values$data, function(parameters) {
    d0    <- parameters$d0
    e0    <- parameters$e0
    b     <- parameters$b
    gamma <- parameters$gamma
    
    x <- seq(0, 60, .01)
    
    f <- function(x, d0, e0, b, gamma) {
      dprime <- d0 * e0^b / (e0^b + x^b)
      pc     <- (pnorm(dprime/2 - gamma) + (1 - pnorm(-dprime/2 - gamma)))/2
    } 
    
    f.1 <- Vectorize(f, vectorize.args = "x")
    
    data.frame(eccentricity  = x, percent_correct = f.1(x, d0, e0, b, gamma))
    
  })
  
  psychometric.line       <- psychometric.values %>%
    unnest(pc_hat) %>%
    filter(BIN == 3) %>%
    filter(observer %in% c("rcw", "sps", "yhb", "optimal"))
  
  psychometric.empirical  <- combined.psychometrics %>% 
    unnest(data) %>% 
    select(TARGET, observer, BIN, eccentricity, pc) %>% 
    rename(percent_correct = pc) %>% 
    filter(BIN == 3) %>% 
    filter(observer %in% c("rcw", "sps", "yhb", "optimal"))
  
  psychometric.line$linetype <- as.factor(ifelse(psychometric.line$observer %in% c("rcw", "sps", "yhb"), "1", "2"))
  psychometric.empirical$linetype <- as.factor(ifelse(psychometric.empirical$observer %in% c("rcw", "sps", "yhb"), "1", "2"))
  
  psychometric.empirical$observer <- factor(psychometric.empirical$observer, levels = c("rcw", "sps", "yhb","optimal"), labels = c("rcw", "sps", "yhb", "Model Observer"))
  psychometric.line$observer <- factor(psychometric.line$observer, levels = c("rcw", "sps", "yhb","optimal"), labels = c("rcw", "sps", "yhb", "Model Observer"))

  # PLOT #
  colour.vals <- RColorBrewer::brewer.pal(4, "GnBu")
  
  plot.colours   <- colour.vals[2:4]
  plot.colours[4] <- "#fc8d59"
  
  plot.colours <- c("#074684", "#0ea5c6", "#a0edf7", "#fc8d59")
  
  axis.text.sz <- 8
  axis.title.sz <- 12
  tick.sz       <- 1
  axis.ticks.length.sz = .075
  legend.key.size <- 15
  legend.text.size <- 10
  panel.border.sz  <- 1
  strip.text.sz    <- 12
  legend.line.sz   <- 1
  width.fig  <- 4
  height.fig <- 4
  
  pt.sz <- 1.25
  line.sz <- .25
  
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
               #legend.key.size = unit(legend.key.size, "pt"),
               panel.border = element_rect(fill=NA, colour = "#19181A", size=panel.border.sz),
               plot.title = element_text(hjust = 0.5),
               #legend.key=element_blank(),
               strip.text = element_text(size = strip.text.sz),
               legend.key.width = unit(legend.line.sz,"cm"))
  

  psychometric.line.1 <- psychometric.line %>% filter(observer != "Model Observer")
  
  base.plot.with.legend <- ggplot(psychometric.empirical, aes(x = eccentricity, y = percent_correct, colour = observer)) +
    geom_line(data = psychometric.line.1, aes(x = eccentricity, y = percent_correct, colour = observer), size = line.sz) +
    geom_point(size = pt.sz) + 
    facet_wrap(~ TARGET) +
    #scale_linetype_manual("", labels = c("rcw", "sps", "esm"), values = c(1,1,1,2,2,2)) +
    scale_colour_manual("", values=rev(plot.colours)) + 
    xlab("Eccentricity (\U00B0)") + ylab("Percent Correct") + 
    #scale_x_continuous(breaks = seq(3, 23, length.out = 5)) +
    expand_limits(x = c(0, 62)) +
    scale_x_continuous(breaks = c(0, 30, 60)) +
    scale_y_continuous(breaks = c(.5, .75, 1))
    
  ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/eccentricity_psychometric/_figures/legend.pdf", plot = base.plot.with.legend, device = "pdf", width = width.fig, height = height.fig, useDingbats = FALSE)
  
  base.plot.no.legend <- base.plot.with.legend + guides(linetype = FALSE, colour= FALSE)
  ggsave(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/eccentricity_psychometric/_figures/eccentricity_psychometric.pdf", plot = base.plot.no.legend, device = "pdf", width = width.fig, height = height.fig, useDingbats = FALSE)
}

