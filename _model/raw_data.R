# Visualize dprime across conditions.
plot_vis_dprime <- function(template.response) {
  library(dplyr)
  library(ggplot2)
  
  dprime.all <- template.response %>% group_by(eccentricity, BIN, TPRESENT, TARGET, function_name, statType, statValue) %>%
    summarize(var = var(TRESP), mu = mean(TRESP)) %>%
    arrange(eccentricity, BIN, TARGET, function_name, TPRESENT) %>%
    group_by(eccentricity, BIN, TARGET, function_name, statType, statValue) %>%
    summarize(sd_avg = sqrt((var[1] + var[2])/2), delta = abs(mu[1] - mu[2]), dprime = delta/sd_avg, pc = pnorm(dprime/2), logdprime = log10(dprime))
  
  
  mk.fig <- function(x) {
    fig <- ggplot(data = dprime.all %>% filter(eccentricity == x), aes(x = statValue, y = logdprime, colour = function_name)) + 
      geom_point() + 
      geom_line() +
      facet_wrap(statType~TARGET, nrow = 3, scales = "free_x") +
      ylim(-3, 4) +
      theme(aspect.ratio = 1)
    
    ggsave(filename = paste0("~/Dropbox/Calen/Dropbox/", x,'.pdf'), plot = fig, height = 40, width = 40, units = "cm")
  }
  lapply(unique(template.response$eccentricity), FUN = function(x) mk.fig(x))
}