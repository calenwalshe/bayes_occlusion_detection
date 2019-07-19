#' Visualize a single psychometric function
plot_psychometric <- function(psychometrics, out_path = "~/Dropbox/Calen/Dropbox/") {
  library(dplyr)
  library(ggplot2)
  library(purrrlyr)
  library(purrr)
  
  psychometrics <- psychometrics %>% arrange(observer, BIN, TARGET)
  
  n_row <- nrow(psychometrics)
  
  plot.fig.data <- lapply(1:n_row, FUN = function(x) {
    data <- psychometrics[x, ]
    
    # model fit
    eccentricity <- seq(0, 25, .1)
    d0 <- data$d0
    e0 <- data$e0
    b  <- data$b
    dprime <- d0 * e0^b/(e0^b + eccentricity^b)
    
    model <- data.frame(TARGET = data$TARGET, BIN = data$BIN, observer = data$observer, eccentricity = eccentricity, dprime = dprime, data = factor("psychometric"))
    
    
    # raw data
    raw.data <- psychometrics[[x, "data"]]
    eccentricity <- raw.data$eccentricity
    dprime       <- raw.data$dprime
    
    response.dat <- data.frame(TARGET = data$TARGET, BIN = data$BIN, observer = data$observer, eccentricity = eccentricity, dprime = dprime, data = factor("response"))
    
    all.dat <- rbind(model, response.dat)
    
  })
  
  plot.data <- do.call(rbind, plot.fig.data) %>% as_tibble()
  
  fig <- ggplot() + 
    geom_point(data = plot.data %>% filter(data == "response"), 
               aes(x = eccentricity, y = dprime, colour = observer)) +
    geom_line(data = plot.data %>% filter(data == "psychometric"), 
              aes(x = eccentricity, y = dprime, colour = observer)) +
    facet_grid(TARGET ~ BIN) +
    theme(aspect.ratio = 1)
  
  ggsave(filename = '~/Dropbox/Calen/Dropbox/psychometrics.pdf', plot = fig, device = 'pdf', width = 40, height = 30, units = "in")
    
}
