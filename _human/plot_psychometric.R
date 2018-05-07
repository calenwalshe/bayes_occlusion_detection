#' Visualize a single psychometric function
plot_psychometric <- function(psychometric.params, empirical.obs, bin = 1, target = "vertical", out_path = "~/Dropbox/Calen/Dropbox/") {
  library(dplyr)
  library(ggplot2)
  library(purrrlyr)
  library(purrr)
  
  
  psychometric.params$e0 <- as.numeric(psychometric.params$e0)
  
  empirical.dat <- empirical.obs %>%
    group_by(BIN, TARGET, SUBJECT) %>%
    nest() %>%
    mutate(human.dat = map(data, function(data) {
      data.frame(dprime = data$dprime, pc = data$percent_correct, eccentricity = data$eccentricity)
    })) %>%
    select(BIN, TARGET, SUBJECT, human.dat)

  dat <- merge(psychometric.params, empirical.dat) %>%
    as_tibble()

  all.dat <- dat %>%
    as_tibble() %>%
    group_by(SUBJECT, TARGET, BIN) %>%
    nest() %>%
    mutate(psychometric.dat = map(data, function(data) {
      
      eccentricity <- seq(0, 25, .1)
      d0 <- data$d0
      e0 <- data$e0
      b  <- data$b
      gamma  <- data$gamma
      dprime <- d0 * e0^b/(e0^b + eccentricity^b) - gamma
      percent_correct     <- pnorm(1/2 * dprime)
      
      psychometric.dat <- data.frame(eccentricity = eccentricity, percent_correct = percent_correct, dprime = dprime)
    })) %>%
    unnest(data) 
  
    by_row(all.dat, function(row) {
      human <- row$human.dat[[1]]
      model <- row$psychometric.dat[[1]]
      
      BIN <- row$BIN
      TARGET <- row$TARGET
      SUBJECT <- row$SUBJECT      
      
      fig <- ggplot(model, aes(x = eccentricity, y = percent_correct)) + 
        geom_line() +
        geom_point(data = human, inherit.aes = F, aes(x = eccentricity, y = pc)) +
        ggtitle(paste('Bin:', BIN, 'TARGET:', TARGET, 'SUBJECT:', SUBJECT))
      
      ggsave(plot = fig, filename = paste0('~/Dropbox/Calen/Dropbox/', SUBJECT, '-', TARGET,'-', BIN, '.pdf'), device = 'pdf')      
    })
}
