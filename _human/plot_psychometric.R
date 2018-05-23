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

# Plot a single bin. Pretty and ready for presentation.
plot.figure.bin <- function(human.psychometrics, human.detect) {
    human.dat <- human.detect %>% filter(BIN == 3)
    human.psychometrics$TARGET <- factor(human.psychometrics$TARGET, levels = c("vertical", "horizontal", "bowtie", "spot"))
    obs.dat   <- human.psychometrics %>%
      filter(BIN == 3) %>%
      group_by(TARGET, BIN, SUBJECT) %>%
      nest() %>%
      mutate(psy.obs = map(data, function(x) {
        e0 <- x$e0
        b  <- x$b
        gamma <- x$gamma
        d0    <- x$d0
        
        ecc <- seq(0, 23, .1)
        
        obs <- pnorm(1/2 * d0 * e0^b/(e0^b + ecc^b))
        
        data.frame(eccentricity = ecc, percent_correct = obs, threshold = x$threshold)
      })) %>%
      unnest(psy.obs)
    
    fig <- ggplot(data = obs.dat, aes(x = eccentricity, y = percent_correct, colour = SUBJECT)) + 
      geom_line(size = 1.5) + 
      geom_point(data = human.dat, aes(x = eccentricity, y = percent_correct), size = 2.25) +
      facet_wrap(~TARGET) +
      theme_bw() +
      theme(aspect.ratio = 1) + 
      theme_set(theme_bw(base_size = 45))  +# pre-set the bw theme.
      scale_color_brewer(name = "Subject", palette = "Dark2") +
      expand_limits(y = c(.5, 1)) +
      xlab("Eccentricity (º)") +
      ylab("Percent Correct")
    
    
    plot(fig)
    
    ggsave(file = '~/Dropbox/Calen/Work/occluding/detection_model_analysis/presentations/vss_2018/center_bin_psychometrics.pdf', fig, scale = 1.35, limitsize = FALSE)

}

