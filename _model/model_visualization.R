plot_model_cormat <- function(file_path = '~/Dropbox/Calen/Work/data_storage/eccentricity_project/') {
  
  f <- function(file_path) {
    
    load(file_path)
    
    melt_cor <- model.data %>% 
      rowwise() %>% 
      mutate(melted_cor = list(data.frame(BIN = BIN, TARGET = TARGET, eccentricity, TPRESENT, unlist(LL))))
    
    melt_cor <- do.call(rbind, melt_cor$melted_cor)
  }
  
  files <- list.files(file_path, full.names = T)
  
  all_cor <- lapply(files, FUN = function(x) f(x)) %>% 
    do.call(rbind, .) %>%
    arrange(BIN, TARGET, eccentricity, TPRESENT)
  
  bin_values <- get_experiment_bin_values()
  
  all_cor <- merge(bin_values, all_cor)
  
  f <- function(x) {
    p.1 <- ggplot(data = all_cor %>% filter(eccentricity == x, TPRESENT == 0), aes(x = Var1, y = Var2, fill = (value))) + 
      geom_tile() + 
      facet_grid(BIN ~ TARGET) + 
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
      theme_minimal()
    
    return(p.1)
  }

  # create figures
  lapply(unique(all_cor$eccentricity),
         FUN = function(x) ggsave(filename = paste0('~/Dropbox/Calen/Dropbox/correlations/', round(x,3), '.pdf'),
                                  f(x), width = 40, height = 20, units = 'in'))
}

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
plot_efficiency <- function(efficiency.df) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  bin_values <- get_experiment_bin_values()
  
  m.eff.all <- merge(bin_values, efficiency.df)

  
  # Efficiency by statistic    
  fig.1 <- ggplot(m.eff.all, aes(x = statValue, y = log10(efficiency), colour = as.factor(eccentricity))) +
    facet_grid(~statType) +
    geom_point() + 
    geom_line() +
    #geom_hline(yintercept = average.eff$efficiency.mean) +
    facet_grid(TARGET ~ statType, scales = "free") +
    theme_gray()

  ggsave(file = "~/Dropbox/Calen/Dropbox/efficiency.fig.statistic.pdf", width = 10, height = 10 * 1200/1600)
  
  # Efficiency by eccentricity
  
  fig.dat <- m.eff.all %>%
    group_by(statType, TARGET) %>%
    nest() %>%
    by_row(., function(row) {

      fig.2 <- ggplot(row$data[[1]], aes(x = eccentricity, y = (efficiency), colour = as.factor(statValue))) +
        geom_point() + 
        theme_gray()
      
      ggsave(file = paste0("~/Dropbox/Calen/Dropbox/efficiency.fig.eccentricity", "-", row$TARGET[[1]], "-", row$statType[[1]], ".pdf"), width = 10, height = 10 * 1200/1600)
    })
    
  plot(fig.2)
  
  ggsave(file = "~/Dropbox/Calen/Dropbox/efficiency.fig.statistic.pdf", width = 10, height = 10 * 1200/1600)
  
}