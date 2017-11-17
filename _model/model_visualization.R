plot_all_single_eccentricity <- function(human.psychometrics, 
    model.dprime) {
    
    sub_params <- human.psychometrics
    
    sub_params <- merge(sub_params, data.frame(eccentricity = unique(model.dprime$eccentricity)))
    
    pc <- lapply(1:nrow(sub_params), FUN = function(x) pc = pnorm(sub_params[x, 
        "d0"]/2 * (sub_params[x, "e0"]^sub_params[x, "b"])/(sub_params[x, 
        "e0"]^sub_params[x, "b"] + sub_params[x, "eccentricity"]^sub_params[x, 
        "b"]))) %>% do.call(rbind, .)
    
    sub_params$pc <- pc
    
    human_dat <- sub_params %>% group_by(TARGET, BIN, eccentricity) %>% 
        summarize(pc = mean(pc)) %>% ungroup() %>% mutate(dprime = 2 * 
        as.numeric(qnorm(pc)), SUBJECT = "human") %>% data.frame
    
    model_dat <- model.dprime %>% select(TARGET, BIN, 
        eccentricity, percent_correct, dprime) %>% rowwise() %>% 
        mutate(pc = percent_correct, dprime = 2 * as.numeric(qnorm(pc)), 
            percent_correct = NULL, SUBJECT = "model") %>% data.frame
    
    dat <- rbind(human_dat, model_dat) %>% group_by(BIN, TARGET, 
        eccentricity) %>% summarize(pc_diff = pc[1] - pc[2])
    
    fig <- ggplot(dat, aes(x = eccentricity, y = pc_diff)) + 
        geom_point() + facet_grid(~TARGET)
    
    ggsave(last_plot(), file = "~/Dropbox/Calen/Dropbox/compare_pc.pdf")
    print(1)
}

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

plot_histograms <- function(template.response) {
  library(ggplot2)
  library(dplyr)
  lapply(unique(template.response$eccentricity), FUN = 
           function(x) {
             p1 <- ggplot(data = template.response %>% filter(eccentricity == x), aes(x = TRESP, fill = factor(TPRESENT))) + geom_histogram(binwidth = .00001) + facet_wrap(BIN ~ TARGET + TPRESENT, ncol = 4, scales = "free_x") + theme(aspect.ratio = 1)
             ggsave(file = paste0('~/Dropbox/Calen/Dropbox/model_histogram_', round(x,3), '.pdf'), p1, height = 40,units = "in")
  }
  )
}

plot_dprime <- function(dprime1, dprime2) {
  library(ggplot2)
  dprime1 <- dprime1 %>% data.frame %>% 
    select(BIN, TARGET, eccentricity, SUBJECT, dprime) %>%
    unique()
  
  dprime2 <- dprime2 %>% data.frame %>% 
    select(BIN, TARGET, eccentricity, SUBJECT, dprime) %>%
    unique()
  
  dprime_combined <- rbind(dprime1, dprime2)
  
  p1 <- ggplot(data = dprime_combined, aes(x = eccentricity, y = dprime, colour = SUBJECT)) + geom_point() + facet_grid(BIN ~ TARGET, scales = "free_y")
  
  ggsave(file = '~/Dropbox/Calen/Dropbox/figure1.pdf', p1, width = 30, height = 30)
  
}

plot_dprime_family <- function(model.dprime) {
  library(ggplot2)
  p1 <- ggplot(data = model.dprime, aes(x = (SUBJECT), y = (dprime), colour = as.factor(eccentricity))) +
    geom_point() +
    facet_grid(BIN ~ TARGET, scales = "free_y")
  ggsave(file = "~/Dropbox/Calen/Dropbox/fig.pdf", p1, width = 30, height = 30, units = "in")
  
}