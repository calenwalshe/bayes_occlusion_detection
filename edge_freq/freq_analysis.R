freq.analysis <- function() {
  
  
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(purrrlyr)
  library(tidyr)
  library(broom)
  
  bin_values <- get_experiment_bin_values() %>%
    filter(TARGET == "vertical", S == 5, statType %in% c("Cvals", "Lvals")) %>%
    mutate(S = 10) %>%
    arrange(statType, L, C, S)
  
  bin_values$BIN <- c(1,2,3,4,5,6,7,3,8,9,10)
  
  
  path.dat <- '~/Dropbox/Calen/Dropbox/edge.txt'
  
  edge.freq <- read.csv(path.dat, sep = '\t')
  
  edge.freq.subset <- edge.freq %>%
    filter((SAMPLE %in% 1:1000 & TPRESENT == 1) | (SAMPLE %in% 1001:2000 & TPRESENT == 2))
  
  edge.freq.nest <- edge.freq.subset %>%
    group_by(SAMPLE, BIN, TARGET, PYRAMIDLVL, TPRESENT) %>%
    nest() %>%
    mutate(fft.edge = map(data, function(x) abs(fft((x$VALUE))))) %>%
    group_by(BIN, TARGET, PYRAMIDLVL, TPRESENT) %>%
    nest() %>%
    mutate(fft_mean = map(data, function(x) data.frame(freq = 1:length(x$fft.edge[[1]]),
                                                       amp.mean = colMeans(do.call(rbind, x$fft.edge)),
                                                       amp.var  = apply(do.call(rbind, x$fft.edge), 2, var)))) %>%
    unnest(fft_mean) %>%
    group_by(BIN, TARGET, PYRAMIDLVL, freq) %>%
    summarize(dprime = abs(diff(amp.mean)/sqrt(mean(amp.var))))
  
  edge.freq.all <- merge(edge.freq.nest, bin_values, by = "BIN") %>%
    mutate(TARGET = TARGET.y) %>%
    select(-TARGET.x, -TARGET.y) %>%
    filter(!statType == "Svals") %>%
    filter(freq < max(freq)/2)
  
  
  fig <- ggplot(edge.freq.all, aes(x = freq, y = dprime)) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(statType ~ statValue, ncol = 6) +
    theme(aspect.ratio = 1)
  
  ggsave(plot = fig,filename = '~/Dropbox/Calen/Work/freq_fig.pdf', width = 40, height = 15)
  
}