pink_template_stat <- function(statType = 'tSigma', bPresent = 0, target = "vertical") {
  
  templateSigmaAbs  <- read.table('~/Dropbox/Calen/Dropbox/tMatchpink_Response_sigma.txt', header = T, sep = '\t') %>% mutate(PRESENT = 0)
  templateMeanAbs   <- read.table('~/Dropbox/Calen/Dropbox/tMatchpink_Response_mean.txt', header = T, sep = '\t') %>% mutate(PRESENT = 0)
  
  templateSigma     <- rbind(templateSigmaAbs, templateMeanAbs)
  
  templateAll       <- templateSigma %>% group_by(L,C,S,PYRAMIDLVL,TARGET) %>% arrange(L,C,S,TARGET, PRESENT)
  
  bin_values        <- get_bin_values()
  
  templateStats     <- templateAll %>% merge(., bin_values) %>% filter(TARGET_NAME == target) %>% filter(PRESENT == bPresent)
  
  fig               <- ggplot(data = templateStats, aes(x = Svals, y = tSigma, shape = TARGET_NAME, colour = as.factor(Cvals))) +
    geom_point() + facet_wrap(~Lvals, ncol = 3, scales = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig)
}

pink_edge_stat <- function(statType = 'tSigma', bPresent = 0, target = "vertical") {
  edgeSigmaAbs <- read.table('~/Dropbox/Calen/Dropbox/Eabspink_Response_sigma.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(0))
  edgeSigmaPres <- read.table('~/Dropbox/Calen/Dropbox/Eprespink_Response_sigma.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(1))
  
  edgeMeanAbs  <- read.table('~/Dropbox/Calen/Dropbox/Eabspink_Response_mean.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(0), tmu = tSigma, tSigma = NULL)
  edgeMeanPres <- read.table('~/Dropbox/Calen/Dropbox/Eprespink_Response_mean.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(1), tmu = tSigma, tSigma = NULL)
  
  edgeAbs  <- merge(edgeSigmaAbs, edgeMeanAbs)
  edgePres <- merge(edgeSigmaPres, edgeMeanPres)
  
  edgeAll  <- rbind(edgeAbs, edgePres) %>% group_by(L,C,S,PYRAMIDLVL,TARGET) %>% arrange(L,C,S,TARGET, PRESENT)
  
  bin_values         <- get_bin_values()
  
  edgeStatistics     <- edgeAll %>% merge(., bin_values) %>% filter(TARGET_NAME == target) %>% filter(PRESENT == bPresent)
  
  # All bins
  if(statType == "tmu") {
    fig <- ggplot(data = edgeStatistics, aes(x = Svals, y = tmu, shape = (TARGET_NAME), colour = as.factor(Cvals))) +
      geom_point() + 
      facet_wrap(~Lvals, ncol = 3) +
      theme(aspect.ratio = 1) +
      scale_color_brewer(palette="Spectral")
  }else{
    fig <- ggplot(data = edgeStatistics, aes(x = Svals, y = tSigma, shape = (TARGET_NAME), colour = as.factor(Cvals))) +
      geom_point() + 
      facet_wrap(~Lvals, ncol = 3) +
      theme(aspect.ratio = 1) +
      scale_color_brewer(palette="Spectral")
  }
  ggsave(fig, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'pink_',target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")  
  
  
  # Only Contrast
  contrastStats <- edgeAll %>% merge(., bin_values)
  if(statType == "tmu") {
    contrastStats <- contrastStats %>% group_by(Cvals, PRESENT, TARGET_NAME) %>% summarize(tmu = mean(tmu, na.rm = T))    
  }else{
    contrastStats <- contrastStats %>% group_by(Cvals, PRESENT, TARGET_NAME) %>% summarize(tSigma = mean(tSigma, na.rm = T))        
  }
  
  fig_contrast <- ggplot(data = contrastStats, aes_string(x = 'Cvals', y = statType, colour = "TARGET_NAME")) +
    geom_point() +
    facet_wrap(~PRESENT, scales = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig_contrast, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'pink_', 'con_', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig_contrast)
  
  # Only Luminance
  
  lumStats <- edgeAll %>% merge(., bin_values)
  
  if(statType == "tmu") {
    lumStats <- lumStats %>% group_by(Lvals, PRESENT, TARGET_NAME) %>% summarize(tmu = mean(tmu, na.rm = T))    
  }else{
    lumStats <- lumStats %>% group_by(Lvals, PRESENT, TARGET_NAME) %>% summarize(tSigma = mean(tSigma, na.rm = T))        
  }
  
  fig_lum <- ggplot(data = lumStats, aes_string(x = 'Lvals', y = statType, colour = "TARGET_NAME")) +
    geom_point() +
    facet_wrap(~PRESENT, scales = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig_lum, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'pink_', 'lum_', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig_lum)
  
  # Only similarity
  simStats <- edgeAll %>% merge(., bin_values)
  if(statType == "tmu") {
    simStats <- simStats %>% group_by(Svals, PRESENT, TARGET_NAME) %>% summarize(tmu = mean(tmu, na.rm = T))    
  }else{
    simStats <- simStats %>% group_by(Svals, PRESENT, TARGET_NAME) %>% summarize(tSigma = mean(tSigma, na.rm = T))        
  }
  fig_sim <- ggplot(data = simStats, aes_string(x = 'Svals', y = statType, linetype = 'PRESENT', colour = "TARGET_NAME")) +
    geom_point() +
    facet_wrap(~PRESENT, scale = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig_sim, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'pink_', 'sim_', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig_sim)  
}

template_stat <- function(statType = 'tSigma', bPresent = 0, target = "vertical") {
  
  library(dplyr)
  library(ggplot2)
  
  templateSigmaAbs  <- read.table('~/Dropbox/Calen/Dropbox/tMatch_Response_sigma.txt', header = T, sep = '\t') %>% mutate(PRESENT = 0)
  templateMeanAbs   <- read.table('~/Dropbox/Calen/Dropbox/tMatch_Response_mean.txt', header = T, sep = '\t') %>% mutate(PRESENT = 0)
  
  templateSigma     <- rbind(templateSigmaAbs, templateMeanAbs)
  
  templateAll       <- templateSigma %>% group_by(L,C,S,PYRAMIDLVL,TARGET) %>% arrange(L,C,S,TARGET, PRESENT)
  
  bin_values        <- get_bin_values()
  templateStats     <- templateAll %>% merge(., bin_values) %>% filter(TARGET_NAME == target) %>% filter(PRESENT == bPresent)
  
  templateStats$Cvals <- as.factor(templateStats$Cvals)
  
  fig               <- ggplot(data = templateStats, aes_string(x = 'Svals', y = statType, shape = 'TARGET_NAME', colour = 'Cvals')) +
    geom_point() + 
    facet_wrap(~Lvals, ncol = 3, scales = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig, file = paste0('~/Dropbox/Calen/Dropbox/', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig)
}

edge_stat <- function(statType = 'tSigma', bPresent = 0, target = "vertical") {
  
  edgeSigmaAbs <- read.table('~/Dropbox/Calen/Dropbox/Eabs_Response_sigma.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(0))
  edgeSigmaPres <- read.table('~/Dropbox/Calen/Dropbox/Epres_Response_sigma.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(1))
  
  edgeMeanAbs  <- read.table('~/Dropbox/Calen/Dropbox/Eabs_Response_mean.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(0), tmu = tSigma, tSigma = NULL)
  edgeMeanPres <- read.table('~/Dropbox/Calen/Dropbox/Epres_Response_mean.txt', header = T, sep = '\t') %>% mutate(PRESENT = as.factor(1), tmu = tSigma, tSigma = NULL)
  
  edgeAbs  <- merge(edgeSigmaAbs, edgeMeanAbs)
  edgePres <- merge(edgeSigmaPres, edgeMeanPres)

  edgeAll  <- rbind(edgeAbs, edgePres) %>% group_by(L,C,S,PYRAMIDLVL,TARGET) %>% arrange(L,C,S,TARGET, PRESENT)
  
  bin_values         <- get_bin_values()
  
  
  edgeStatistics     <- edgeAll %>% merge(., bin_values) %>% filter(TARGET_NAME == target) %>% filter(PRESENT == bPresent)
  # All bins
  if(statType == "tmu") {
    fig <- ggplot(data = edgeStatistics, aes(x = Svals, y = tmu, shape = (TARGET_NAME), colour = as.factor(Cvals))) +
      geom_point() + 
      facet_wrap(~Lvals, ncol = 3) +
      theme(aspect.ratio = 1) +
      scale_color_brewer(palette="Spectral")
  }else{
    fig <- ggplot(data = edgeStatistics, aes(x = Svals, y = tSigma, shape = (TARGET_NAME), colour = as.factor(Cvals))) +
      geom_point() + 
      facet_wrap(~Lvals, ncol = 3) +
      theme(aspect.ratio = 1) +
      scale_color_brewer(palette="Spectral")
  }
  ggsave(fig, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")  


  # Only Contrast
  contrastStats <- edgeAll %>% merge(., bin_values)
  if(statType == "tmu") {
    contrastStats <- contrastStats %>% group_by(Cvals, PRESENT, TARGET_NAME) %>% summarize(tmu = mean(tmu, na.rm = T))    
  }else{
    contrastStats <- contrastStats %>% group_by(Cvals, PRESENT, TARGET_NAME) %>% summarize(tSigma = mean(tSigma, na.rm = T))        
  }

  fig_contrast <- ggplot(data = contrastStats, aes_string(x = 'Cvals', y = statType, colour = "TARGET_NAME")) +
    geom_point() +
    facet_wrap(~PRESENT, scales = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig_contrast, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'con_', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
    
  plot(fig_contrast)
  
  # Only Luminance
  
  lumStats <- edgeAll %>% merge(., bin_values)
  
  if(statType == "tmu") {
    lumStats <- lumStats %>% group_by(Lvals, PRESENT, TARGET_NAME) %>% summarize(tmu = mean(tmu, na.rm = T))    
  }else{
    lumStats <- lumStats %>% group_by(Lvals, PRESENT, TARGET_NAME) %>% summarize(tSigma = mean(tSigma, na.rm = T))        
  }
  
  fig_lum <- ggplot(data = lumStats, aes_string(x = 'Lvals', y = statType, colour = "TARGET_NAME")) +
    geom_point() +
    facet_wrap(~PRESENT, scales = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig_lum, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'lum_', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig_lum)

  # Only similarity
  simStats <- edgeAll %>% merge(., bin_values)
  if(statType == "tmu") {
    simStats <- simStats %>% group_by(Svals, PRESENT, TARGET_NAME) %>% summarize(tmu = mean(tmu, na.rm = T))    
  }else{
    simStats <- simStats %>% group_by(Svals, PRESENT, TARGET_NAME) %>% summarize(tSigma = mean(tSigma, na.rm = T))        
  }
  fig_sim <- ggplot(data = simStats, aes_string(x = 'Svals', y = statType, linetype = 'PRESENT', colour = "TARGET_NAME")) +
    geom_point() +
    facet_wrap(~PRESENT, scale = "free_y") +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette="Spectral")
  
  ggsave(fig_sim, file = paste0('~/Dropbox/Calen/Dropbox/edge_stats/', 'sim_', target, '_', statType, '_', as.character(bPresent), '.pdf'), width = 50, height = 50, units = "cm")
  
  plot(fig_sim)  
}



