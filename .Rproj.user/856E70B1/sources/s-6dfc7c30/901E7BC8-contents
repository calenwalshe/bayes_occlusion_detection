#' Import the model responses from the raw data.
get_template_response <- function(file_path) {
    if (missing(file_path)) {
        stop('Missing file path to raw data')
    }
  
    library(tidyr)
    library(dplyr)
    
    template_response <- get_raw_data(file_path) %>%
      as_tibble()
    
    template_response <- template_response %>% mutate(TARGET = ifelse(TARGET == 
        1, "vertical", ifelse(TARGET == 2, "horizontal", ifelse(TARGET == 
        3, "bowtie", ifelse(TARGET == 4, "spot", "error"))))) %>%
      filter(PYRAMIDLVL %in% c(1,2,3,4,5))
    
    template_response$TARGET <- as.factor(template_response$TARGET)
    
    bin_values <- get_experiment_bin_values()
    
    template_response <- merge(template_response, bin_values)
    
    pyramidLvl <- data.frame(PYRAMIDLVL = c(1, 2, 3, 4, 5), 
        eccentricity = round(c(0.000000,1.261062,4.172228,12.255842,54.468623),2), 
        downsampleR = c(1, 2, 4, 8, 16))
    
    template_response <- merge(template_response, pyramidLvl) %>% 
        arrange(TARGET, PYRAMIDLVL, L, C, S)
    
    template_response$TPRESENT <- ordered(template_response$TPRESENT, 
        labels = c("absent", "present"))
    
    # Ordering the grouping factor to display nicely.
    template_response$statType <- factor(template_response$statType, 
        levels <- c("Lvals", "Cvals", "Svals"))
    
    template_response <- template_response %>% arrange(SAMPLE, 
        BIN, L, C, S, eccentricity, downsampleR, TARGET, statType, 
        statValue)
    
    template_response <- template_response %>% filter((statType == 
        "Lvals" & BIN %in% c(1, 2, 3, 4, 5, 6)) | (statType == 
        "Cvals" & BIN %in% c(3, 7, 8, 9, 10)) | (statType == 
        "Svals" & BIN %in% c(3, 11, 12, 13, 14, 15)))
    
    template_response$BIN <- as.factor(template_response$BIN)
    
    # Import multiple response function types.
    template_response$RESPONSEFUNC_NAME <- str_extract(template_response$RESPONSEFUNC_NAME, "edge|TR")
    
    template_response <- template_response %>% 
      rename(function_name = RESPONSEFUNC_NAME) %>%
      arrange(TARGET, 
        BIN, statType, eccentricity) %>%
      as_tibble()
    
    template_response$function_name[template_response$function_name == "edge"] <- "edge_cos"
    template_response$function_name[template_response$function_name == "TR"] <- "pattern_only"
    
    return(template_response)
}

#' Import model responses from file
get_raw_data <- function(file_path) {
    if (missing(file_path)) {
        stop('Missing file path')
    }
    model_data <- read.table(file_path, sep = "\t", header = T)
}

#' Return bin values used in the experiment
get_experiment_bin_values <- function() {
    BIN        <- seq(1, 15, 1)
    
    luminance  <- c(1, 3, 5, 7, 9, 10, 5, 5, 5, 5, 5, 5, 5, 5, 5)
    contrast   <- c(5, 5, 5, 5, 5, 5, 3, 7, 9, 10, 5, 5, 5, 5, 5)
    similarity <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 3, 7, 9, 10)
    
    experiment_bins <- data.frame(BIN = BIN, L = luminance, C = contrast, S = similarity)
    
    bin_values <- get_bin_values() %>%
      mutate(TARGET = TARGET_NAME) %>% 
      select(-TARGET_NAME)
    
    
    experiment_bin_values <- merge(experiment_bins, bin_values)
    
    experiment_bin_values <- experiment_bin_values %>% 
      gather(key = statType, value = statValue, Lvals, Cvals, Svals)
    
    experiment_bin_values <- experiment_bin_values %>% filter((statType == 
        "Lvals" & BIN %in% c(1, 2, 3, 4, 5, 6)) | (statType == 
        "Cvals" & BIN %in% c(3, 7, 8, 9, 10)) | (statType == 
        "Svals" & BIN %in% c(3, 11, 12, 13, 14, 15))) %>% arrange(TARGET, 
        statType, statValue)
    
    # Ordering the grouping factor to display nicely.
    experiment_bin_values$statType <- factor(experiment_bin_values$statType, 
        levels <- c("Lvals", "Cvals", "Svals"))
    
    
    experiment_bin_values$BIN <- as.factor(experiment_bin_values$BIN)
    return(experiment_bin_values)
}

#' Return a dataframe that contains all the statistics of the bins that are used.
get_bin_values <- function() {
    
    Lvals <- (c(2.9627, 4.1447, 5.7984, 8.1117, 11.348, 15.8755, 22.2093, 31.07, 43.4659, 60.8073))
    Cvals <- (c(0.0247, 0.0364, 0.0536, 0.079, 0.1163, 0.1713, 0.2523, 0.3716, 0.5472, 0.8059))
    Sa1vals <- (c(0.4255, 0.4577, 0.4923, 0.5295, 0.5696, 0.6127, 0.659, 0.7088, 0.7625, 0.8201))
    Sa2vals <- (c(0.4454, 0.4774, 0.5117, 0.5484, 0.5878, 0.6299, 0.6751, 0.7236, 0.7755, 0.8312))
    Sa3vals <- (c(0.5394, 0.5584, 0.5781, 0.5986, 0.6197, 0.6416, 0.6643, 0.6877, 0.712, 0.7372))
    Sa4vals <- (c(0.6973, 0.712, 0.727, 0.7423, 0.7579, 0.7739, 0.7902, 0.8069, 0.8239, 0.8412))
    
    L   <- data.frame(Lvals, L = 1:10)
    C   <- data.frame(Cvals, C = 1:10)
    Sa1 <- data.frame(Svals = Sa1vals, S = 1:10, TARGET = factor(rep("1", 10)), TARGET_NAME = factor(rep("vertical", 10)))
    Sa2 <- data.frame(Svals = Sa2vals, S = 1:10, TARGET = factor(rep("2", 10)), TARGET_NAME = factor(rep("horizontal", 10)))
    Sa3 <- data.frame(Svals = Sa3vals, S = 1:10, TARGET = factor(rep("3", 10)), TARGET_NAME = factor(rep("bowtie", 10)))
    Sa4 <- data.frame(Svals = Sa4vals, S = 1:10, TARGET = factor(rep("4", 10)), TARGET_NAME = factor(rep("spot", 10)))
    S   <- rbind(Sa1, Sa2, Sa3, Sa4)
    
    bin_values <- merge(merge(L, C), S)
}

#' Get eccentricities
#' Get receptive field spacing based on Drasdo (2007).
get_eccentricity <- function(model_space, direction = "temporal") {
  library(Hmisc)
  pyramid_space    <- 1/(120/(2^seq(0,4,1)))
  degrees          <- c(0, .25, .5, 1, 2, 5, 10, 20, 25, 30)
  gc_count         <- c(27930, 19438, 14614, 9365, 4986, 1633, 578, 231, 214, 182)
  spacing          <- 1/sqrt(gc_count/2)
  spacing[1]          <- 1/120
  
  # michaelis menten
  
  f <- function(xx) {
    s_max <- xx[1]
    a     <- xx[2]
    b     <- 1/120
    
    pred.f <- function(degrees) {
      r.val <- s_max * degrees / (a + degrees) + b
    }
    
    sqrt(sum((spacing - pred.f(degrees))^2))
  }
  
  f.1 <- function(s_max, a, b = 1/120, x) {
    r.val <- s_max * x / (a + x) + b
  }
  
  start.val <- c(1,1)
  params <- optim(start.val, f)
  
  x.1 <- params$par[1]
  x.2 <- params$par[2]
  #x.3 <- params$par[3]
  
  output <- f.1(s_max = x.1, a = x.2, x = seq(0, 60, .01))
  
  my.dat <- data.frame(degrees = degrees, data = spacing) %>% mutate(type = "obs")
  my.dat.1 <- data.frame(degrees = seq(0, 60, .01), data = output) %>% mutate(type = "pred")
  
  ggplot(my.dat, aes(x = degrees, y = data, colour = type)) + geom_point() + geom_line(data = my.dat.1, aes(x = degrees, y = data), inherit.aes = TRUE) + coord_cartesian(ylim = c(0, .14)) + geom_hline(yintercept = 1/c(120 / 2^c(0, 1, 2, 3, 4)))
  
  s_max <- x.1
  a <- x.2
  b <- 1/120
  
  c <- (pyramid_space - 1/120)/s_max
  deg.interp <- c*a / (1 - c)

  return(deg.interp)
}
