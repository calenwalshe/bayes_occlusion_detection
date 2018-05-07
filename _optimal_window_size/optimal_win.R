library(dplyr)
library(stringr)
library(ggplot2)

path.dat <- '~/Dropbox/Calen/Work/occluding/detection_model_analysis/_optimal_window_size/_data'
files    <- files <- list.files(path = path.dat, full.names = T)
template.response.list <- lapply(files, FUN = function(x) get_template_response(x) %>% mutate(window_code = x))
template.response      <- do.call(rbind, template.response.list) %>% 
  mutate(window_code = as.factor(as.numeric(str_extract(window_code, "(?<=data/)[:digit:]+"))))

levels(template.response$window_code) <- c(0, 1, 2, .5, 1.5)

template.response$window_code <- factor(template.response$window_code, levels(template.response$window_code)[order(levels(template.response$window_code))])

template.response.window_code <- template.response %>% 
  filter(function_name == "pattern_only") %>%
  group_by(window_code, PYRAMIDLVL, BIN, TPRESENT) %>% 
  summarize(tMean = mean(TRESP), tVar = var(TRESP)) %>%
  summarize(dprime = diff(tMean)/(sqrt(mean(tVar)))) %>%
  arrange(window_code)

template.response.window_code %>%
  group_by(PYRAMIDLVL, BIN) %>%
  filter(dprime == max(dprime)) %>%
  arrange(PYRAMIDLVL) %>%
  data.frame()

ggplot(template.response.window_code, aes(x = PYRAMIDLVL, y = (dprime), colour = window_code)) + 
  geom_jitter(height = 0)
  

