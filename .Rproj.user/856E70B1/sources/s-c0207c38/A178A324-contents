library(dplyr)
library(ggplot2)
library(purrr)
library(purrrlyr)
library(tidyr)
summarize <- dplyr::summarize

summary.raw.data <- raw.data %>%
  group_by(SUBJECT, BIN, TARGET, ECCENTRICITY) %>%
  summarize(pc = sum(CORRECT) / n())

lapply(c("vertical", "horizontal", "bowtie", "spot"), FUN = function(x) {
  
  fig <- ggplot(summary.raw.data %>% filter(TARGET == x), aes(x = ECCENTRICITY, y = pc, colour = SUBJECT)) +
    facet_wrap(~BIN, ncol = 4, scales = "free_x") + 
    geom_point() +
    geom_line() + 
    theme(aspect.ratio = 1) +
    coord_cartesian(ylim = c(.5, 1))
  
  ggsave(filename = paste0("~/Dropbox/Calen/Dropbox/psychometrics.raw-", x, '.pdf'), device = "pdf", plot = fig)
})
