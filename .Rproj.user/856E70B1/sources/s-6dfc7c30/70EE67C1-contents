rm(list = ls())

library(R.matlab)

load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model_wide.rdata")

dprime.vals.matlab <- read.csv("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/matlab_export/classify_output.txt", sep = '\t', header = T) %>% 
  as_tibble()

dprime.vals.matlab <- dprime.vals.matlab %>% rename(eccentricity.round = ECCENTRICITY)
dprime.vals.matlab$BIN <- factor(dprime.vals.matlab$BIN)

model.wide$eccentricity.round <- round(model.wide$eccentricity)

model.wide <- left_join(dprime.vals.matlab, model.wide, by = c("BIN", "TARGET", "eccentricity.round"))

model.wide <- model.wide %>% rename(dprime = DPRIME) %>% select(-eccentricity.round)

model.error <- model.wide

model.error$observer <- "optimal"
model.error$SUBJECT <- "optimal"

save(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.error2d.rdata", model.error)
