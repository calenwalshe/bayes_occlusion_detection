source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R")
source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/model_psychometrics.R")


load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/human.psychometrics.rdata")
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.error2d.rdata")
#load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.mahal.rdata")

model.all <- model.error %>% mutate(pc = pnorm(dprime/2)) %>% select(BIN, TARGET, observer, SUBJECT, eccentricity,dprime) %>% filter(SUBJECT %in% c("optimal", "nocov"))
bin.values <- get_experiment_bin_values()


model.all$observer <- model.all$SUBJECT
model.all$SUBJECT <- NULL

model.mahal <- model.mahal %>% select(BIN, TARGET, observer, eccentricity,dprime)

model.all <- rbind(model.all, model.mahal)

model.responses <- model.all %>% filter(observer %in% c("optimal", "nocov", "mahalanobis")) %>% group_by(TARGET, BIN, observer) %>% nest(eccentricity, dprime) %>% as_tibble()

scale.par <- by_row(model.responses, function(x) {
  dat <- data.frame(TARGET = x$TARGET, BIN = x$BIN, observer = x$observer, eccentricity = x$data[[1]]$eccentricity, dprime = x$data[[1]]$dprime)
  human.psychometrics.ave <- human.psychometrics %>% group_by(TARGET, BIN) %>% summarize(threshold = mean(threshold)) %>% mutate(observer = "ave")
  f.optim <- function(scale) {
    psy.tmp <- get_model_psychometric(dat, scale)
    dprime.at.t <- get_dprime_at_eccentricity(psy.tmp, human.psychometrics.ave)    
    
    if( (min(psy.tmp$data[[1]]) < 1) & (max(psy.tmp$data[[1]]) > 1) ) {
      optim.val <- sum((dprime.at.t$dprime_at_threshold - 1)^2)
    } else {
      optim.val <- NA
    }
    return(optim.val)
  }
  optim <- optim(.1, f.optim, method = "Brent", lower = 0, upper = 1)$par
  
  psy.tmp <- get_model_psychometric(dat, optim)
  fitted.dprime <- get_dprime_at_eccentricity(psy.tmp, human.psychometrics.ave)
  
  return <- data.frame(optim = optim, fitted.dprime = fitted.dprime)
})



model.responses$efficiency <- unlist(map(scale.par$.out, 1))
model.responses$fitted.dprime <- unlist(map(scale.par$.out, function(x) {
  x$fitted.dprime.dprime_at_threshold
}))


model.responses <- merge(model.responses, bin.values) %>% as_tibble() %>% arrange(TARGET, BIN)

model.efficiency <- model.responses
save(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/efficiency.rdata", model.efficiency)


## Alternative Method
summarize <- dplyr::summarise
source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R")
source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/model_psychometrics.R")


load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/human.psychometrics.rdata")
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.error2d.rdata")
#load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.mahal.rdata")

model.all <- model.error %>% mutate(pc = pnorm(dprime/2)) %>% 
  select(BIN, TARGET, observer, SUBJECT, eccentricity,dprime)
bin.values <- get_experiment_bin_values()


model.all$observer <- model.all$SUBJECT
model.all$SUBJECT <- NULL

#model.mahal <- model.mahal %>% select(BIN, TARGET, observer, eccentricity,dprime)

#model.all <- rbind(model.all, model.mahal)

model.responses <- model.all %>% 
  group_by(TARGET, BIN, observer) %>% 
  nest(eccentricity, dprime) %>% as_tibble()

human.psychometrics.ave <- human.psychometrics %>% group_by(TARGET, BIN) %>% summarize(threshold = mean(threshold)) %>% mutate(observer = "ave")
model.psy.tmp <- get_model_psychometric(model.all, 1)
dprime.at.t <- get_dprime_at_eccentricity_interpolation(model.error, human.psychometrics.ave) %>% select(-subject_name) %>% mutate(efficiency = 1/dprime_at_threshold)

model.efficiency <- dprime.at.t
save(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/efficiency.rdata", model.efficiency)


