load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/human.psychometrics.rdata")
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/optimal.model.all.rdata")
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/model.psychometrics.unscaled.rdata")
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/optim.scale.rdata")

source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R")

optimal.all.1 <- optimal.all %>% filter(sub_type == "optimal") %>% select(-SUBJECT)

efficiency <- get_dprime_at_eccentricity(model.psychometrics.unscaled, human.psychometrics)

plot_efficiency(efficiency)
