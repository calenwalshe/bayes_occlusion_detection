#### 
# Human psychometrics are fit.
# Eccentricity at a dprime of 1 is computed.
# Dprime for the ideal observer classes is computed for that eccentricity.
# Displayed is the ratio of dprimes at the eccentricity for the human and observers.
####

# Comparison of human and model efficiency
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/model.psychometrics.optimal.all")
load(file = "~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/human.psychometrics.rdata")

human.threshold     <- get_threshold(human.psychometrics)

dprime.efficiency  <- get_dprime_at_eccentricity(model.psychometrics.optimal.all, human.psychometrics)

plot_efficiency(dprime.efficiency)
