# Plot all psychometrics from the scaled eccentricity functions

load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.psychometrics.scaled.rdata")

source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/model_visualization.R")

path <- "~/Dropbox/Calen/Work/occluding/occlusion_detect/_paper1figures/eccentricity_psychometric/_figures/psychometric.pdf"

plot.model.psychometric(model.psychometrics.scaled, path)
