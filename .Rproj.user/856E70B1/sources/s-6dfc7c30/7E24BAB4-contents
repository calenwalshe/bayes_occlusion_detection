rm(list = ls())

library(R.matlab)

# Export Data to matlab
# *crazy because I already imported it from matlab, but here we go. 


source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/import_model.R")
source("~/Dropbox/Calen/Work/occluding/occlusion_detect/_model/model_psychometrics.R")

load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model_wide.rdata")

by_row(model.wide, function(x) {
  dat.absent  <- as.matrix(x$data_response_vec_0[[1]])
  dat.present <- as.matrix(x$data_response_vec_1[[1]])
  
  file.path <- '~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/matlab_export/'
  
  file.name.absent  <- paste0(file.path, paste(x$BIN, as.character(x$TARGET), round(as.numeric(x$eccentricity)), 'absent', sep = '_'), '.mat')
  file.name.present <- paste0(file.path, paste(x$BIN, as.character(x$TARGET), round(as.numeric(x$eccentricity)), 'present', sep = '_'), '.mat')

  writeMat(file.name.absent, dat.absent  = dat.absent)  
  writeMat(file.name.present, dat.present = dat.present)  
})

