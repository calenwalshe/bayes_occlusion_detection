#######
# Add a model that simplifies the noise to the average of the variance of the two classes.

load('~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/model_wide.rdata')

model.wide.reduced <- model.wide

model.wide.reduced$data_cov_mat_0 <- map(model.wide.reduced$data_cov_mat_0, function(x) diag(x))
model.wide.reduced$data_cov_mat_1 <- map(model.wide.reduced$data_cov_mat_1, function(x) diag(x))

model.wide.reduced$sd_ave <- map2(model.wide.reduced$data_cov_mat_0, model.wide.reduced$data_cov_mat_1, function(x,y) {sqrt((x + y)/2)})

model.wide.reduced$delta_mean <- map2(model.wide.reduced$data_mean_vec_0, model.wide.reduced$data_mean_vec_1, function(x,y) {abs(x-y)})

model.wide.reduced$dprime <- map2(model.wide.reduced$delta_mean, model.wide.reduced$sd_ave, function(x,y) {
  sqrt(sum((x/y)^2))
})

model.wide.reduced <- model.wide.reduced %>% mutate(observer = "model", sub_type = "std_ave", SUBJECT = "model", model.type = "std_ave")


#######

#######
# Add a model that simplifies using only the variance of target absent distribution

# Add a model that simplifies the noise to the average of the variance of the two classes.

model.wide.reduced.abs <- model.wide

model.wide.reduced.abs$data_cov_mat_0 <- map(model.wide.reduced.abs$data_cov_mat_0, function(x) diag(x))
model.wide.reduced.abs$data_cov_mat_1 <- map(model.wide.reduced.abs$data_cov_mat_1, function(x) diag(x))

model.wide.reduced.abs$sd_abs <- map2(model.wide.reduced.abs$data_cov_mat_0, model.wide.reduced.abs$data_cov_mat_1, function(x,y) {sqrt(x)})

model.wide.reduced.abs$delta_mean <- map2(model.wide.reduced.abs$data_mean_vec_0, model.wide.reduced.abs$data_mean_vec_1, function(x,y) {
  abs(x-y)
})

model.wide.reduced.abs$dprime <- map2(model.wide.reduced.abs$delta_mean, model.wide.reduced.abs$sd_abs, function(x,y) {
  sqrt(sum((x/y)^2))
})

model.wide.reduced.abs <- model.wide.reduced.abs %>% mutate(observer = "model", sub_type = "std_abs", SUBJECT = "model", model.type = "std_abs")

intersect.names <- intersect(names(model.wide.reduced.abs), names(optimal.model))

optimal.all <- rbind(optimal.model[, intersect.names], model.wide.reduced[, intersect.names], model.wide.reduced.abs[, intersect.names]) %>% unnest(dprime)

model.psychometrics <- get_model_psychometric(optimal.all)

#######

efficiency.df <- get_dprime_at_eccentricity(model.psychometrics, human.psychometrics)
plot_efficiency(efficiency.df)
