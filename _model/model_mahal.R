# Mahalobis Distance Predicted From the model.

load('~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model_error.rdata')

model.mahal <- model.error %>% filter(sub_type == "optimal")

model.mahal$cov_ave  <- map2(model.mahal$data_cov_mat_0,model.mahal$data_cov_mat_1, function(x,y) (x + y)/2)
model.mahal$cov_ave  <- map2(model.mahal$data_cov_mat_0,model.mahal$data_cov_mat_1, function(x,y) (x + y)/2)
model.mahal          <- model.mahal %>% rowwise() %>% mutate(dprime = sqrt(mahalanobis(data_mean_vec_1, data_mean_vec_0, cov_ave)))
model.mahal$observer <- "mahalanobis"
model.mahal$SUBJECT <- "mahalanobis"
model.mahal$pc       <- pnorm(model.mahal$dprime/2)

save(file = '~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model.mahal.rdata', model.mahal)

