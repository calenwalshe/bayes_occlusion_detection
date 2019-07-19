# Merge human responses with model responses
#
#
# Returns:
#   Save dataframe to disk.
human.dvr <- function(model_wide) {
library(Hmisc)
library(dplyr)
library(mvtnorm)
summarize <- dplyr::summarize

#Format precomputed covariances
model.lik <- model_wide %>% rowwise() %>%
  mutate(lik_1 = list(
           mvtnorm::dmvnorm(data_response_vec_1, data_mean_vec_1, data_cov_mat_1, log = T) - 
                 mvtnorm::dmvnorm(data_response_vec_1, data_mean_vec_0, data_cov_mat_0, log = T)
           ),
         lik_0 = list(
           mvtnorm::dmvnorm(data_response_vec_0, data_mean_vec_1, data_cov_mat_1, log = T) - 
                 mvtnorm::dmvnorm(data_response_vec_0, data_mean_vec_0, data_cov_mat_0, log = T))
         ) %>%
  select(BIN, TARGET, eccentricity, lik_1, lik_0, data_PATCHID_0, data_PATCHID_1)

model.lik.present <- model.lik %>% 
  select(-matches('_0')) %>%
  unnest() %>%
  rename(PATCHID = data_PATCHID_1, TRESP = lik_1) %>%
  mutate(TPRESENT = "present")
  
model.lik.absent <- model.lik %>% 
  select(-matches('_1')) %>%
  unnest() %>%
  rename(PATCHID = data_PATCHID_0, TRESP = lik_0) %>%
  mutate(TPRESENT = "absent")

dvr.model <- rbind(model.lik.present, model.lik.absent) 

# Format model
dvr.model <-
  dvr.model %>% rename(eccentricity_model = eccentricity) %>%
  distinct()

# Import human responses
human.responses <- export.responses()

# Format human responses
human.responses <- human.responses %>%
  filter(TRIAL != 1) %>%
  rename(eccentricity_human = ECCENTRICITY)

human.responses$TPRESENT <-
  ifelse(human.responses$TPRESENT == 1, "present", "absent")


# Merge dataframes
human.response.template <-
  inner_join(human.responses,
            dvr.model,
            by = c("BIN", "PATCHID", "TARGET", "TPRESENT"))

# Interpolate model responses at human eccentricities
interpolated.tresp <- human.response.template %>% 
  group_by(TRIAL, TPRESENT, SUBJECT, BIN, TARGET, PATCHID, eccentricity_human) %>%
  nest() %>%
  mutate(TRESP_interpolated = map2(eccentricity_human,data, function(x,y) {
                TRESP_interpolated = approxExtrap(y$eccentricity_model, y$TRESP, xout = x)$y
  }
  )) %>%
  unnest(TRESP_interpolated, .preserve = data)
  
# Merge human and model (interpolated values)
human.dvr <- inner_join(human.responses, interpolated.tresp, by = c("TRIAL", "BIN", "TPRESENT", "eccentricity_human", "PATCHID", "SUBJECT", "TARGET")) %>%
  arrange(SUBJECT, TARGET, BIN, LEVEL, BIN, TARGET, TRIAL) %>%
  select(-data)

human.dvr <- human.dvr %>% mutate(response = ifelse(HIT | FALSEALARM, 1, 0))

save(file = '~/Dropbox/Calen/Work/occluding/detection_model_analysis/_data/human.dvr.rdata', human.dvr)

return(human.dvr)
}

