interpolated.model <- function(template_response, human.responses) {
# Interpolate the model response for measured human eccentricities.
  human.eccentricities <- human.responses %>%
    select(TRIAL, TARGET, BIN, SUBJECT, SESSION, ECCENTRICITY, PATCHID) %>%
    unique() %>%
    rename(eccentricity = ECCENTRICITY)
  
}

#' Title
#'
#' @param conditions
#'
#' @return
#' @export
#'
#' @import broom dplyr
#' @examples
get_all_interpolated <- function(template_response) {
    library(broom)
    df <- expand.grid(c(as.character(unique(template_response$TARGET))), 
        c(seq(0.5, 24, length.out = 16)), in_bin = 1:15, in_function_name = c("edge_cos", 
            "mean_only", "pattern_only"), in_tpresent = c("present", 
            "absent"))
    
    df_len <- nrow(df)
    
    interpolated_responses <- mclapply(1:df_len, FUN = function(x) get_interpolated_model_responses(template_response, 
        as.character(df[x, ][[1]]), df[x, ][[2]], df[x, ][[3]], 
        as.character(df[x, ][[4]]), as.character(df[x, ][[5]])), 
        mc.cores = 16)
    
    return(interpolated_responses)
}

interpolate_dprime <- function(model.dprime, human.dprime) {
  
  model.dprime <- model.dprime %>% rename(type = SUBJECT)
  
  human.dprime$BIN <- as.factor(human.dprime$BIN)
  model.dprime$BIN <- as.factor(model.dprime$BIN)
  
  human.dprime$TARGET <- as.factor(human.dprime$TARGET)
  model.dprime$TARGET <- as.factor(model.dprime$TARGET)
  
  human.dprime <- human.dprime %>%
    select(SUBJECT, TARGET, BIN, eccentricity, dprime) %>% 
    distinct()
  
  human.measure <- human.dprime %>% select(TARGET, BIN, SUBJECT, eccentricity)
  
  dprime.interpolate <- model.dprime %>% group_by(type, TARGET, BIN) %>%
    nest() %>%
    left_join(., human.measure, by = c("TARGET", "BIN")) %>%
    mutate(dprime_hat = map2(eccentricity, data, function(x,y) {
      x_in <- y$eccentricity
      y_in <- y$dprime
      
      approxExtrap(x_in, y_in, x)$y
      
    })) %>%
    unnest(dprime_hat)
  
  dprime.interpolate %>% select(-data) %>% left_join(., human.dprime, by = c("TARGET", "BIN", "SUBJECT", "eccentricity"))
}
