get_interpolated_model_responses <- function(template_response, 
    in_target, in_eccentricity, in_bin, in_function_name, in_tpresent) {
    
    # print(in_target) print(in_eccentricity) print(in_bin)
    # print(in_function_name) print(in_tpresent)
    
    pyr_ecc <- unique(template_response$eccentricity)
    pyr_ecc_low <- c(max(which((in_eccentricity - pyr_ecc) > 
        0)), min(which((in_eccentricity - pyr_ecc) < 0)))
    
    ecc_vals <- pyr_ecc[pyr_ecc_low]
    
    d.1 <- template_response %>% filter(eccentricity %in% ecc_vals, 
        TARGET == in_target, BIN == in_bin, function_name == 
            in_function_name, TPRESENT == in_tpresent) %>% select(-statType, 
        -statValue, -L, -C, -S) %>% distinct() %>% data.frame()
    
    lm.1 <- d.1 %>% mutate(TRESP = TRESP) %>% group_by(SAMPLE) %>% 
        do(lin_mod = augment(lm(TRESP ~ eccentricity, data = .), 
            newdata = data.frame(eccentricity = in_eccentricity))) %>% 
        select(lin_mod)
    
    new.dat <- unlist(lapply(lm.1$lin_mod, FUN = function(x) x$.fitted))
    
    
    d.2 <- d.1 %>% filter(eccentricity == ecc_vals[1])
    d.2$eccentricity <- in_eccentricity
    d.2$TRESP <- new.dat
    
    
    return(d.2)
    
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
    source(broom)
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
