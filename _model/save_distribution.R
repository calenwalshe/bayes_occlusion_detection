#' Compute mean vector and covariance matrices for all conditions in the experiment.
save_gaussian_distributions <- function(template_response) {
    template_response[, c("statType", "statValue", "L", "C", 
        "S")] <- NULL
    
    f <- function(x) {
        model_distribution <- template_response %>% unique() %>% 
            filter(eccentricity == x) %>% dplyr::select(eccentricity, 
            BIN, function_name, SAMPLE, TARGET, TPRESENT, TRESP) %>% 
            mutate(TARGET = factor(TARGET, levels = c("vertical", 
                "horizontal", "bowtie", "spot"))) %>% arrange(SAMPLE, 
            TARGET, BIN, eccentricity, TPRESENT, function_name) %>% 
            group_by(BIN, TARGET, eccentricity, function_name, 
                TPRESENT) %>% spread(function_name, TRESP) %>% 
            summarize(edge_mean = mean(edge_cos), lum_mean = mean(mean_only), 
                pattern_mean = mean(pattern_only), edge_sd = sd(edge_cos), 
                lum_sd = sd(mean_only), pattern_sd = sd(pattern_only), 
                edge_lum_cor = cor(edge_cos, mean_only), edge_pattern_cor = cor(edge_cos, 
                  pattern_only), lum_pattern_cor = cor(mean_only, 
                  pattern_only)) %>% rowwise() %>% mutate(cor_mat = list(matrix(c(1, 
            edge_lum_cor, edge_pattern_cor, edge_lum_cor, 1, 
            lum_pattern_cor, edge_pattern_cor, lum_pattern_cor, 
            1), nrow = 3, ncol = 3)), mean_mat = list(c(edge_mean, 
            lum_mean, pattern_mean)), sd_mat = list(c(edge_sd, 
            lum_sd, pattern_sd))) %>% arrange(BIN, TARGET, TPRESENT, 
            eccentricity) %>% mutate(TARGET = factor(TARGET)) %>% 
            dplyr::select(BIN, TARGET, eccentricity, TPRESENT, 
                cor_mat, sd_mat, mean_mat)
    }
    
    g <- function(x) {
        model_responses <- template_response %>% unique() %>% 
            filter(eccentricity == x) %>% dplyr::select(eccentricity, 
            BIN, function_name, SAMPLE, TARGET, TPRESENT, TRESP) %>% 
            mutate(TARGET = factor(TARGET, levels = c("vertical", 
                "horizontal", "bowtie", "spot"))) %>% mutate(TPRESENT = ifelse(TPRESENT == 
            "absent", 0, 1)) %>% dplyr::select(SAMPLE, BIN, TPRESENT, 
            TARGET, function_name, eccentricity, TRESP) %>% unique() %>% 
            spread(function_name, TRESP) %>% rowwise() %>% mutate(response_vec = list(c(edge_cos, 
            mean_only, pattern_only))) %>% mutate(TARGET = factor(TARGET)) %>% 
            dplyr::select(SAMPLE, BIN, TPRESENT, TARGET, response_vec)
    }
    
    mclapply(unique(template_response$eccentricity), FUN = function(x) {
        model_distribution <- f(x)
        save(file = paste0("~/Dropbox/Calen/Work/data_storage/eccentricity_project/model_distribution_", 
            as.character(floor(x)), ".rdata"), model_distribution)
    }, mc.cores = 16)
    
    mclapply(unique(template_response$eccentricity), FUN = function(x) {
        model_responses <- g(x)
        save(file = paste0("~/Dropbox/Calen/Work/data_storage/eccentricity_project/model_responses_", 
            as.character(floor(x)), ".rdata"), model_responses)
    }, mc.cores = 16)
}
