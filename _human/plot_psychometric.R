#' Psychometric functions contained in the supplied dataframe will be visualized with the human data overlaid.
#'
#' @return
#' @export
#'
#' @import dplyr ggplot2
#'
#' @examples
plot_human_psychometrics <- function(params, response_data, out_path, 
    target = "vertical") {
    if (missing(params)) {
        params <- load_cached_parameters()
    }
    if (missing(out_path)) {
        out_path <- tempfile(pattern = "file", tmpdir = tempdir(), 
            fileext = ".pdf")
    }
    if (missing(response_data)) {
        stop("Missing data")
    }
    if (!any(names(params) == "gamma_bias")) {
        params$gamma_bias <- 0
    }
    
    params <- params %>% filter(TARGET %in% target)
    response_data <- response_data %>% filter(TARGET %in% target)
    
    interpolated_responses <- generate_epf_observations(params, 
        gamma_bias = TRUE)
    
    interpolated_responses <- interpolated_responses %>% filter(psychometric_type %in% 
        c("percent_correct", "no_bias")) %>% mutate(has_bias = ifelse(gamma_bias == 
        0, FALSE, TRUE)) %>% mutate(percent_correct = PC)
    
    interpolated_responses <- as.data.frame(lapply(interpolated_responses, 
        function(x) if (is.factor(x)) 
            factor(x) else x)) %>% filter(has_bias == FALSE)
    response_data <- as.data.frame(lapply(response_data, function(x) if (is.factor(x)) 
        factor(x) else x))
    
    interpolated_responses$BIN <- factor((interpolated_responses$BIN), 
        levels = seq(1, 15, 1))
    
    psychometrics_percentcorrect <- ggplot(interpolated_responses, 
        aes(x = eccentricity, y = percent_correct)) + geom_line() + 
        geom_point(data = response_data, aes(x = eccentricity, 
            y = percent_correct)) + facet_wrap(statType ~ statValue + 
        BIN, scales = "free") + scale_x_continuous(limits = c(0, 
        30)) + scale_y_continuous(limits = c(0.25, 1)) + theme(axis.line = element_line()) + 
        ggtitle(target)
    
    
    ggsave("~/Dropbox/Calen/tmp_fig/tmp_im.pdf", height = 40, 
        width = 40, units = "in", device = "pdf")
}

#' Returns interpolated observations for the eccentricity psychometric functions stored in params
#'
#' @param params
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples

generate_epf_observations <- function(params, gamma_bias = FALSE) {
    if (missing(params)) {
        error("Parameters missing")
    }
    
    x <- seq(0, 30, 0.01)
    interpolated_data <- data.frame(BIN = NULL, TARGET = NULL, 
        SUBJECT = NULL, eccentricity = NULL, PC = NULL, gamma_bias = NULL, 
        psychometric_type = NULL)
    
    for (i in unique(params$BIN)) {
        for (j in unique(params$TARGET)) {
            for (k in unique(params$SUBJECT)) {
                sub_params <- params %>% filter(BIN == i, TARGET == 
                  j, SUBJECT == k) %>% .[1, ]
                
                d0 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                  eccentricity = x, PC = with(sub_params, pnorm(0.5 * 
                    plo * (e0^b)/(e0^b + x^b), mean = 0, sd = 1)), 
                  gamma_bias = 0, psychometric_type = "no_bias")
                d1 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                  eccentricity = x, PC = with(sub_params, pnorm(0.5 * 
                    dmax * (e0^b)/(e0^b + x^b) - gamma_bias, 
                    mean = 0, sd = 1)), gamma_bias = sub_params$gamma_bias, 
                  psychometric_type = "hit")
                d2 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                  eccentricity = x, PC = with(sub_params, pnorm(-0.5 * 
                    dmax * (e0^b)/(e0^b + x^b) - gamma_bias, 
                    mean = 0, sd = 1)), gamma_bias = sub_params$gamma_bias, 
                  psychometric_type = "falsealarm")
                d3 <- data.frame(BIN = i, TARGET = j, SUBJECT = k, 
                  eccentricity = x, PC = with(sub_params, 0.5 * 
                    pnorm(0.5 * dmax * (e0^b)/(e0^b + x^b) - 
                      sub_params$gamma_bias) + 0.5 * (1 - pnorm(-0.5 * 
                    dmax * (e0^b)/(e0^b + x^b) - sub_params$gamma_bias))), 
                  gamma_bias = sub_params$gamma_bias, psychometric_type = "percent_correct")
                
                interpolated_data <- rbind(interpolated_data, 
                  d0, d1, d2, d3)
            }
        }
        
    }
    
    experiment_bin_values <- get_experiment_bin_values()
    interpolated_data <- merge(interpolated_data, experiment_bin_values, 
        by = c("BIN", "TARGET"))
    return(interpolated_data)
}


