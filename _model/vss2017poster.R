plot_vss_psychometrics <- function(human.psychometrics, performance_measures) {
    sub_params <- psychometric_parameters %>% filter(BIN == 3) %>% data.frame()
    
    sub_response <- performance_measures %>% filter(BIN == 3) %>% data.frame()
    
    responses <- lapply(1:nrow(sub_params), FUN = function(x) data.frame(TARGET = sub_params[x, "TARGET"], SUBJECT = sub_params[x, "SUBJECT"], eccentricity = seq(0, 30, 0.1), pc = pnorm(sub_params[x, "d0"]/2 * (sub_params[x, "e0"]^sub_params[x, "b"])/(sub_params[x, 
        "e0"]^sub_params[x, "b"] + seq(0, 30, 0.1)^sub_params[x, "b"])))) %>% do.call(rbind, .)
    
    responses <- responses %>% mutate(BIN = 3)
    
    responses$TARGET <- factor(responses$TARGET, levels(responses$TARGET)[c(4, 2, 1, 3)])
    sub_response$TARGET <- factor(sub_response$TARGET, levels(sub_response$TARGET)[c(4, 2, 1, 3)])
    
    p <- ggplot(responses, aes(x = eccentricity, y = pc, colour = SUBJECT)) + geom_line() + geom_point(data = sub_response, aes(x = eccentricity, y = percent_correct), size = 2) + scale_y_continuous(limits = c(0.4, 1)) + facet_wrap(~TARGET, ncol = 2) + theme_light() + 
        theme(legend.title = element_blank(), strip.text.x = element_text(size = 20), aspect.ratio = 1, axis.title = element_text(size = 20), axis.text = element_text(size = 17), legend.text = element_text(size = 20)) + ylab("Percent Correct") + xlab("Eccentricity")
    
    ggsave(last_plot(), file = paste0("~/Dropbox/Calen/Dropbox/tmp_images/vss_psychometrics.pdf"))
    
    
}
