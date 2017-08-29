plot_all_single_eccentricity <- function(human.psychometrics, 
    model.dprime) {
    
    sub_params <- human.psychometrics
    
    sub_params <- merge(sub_params, data.frame(eccentricity = unique(model.dprime$eccentricity)))
    
    pc <- lapply(1:nrow(sub_params), FUN = function(x) pc = pnorm(sub_params[x, 
        "d0"]/2 * (sub_params[x, "e0"]^sub_params[x, "b"])/(sub_params[x, 
        "e0"]^sub_params[x, "b"] + sub_params[x, "eccentricity"]^sub_params[x, 
        "b"]))) %>% do.call(rbind, .)
    
    sub_params$pc <- pc
    
    human_dat <- sub_params %>% group_by(TARGET, BIN, eccentricity) %>% 
        summarize(pc = mean(pc)) %>% ungroup() %>% mutate(dprime = 2 * 
        as.numeric(qnorm(pc)), SUBJECT = "human") %>% data.frame
    
    model_dat <- performance_measures %>% select(TARGET, BIN, 
        eccentricity, percent_correct, dprime) %>% rowwise() %>% 
        mutate(pc = percent_correct, dprime = 2 * as.numeric(qnorm(pc)), 
            percent_correct = NULL, SUBJECT = "model") %>% data.frame
    
    dat <- rbind(human_dat, model_dat) %>% group_by(BIN, TARGET, 
        eccentricity) %>% summarize(pc_diff = pc[1] - pc[2])
    
    fig <- ggplot(dat, aes(x = eccentricity, y = pc_diff)) + 
        geom_point() + facet_grid(~TARGET)
    
    ggsave(last_plot(), file = "~/Dropbox/Calen/Dropbox/tmp_images/compare_pc.pdf")
    print(1)
}


