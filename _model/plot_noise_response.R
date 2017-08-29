plot_noise_response <- function(noise.response) {
    out_path = "~/Dropbox/Calen/Dropbox/tmp_images/"
    plot.1 <- ggplot(noise.response, aes(x = noise, y = percent_correct, 
        colour = as.factor(noise))) + geom_point() + facet_grid(TARGET ~ 
        BIN)
    
    ggsave(last_plot(), file = paste0(out_path, "noise_response_fig", 
        ".pdf"), width = 50, height = 50/(1920/1080), units = c("cm"))
    
    
}
