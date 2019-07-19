# Generate Bivariate Distributions for Middle Bin

library(lattice)
library(hextri)

load('~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/_model_response/model_wide.rdata')

raw.data <- model.wide %>%
  filter(TARGET == "vertical", BIN == 3) %>%
  select(BIN, TARGET, eccentricity, data_response_vec_0, data_response_vec_1)

all.responses <- lapply(1:5, function(x) {

  absent.e.l <- as_data_frame(raw.data$data_response_vec_0[[x]])[, c("edge", "mean")]
  absent.e.l$response_1_label <- names(absent.e.l)[1]
  absent.e.l$response_2_label <- names(absent.e.l)[2]
  
  absent.e.p <- as_data_frame(raw.data$data_response_vec_0[[x]])[, c("edge", "pattern")]
  absent.e.p$response_1_label <- names(absent.e.p)[1]
  absent.e.p$response_2_label <- names(absent.e.p)[2]
  
  absent.p.l <- as_data_frame(raw.data$data_response_vec_0[[x]])[, c("mean", "pattern")]
  absent.p.l$response_1_label <- names(absent.p.l)[1]
  absent.p.l$response_2_label <- names(absent.p.l)[2]
  
  present.e.l <- as_data_frame(raw.data$data_response_vec_1[[x]])[, c("edge", "mean")]
  present.e.l$response_1_label <- names(present.e.l)[1]
  present.e.l$response_2_label <- names(present.e.l)[2]
  
  present.e.p <- as_data_frame(raw.data$data_response_vec_1[[x]])[, c("edge", "pattern")]
  present.e.p$response_1_label <- names(present.e.p)[1]
  present.e.p$response_2_label <- names(present.e.p)[2]
  
  present.p.l <- as_data_frame(raw.data$data_response_vec_1[[x]])[, c("mean", "pattern")]
  present.p.l$response_1_label <- names(present.p.l)[1]
  present.p.l$response_2_label <- names(present.p.l)[2]
  
  responses.frame.present <- list(present.e.l, present.e.p, present.p.l)
  responses.frame.absent <- list(absent.e.l, absent.e.p, absent.p.l)

  responses.frame.present <- map(responses.frame.present, function(x) {
    names(x)[c(1,2)] <- c("response_1", "response_2");
    return(x)
  }
  )

  responses.frame.absent <- map(responses.frame.absent, function(x) {
    names(x)[c(1,2)] <- c("response_1", "response_2");
    return(x)
  }
  )  
  
  response.frame.present <- do.call(rbind, responses.frame.present)
  response.frame.absent <- do.call(rbind, responses.frame.absent)
  
  response.frame.present$present <- "present"
  response.frame.absent$present  <- "absent"

  response.frame.present$TARGET       <- "vertical"
  response.frame.present$BIN          <- 3
  response.frame.present$eccentricity  <- raw.data$eccentricity[x]
  
  response.frame.absent$TARGET       <- "vertical"
  response.frame.absent$BIN          <- 3
  response.frame.absent$eccentricity  <- raw.data$eccentricity[x]
  
  response.frame <- rbind(response.frame.present, response.frame.absent)
  
  return(response.frame)
})

plot.responses <- do.call(rbind, all.responses)

plot.responses$bivariate.factor <- factor(paste0(plot.responses$response_1_label, '.', plot.responses$response_2_label))

plot.frames <- plot.responses %>% group_by(eccentricity, bivariate.factor,present) %>% nest()

plot.frames$eccentricity <- factor(plot.frames$eccentricity, levels = sort(unique(plot.frames$eccentricity)), labels = round(sort(unique(plot.frames$eccentricity)), 2))
plot.frames$present <- factor(plot.frames$present, levels = c("absent", "present"))

plot.frames$fig <- map(plot.frames$data, function(data){
  names(data)[c(1,2)] <- c(data$response_1_label[1], data$response_2_label[1])

  var1 <- first(data$response_1_label)
  var2 <- first(data$response_2_label)
  
  fig <- ggplot(data, aes_string(x = var1, y = var2)) + 
    geom_hex() +
    theme(aspect.ratio = 1) +
    expand_limits(fill = c(0)) +
    guides(fill = FALSE)
})

edge.mean <- plot.frames %>% filter(bivariate.factor == "edge.mean") %>% arrange(eccentricity, present)


gs <- lapply(edge.mean$fig, function(x) arrangeGrob(x))

grobframe.1 <- arrangeGrob(grobs = gs,nncol = 2, nrow = 5,
                         main = textGrob("Plots", gp = gpar(fontsize=12, fontface="bold.italic", fontsize=12)))

plot(grobframe.1)

