library(bbmle)
library(ggplot2)
library(knitr)
library(kableExtra)


statistics <- read.table('~/Dropbox/Calen/Work/search/scene_statistics/_data/tMatchSigmaNoise.txt', header = T, sep = "\t")

source('~/Dropbox/Calen/Work/occluding/detection_model_analysis/_human/')

LVals <- data.frame(L = 1:10, LVals = c(45.7815,   46.6627,   47.5440,   48.4252,   49.3064,   50.1877,   51.0689,   51.9501,   52.8314,   53.7126))
CVals <- data.frame(C = 1:10, CVals = c(0.0336,    0.0419 ,   0.0502,    0.0586 ,   0.0669,    0.0752 ,   0.0836,    0.0919 ,   0.1002,    0.1086))
SVals <- data.frame(S = 1:10, SVals = c(0.5693,    0.5966 ,   0.6240,    0.6514 ,   0.6788,    0.7062 ,   0.7335,    0.7609 ,   0.7883,    0.8157))

stats.all <- statistics %>% 
  select(-eccentricity) %>% 
  rename(response = tSigma) %>%
  as_tibble() %>%
  mutate(N.norm = N / max(N))

stats.all <- full_join(stats.all, LVals, by = c("L"))
stats.all <- full_join(stats.all, CVals, by = c("C"))
stats.all <- full_join(stats.all, SVals, by = c("S"))

# Fit the separable model to the range ranges in Sebastian 2016.
stats.fit <- stats.all %>%
  filter(!(is.nan(response))) %>%
  group_by(TARGET) %>%
  nest() %>%
  mutate(models = map(data, function(x) {
    m1 <- mle2(response ~ dnorm(mean = k0 * (I(CVals) + b) * ((I(SVals) + c)^d), sd = 1/sqrt(N.norm)),
               start = list(k0 = 0, b = 0, c = 0, d = 0),
               data = x)
  })) %>%
  mutate(prediction = map(models, predict)) 

stats.fit$models

# Compute predictions for the data from the fitted model
stats.fit.predict <- stats.fit %>%
  unnest(data, prediction) %>%
  gather(type, value, c("response", "prediction"))

# Extract parameters from fitted model.
stats.fit.params <- stats.fit %>%
  mutate(parameters = map(models, . %>% coef %>% as.list() %>% as_tibble())) %>%
  unnest(parameters)

# Print the model parameters to disk.
stats.fit.params %>% 
  select(TARGET, k0, b, c, d) %>%
  kable(format = 'latex') %>%
  kable_as_image(filename = '~/Dropbox/Calen/Dropbox/gaussianized_params', file_format = 'png')  

# Create figure with scene statistics and predictions.
stat.fig <- ggplot(stats.fit.predict, aes(x = SVals, y = value, colour = as.factor(CVals))) +
  geom_point(data = . %>% filter(type == "response"), size = 2) +
  geom_line(data = . %>% filter(type == "prediction")) +
  facet_wrap(~LVals, ncol= 4) +
  theme_set(theme_gray(base_size = 30)) +
  list(theme(legend.title = element_text(size=10),
             plot.title = element_text(size = 10),
             axis.title = element_text(size=10),
             axis.text = element_text(size=10),
             legend.text = element_text(size=10))) +
  coord_cartesian(ylim = c(0, 10)) +
  theme(aspect.ratio = 1)

plot(stat.fig)

ggsave(filename = '~/Dropbox/Calen/Dropbox/gaussianized_weight_fig.pdf', width = 20, height = 20, device = 'pdf')
