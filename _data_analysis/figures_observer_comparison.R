load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/human.psychometrics.rdata")
load("~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/model.psychometrics.optimal.all.scaled.rdata")


h1.dat <- human.psychometrics %>% select(TARGET, SUBJECT, BIN, data)
m1.dat <- model.psychometrics.optimal.all.scaled %>% select(TARGET, BIN, sub_type, d0, e0, b)

c.dat <- right_join(h1.dat, m1.dat, by = c("TARGET", "BIN"))

c.dat.1 <- c.dat %>% unnest(data)


dprime.prediction <- c.dat.1 %>% mutate(dprime_hat = d0 * e0^b / (e0^b + eccentricity^b)) %>% filter(!is.infinite(dprime) & !is.infinite(dprime_hat))

dprime.prediction %>% group_by(TARGET, sub_type) %>% summarize(cor(dprime, dprime_hat)^2)

ggplot(dprime.prediction, aes(x = scale(dprime), y = scale(dprime_hat), colour = eccentricity)) + 
  geom_point() + 
  theme(aspect.ratio = 1) + 
  facet_grid(TARGET~sub_type)



scale.dprime <- optimal.all %>% group_by(sub_type) %>% select(BIN, TARGET, eccentricity, dprime, sub_type) %>% mutate(dprime_scaled = scale(dprime))

spread.scale <- scale.dprime %>% select(-dprime) %>% spread(sub_type, dprime_scaled) %>% gather("type", "value", c(4,6,7))

compute.cor <- spread.scale %>% group_by(type) %>% summarize(corr.val = cor(optimal, value))

spread.scale.1 <- spread.scale %>% left_join(., compute.cor, by = c("type"))

ggplot(spread.scale.1, aes(x = optimal, y = value, colour = type)) + geom_point() + expand_limits(x = c(0,1), y = c(0,1)) + theme(aspect.ratio = 1) + geom_abline()
