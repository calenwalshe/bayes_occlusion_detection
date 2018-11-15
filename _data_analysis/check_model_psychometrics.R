

model.dprime.1 <- model.dprime %>% mutate(dprime = dprime * .1)

model.psychometrics <- get_model_psychometric(model.dprime, 1)

plot.model.psychometric(model.psychometrics, model.dprime.1)
