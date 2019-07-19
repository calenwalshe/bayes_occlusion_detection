load('~/Dropbox/Calen/Dropbox/template.response.rdata')
pink.templates <- get_template_response("~/Dropbox/Calen/Dropbox/pink_templates.txt")

all.templates.1 <- template.response %>% group_by(PYRAMIDLVL, BIN, TARGET, TPRESENT, function_name) %>% summarize(TRESP.natural = mean(TRESP))
pink.templates.1 <- white.templates %>% group_by(PYRAMIDLVL, BIN, TARGET, TPRESENT, function_name) %>% summarize(TRESP.white = mean(TRESP))

joined.templates <- left_join(all.templates.1, white.templates.1, by = c("PYRAMIDLVL", "BIN", "TARGET", "TPRESENT", "function_name"))

joined.templates.1 <- joined.templates %>% filter(TPRESENT == "present", function_name == "edge_cos")

ggplot(joined.templates.1 %>% filter(TPRESENT == "present", function_name == "edge_cos"), aes(x = TRESP.white, y = TRESP.natural, colour = BIN)) + geom_point() + theme(aspect.ratio = 1) + facet_grid(PYRAMIDLVL~TARGET) + expand_limits(x = c(.65, 1), y = c(.65, 1))
