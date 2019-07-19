library(dplyr)

load('~/Dropbox/Calen/Work/occluding/occlusion_detect/_data/model_wide.rdata')

scene_stats <- import_stats()
bin_values <- get_experiment_bin_values()
human_detect <- get_human_detect(human.responses) %>% select(TARGET, SUBJECT, eccentricity, dprime, BIN) %>% left_join(., bin_values, by = c("BIN", "TARGET")) %>% select(-statType, -statValue, -BIN)


scene_stats_tMatch <- scene_stats %>% filter(response_type %in% c("tMatch"), statistic == "sigma")
scene_stats_edge   <- scene_stats %>% filter(response_type %in% c("Eabs"), statistic == "sigma")
scene_stats_L   <- scene_stats %>% filter(response_type %in% c("L"), statistic == "sigma")

# Fit template stats
template.stats <- fit.template.stats(scene_stats_tMatch, c(1, 3, 5)) 

template.stats.1 <- template.stats %>% 
  rename(template_std = value) %>%
  filter(response == "prediction")

# Fit L stats
L.stats <- fit.L.stats(scene_stats_L, c(1, 3, 5)) 

L.stats.1 <- L.stats %>% 
  rename(L_std = value) %>%
  select(TARGET, L, C, S, Lvals, Cvals, Svals, response, ecc_deg, L_std) %>%
  filter(response == "prediction") %>% 
  select(-response)

# Fit edge stats
edge.stats <- fit.edge.stats(scene_stats, c(1,3,5))

edge.stats.1 <- edge.stats %>% 
  rename(edge_std = value) %>%
  select(TARGET, L, C, Lvals, Cvals, response, ecc_deg, edge_std) %>%
  filter(response == "prediction")


stats.joined <- left_join(template.stats.1, edge.stats.1, by = c("TARGET", "L", "C", "Cvals", "Lvals", "ecc_deg")) %>% select(-eccentricity) %>% select(TARGET, ecc_deg, L, C, S, Lvals, Svals, Cvals, template_std, edge_std)

stats.joined.1 <- left_join(stats.joined, L.stats.1, by = c("TARGET", "L", "C", "S", "Cvals", "Lvals", "Svals", "ecc_deg")) %>% select(TARGET, ecc_deg, L, C, S, Lvals, Svals, Cvals, template_std, edge_std, L_std)

stats.joined.2 <- stats.joined.1 %>% rename(eccentricity = ecc_deg) %>% group_by(TARGET, L, C,S,Lvals, Cvals, Svals) %>% nest()

human.stats.join <- inner_join(stats.joined.2, human_detect, by = c("TARGET", "L", "C", "S"))


yesno <- map(human.stats.join$data, function(x) {any(nrow(x) < 2)}) # cannot interpolate
idx <- which(yesno == TRUE)

human.stats.join <- human.stats.join[-idx, ]
  
estimates <- map2(human.stats.join$data, human.stats.join$eccentricity, function(x, y) {
 template_std_hat <- approxExtrap(x$eccentricity, y = x$template_std, y)$y
 edge_std_hat     <- approxExtrap(x$eccentricity, y = x$edge_std, y)$y
 L_std_hat     <- approxExtrap(x$eccentricity, y = x$L_std, y)$y
 
 data.frame(t_hat = template_std_hat, edge_hat = edge_std_hat, L_hat = L_std_hat)
}
)

human.stats.join$estimates <- estimates

human.stats.join.1 <- human.stats.join %>% unnest(estimates)



