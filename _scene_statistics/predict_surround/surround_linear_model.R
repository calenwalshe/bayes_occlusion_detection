# Import Libraries
library(purrrlyr)
library(tidyr)
library(purrr)
library(multidplyr)
library(ggplot2)
#

# Parse File Names for Conditions and Target information.
setwd('/main/calen/occluding/surround_prediction/')
files <- list.files('./', full.names = T)

parse.files <- lapply(files, FUN = function(x) {
  stat <- regmatches(x, regexpr("(?<=stat_)[[:alnum:]]+", x, perl = TRUE))
  target <- regmatches(x, regexpr("(?<=target_)[[:alnum:]]+", x, perl = TRUE))
  L_bin <- regmatches(x, regexpr("(?<=L)[[:alnum:]]+", x, perl = TRUE))
  C_bin <- regmatches(x, regexpr("(?<=C)[[:alnum:]]+", x, perl = TRUE))
  S_bin <- regmatches(x, regexpr("(?<=S)[[:alnum:]]+", x, perl = TRUE))
  
  data.frame(stat = stat, target = target, L = L_bin, C = C_bin, S = S_bin, path = as.character(x))
})

parse.files     <- do.call(rbind, parse.files)
parse.files$path <- as.character(parse.files$path)
#

# Perform regression by bin and target
cluster <- create_cluster(12) %>% cluster_library(list("purrr", "dplyr", "tidyr", "purrrlyr"))

regression.bin <- parse.files %>%
  filter(target == 1) %>%
  partition(cluster = cluster,target, stat, L) %>%
  mutate(model = map(path, function(x) {
    bin.dat        <- read.csv(file = x)
    names(bin.dat) <- letters[1:ncol(bin.dat)]
    m.1            <- lm(a ~ 0 + b + c + d + e + f + g + h + i, data = bin.dat)
    })) %>%
  collect()

regression.bin.rsq <- regression.bin %>% 
  mutate(rsq = map(model, function(x) summary(x)$r.squared)) %>%
  ungroup %>%
  mutate(L = factor(L, levels = 1:10), C = factor(C, levels = 1:10), S = factor(S, levels = 1:10)) %>%
  unnest(rsq)

#

# Temporary Figure to look at natural scenes
target.1 <- regression.bin.rsq %>% 
  filter(target == 1, stat == "L")

fig.1 <- ggplot(target.1, aes(x = L, y = rsq, colour = C, linetype = stat)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~S, ncol = 4) +
  theme(aspect.ratio = 1)

ggsave(filename = '~/Dropbox/Calen/Dropbox/rsq.pdf', plot = fig.1, device = 'pdf')

plot(fig.1)
#
