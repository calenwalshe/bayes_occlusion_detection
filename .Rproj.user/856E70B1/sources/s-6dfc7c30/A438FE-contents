library(dplyr)
library(data.table)
library(purrrlyr)
library(purrr)
library(parallel)
library(tidyr)
library(Rfast)

px.max <- 2^14 -1

response1d <- fread("/main/calen/occluding/1dresponse/1Dresponse.txt", header = T)

TRESP.split <- stringr::str_split(response1d$TRESP, "\\s+")
TRESP.split <- map(TRESP.split, function(x) as.numeric(x))
TRESP.split <- map(TRESP.split, function(x) x[(512/2 - 80):(512/2 + 80)])

response1d$TRESP <- TRESP.split

response1d <- response1d %>% select(BIN, TARGET, PYRAMIDLVL, TPRESENT, TRESP)

response1d <- response1d %>% group_by(BIN, TARGET, PYRAMIDLVL, TPRESENT) %>% nest()

response1d$TRESP <- map(response1d$data, function(x) do.call(rbind, x[[1]])/px.max)

viz <- response1d %>% filter(BIN == 1, TARGET == 2, PYRAMIDLVL == 1, TPRESENT == 2)


plot(colMeans(viz$TRESP[[1]]))
plot(colVars(viz$TRESP[[1]]))

save(file = '/main/calen/occluding/response1d.rdata', response1d)
