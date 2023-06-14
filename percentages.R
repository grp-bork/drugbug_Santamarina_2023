
pacman::p_load(tidyverse)

source("config.R")

effects <- read_tsv("data/effects.tsv") %>% filter(flavor == "absolute", hit_def == "AUC")

effects <- effects %>% filter(!is.na(hit)) %>% 
  mutate( reduced_in_community = abundance_ratio < 1/2 ) %>% 
  inner_join(classification, by = c("hit", "reduced_in_community")) 

sum(effects$expected) / nrow(effects)
sum(!effects$expected & effects$reduced_in_community) / nrow(effects)
sum(!effects$expected & !effects$reduced_in_community) / nrow(effects)
