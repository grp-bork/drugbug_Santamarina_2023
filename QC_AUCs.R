#!/usr/bin/env Rscript

pacman::p_load(tidyverse, stringi, reshape2, ggforce, cowplot)

d <- read_tsv("Single_spp_curves_new_Oct2019/aucs.tsv") %>% rename(concentration = conc)

amplicon_drugs <- read_tsv("data/16S_counts.tsv") %>% select(drug, concentration) %>% unique() %>% filter(concentration > 0)
tested_conditions <- read_tsv("Plates_Layout.tsv") %>% rename_all(stri_trans_tolower)  %>% select(drug, concentration) %>% unique() %>% semi_join(amplicon_drugs, by = "drug")
auc_conditions <- read_tsv("data/combined_monoculture_aucs.tsv") %>% select(drug, concentration = conc) %>% unique() %>% semi_join(tested_conditions)

d <- d %>% semi_join(auc_conditions, by = c("drug", "concentration"))

AUC_corr <- function(d) {
  stop()
  m <- d %>% acast(treatment ~ repl, value.var = "AUC") 
  ms <- cor(m, method = "s")
  mp <- cor(m, method = "p")
  tibble(corr_s = ms[ upper.tri(ms) ], corr_p = mp[ upper.tri(mp) ])
}


rmsd <- function(m) {
  sqrt(mean((m[1,] - m[2,])**2))  
}

AUC_rmsd <- function(d) {
  tibble(RMSD = rmsd(combn(d$AUC, 2)))
}

# d_biol <- 

d_biol <- d %>% mutate(repl = paste(experiment_number, biol_repl), treatment = paste(drug, concentration)) %>% 
  filter(!is.na(AUC)) %>% 
  group_by(repl, NT_code, treatment) %>% 
  summarise(AUC = median(AUC)) %>% 
  group_by(NT_code, treatment) %>% 
  filter(n()>1) %>% 
  group_modify(~ AUC_rmsd(.x))

median(d_biol$RMSD)
mean(d_biol$RMSD)

d_tech <- d %>% mutate(repl = paste(experiment_number, biol_repl), treatment = paste(drug, concentration)) %>% 
  filter(!is.na(AUC)) %>% 
  group_by(repl, NT_code, treatment) %>% 
  filter(n()>1) %>% 
  group_modify(~ AUC_rmsd(.x))

median(d_tech$RMSD)
mean(d_tech$RMSD)

