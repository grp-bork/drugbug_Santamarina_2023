#!/usr/bin/env Rscript

pacman::p_load(tidyverse, stringi, reshape2, ggforce, cowplot)

counts <- read_tsv("data/16S_counts.tsv")
effects <- read_tsv("data/effects.tsv") %>% filter(flavor == "absolute", !is.na(hit))

rel_abundance_corr <- function(d) {
  m <- d %>% acast(NT_code ~ repl, value.var = "rel_abundance") 
  ms <- cor(m, method = "s")
  mp <- cor(log(m+1e-6), method = "p")
  tibble(corr_s = ms[ upper.tri(ms) ], corr_p = mp[ upper.tri(mp) ])
}

abundance_ratio_corr <- function(d) {
  m <- d %>% acast(NT_code ~ repl, value.var = "abundance_ratio") 
  ms <- cor(m, method = "s")
  mp <- cor(m, method = "p")
  tibble(corr_s = ms[ upper.tri(ms) ], corr_p = mp[ upper.tri(mp) ])
}


d_control_biol <- counts %>% filter(!rare_species, control) %>% 
  mutate(repl = paste(run, plate_number, sample_id)) %>% 
  group_by() %>% 
  group_modify(~ rel_abundance_corr(.x))

d_control_tech <- counts %>% filter(!rare_species, control) %>% group_by(run) %>% 
  mutate(repl = paste(plate_number, sample_id)) %>% 
  group_modify(~ rel_abundance_corr(.x))

d_biol <- effects %>% filter(!rare_species) %>% group_by(drug, concentration) %>% 
  mutate(repl = paste(run, tech_repl)) %>% 
  group_modify(~ rel_abundance_corr(.x))
  
d_tech <- effects %>% filter(!rare_species) %>% group_by(drug, concentration, run) %>% 
  mutate(repl = tech_repl) %>% 
  group_modify(~ rel_abundance_corr(.x))

communityAUCs <- effects %>% filter(!control) %>% select(drug, concentration, run, biol_repl, tech_repl, communityAUC) %>% unique() %>% 
  group_by(drug, concentration) %>% summarise(communityAUC = mean(communityAUC))

median(d_biol$corr_s)
median(d_tech$corr_s)

median(d_control_biol$corr_s)
median(d_control_tech$corr_s)

dcorr <- bind_rows(
  d_biol %>% mutate(kind = "Drug treatment", repl = "Biological replicates"),
  d_control_biol %>% mutate(kind = "Control", repl = "Biological replicates"),
  d_tech %>% mutate(kind = "Drug treatment", repl = "Technical replicates"),
  d_control_tech %>% mutate(kind = "Control", repl = "Technical replicates")
)

dcorr %>% mutate(repl = stri_replace_all_fixed(repl, " ", "\n")) %>% ggplot(aes(repl, corr_s)) + 
  geom_boxplot(outlier.size = 0.5) + 
  facet_wrap(~kind, ncol = 2) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Spearman correlation between relative species abundances", limits = c(0.4, 1)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), axis.title = element_text(size = 6),
  ) 

ggsave("suppl_figures/QC_corrs.pdf", width = 5, height = 10, units = "cm")

