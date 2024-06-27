#!/usr/bin/env Rscript

pacman::p_load(tidyverse, stringi, reshape2, glue)

counts <- read_tsv("data/16S_counts.tsv")
effects <- read_tsv("data/effects.tsv") %>% filter(flavor == "absolute", !is.na(hit))

rel_abundance_corr <- function(d) {
  m <- d %>% acast(NT_code ~ repl, value.var = "rel_abundance") 
  ms <- cor(m, method = "s")
  mp <- cor(log(m+1e-6), method = "p")
  tibble(corr_s = ms[ upper.tri(ms) ], corr_p = mp[ upper.tri(mp) ])
}

d_control_biol <- counts %>% filter(!rare_species, control) %>% 
  group_by(run, plate_number, NT_code) %>% 
  summarise(rel_abundance = median(rel_abundance)) %>% 
  mutate(repl = paste(run, plate_number)) %>% 
  group_by() %>% 
  group_modify(~ rel_abundance_corr(.x))

d_control_tech <- counts %>% filter(!rare_species, control) %>% 
  mutate(repl = paste(plate_number, sample_id)) %>% 
  group_by(run) %>% 
  group_modify(~ rel_abundance_corr(.x))

d_biol <- effects %>% filter(!rare_species) %>% 
  group_by(drug, concentration, run, NT_code) %>% 
  summarise(rel_abundance = median(rel_abundance)) %>% 
  group_by(drug, concentration) %>% 
  mutate(repl = run) %>% 
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

dcorr <- dcorr %>% left_join(
  dcorr %>% group_by(kind, repl) %>% count() %>% mutate(label = glue("{repl} (N={n})"), label = stri_replace_all_fixed(label, " ", "\n"))
)

dcorr %>% 
  ggplot(aes(label, corr_s)) + 
  geom_boxplot(outlier.size = 0.5, varwidth = T) + 
  facet_wrap(~kind, ncol = 2, scales="free_x") +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Spearman correlation between relative species abundances", limits = c(0.4, 1)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), axis.title = element_text(size = 6),
  ) 

ggsave("suppl_figures/QC_corrs.pdf", width = 5, height = 7, units = "cm")
ggsave("suppl_figures/S1C_QC_corrs.png", width = 7, height = 7, units = "cm", dpi = 1000)

