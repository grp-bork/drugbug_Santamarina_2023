#!/usr/bin/env Rscript

pacman::p_load(tidyverse)
pacman::p_load(ggtext)
pacman::p_load(stringi)
pacman::p_load(glue)

source("config.R")

d_ <- read_tsv("data/ecoli_metabolomics.tsv")
d_ <- d_ %>% mutate(sample_type = stri_replace_all_fixed(sample_type, "MGAM", "mGAM"))

addOrUpdateWorksheet("Supplementary Table 3.xlsx", "E. coli metabolomics", d_, "Data underlying Fig. 4D: Bioaccumulation of niclosamide by E. coli")


d <- d_ %>% select(-starts_with("sd")) %>% pivot_longer(ends_with("uM")) %>% 
  mutate(sample_type = fct_rev(stri_replace_all_fixed(sample_type, "E. coli", "<i>E. coli</i>"))) %>% 
  mutate(name = stri_trans_totitle(stri_replace_all_fixed(name, "_uM", ""))) %>% 
  mutate(compartment_combined = ifelse(compartment == "Bioaccumulated", "bioaccumulated", "in supernatant/medium")) %>% 
  mutate(name_compartment = glue("{name} ({compartment_combined})")) %>% 
  mutate(name_compartment = factor(name_compartment, 
    levels = c("Niclosamide (in supernatant/medium)", "Niclosamide (bioaccumulated)", "Aminoniclosamide (bioaccumulated)", "Aminoniclosamide (in supernatant/medium)"))) %>% 
  mutate(Time_h = fct_rev(factor(Time_h))) 

d <- d %>% arrange(name_compartment) %>% group_by(sample_type, Time_h) %>% 
  mutate(x1 = cumsum(value), x0 = x1 - value)



dr <- d %>% filter(compartment == "Bioaccumulated") %>% group_by(sample_type, Time_h) %>% 
  summarise(x0 = min(x0), x1 = max(x1)) %>% 
  filter(x0!=x1)

dl <- bind_rows(
  d %>% filter(sample_type == "mGAM") %>% group_by(name_compartment) %>% filter(value == max(value)),
  dr %>% filter(x1-x0 == max(x1-x0)) %>% mutate(name = "Bioaccumulated")
) %>% mutate(value = x0+(x1-x0)/2)
  
d %>% ggplot(aes(value, as.integer(Time_h))) + 
  geom_rect(aes(xmin = x0, xmax = x1, ymin = as.integer(Time_h)-0.4, ymax = as.integer(Time_h)+0.4, fill=name)) + 
  geom_rect(data = dr, aes(xmin = x0, xmax = x1, ymin = as.integer(Time_h)-0.4, ymax = as.integer(Time_h)+0.4), alpha=0, fill="white", linewidth = 0.25, color = "black", inherit.aes = FALSE) + 
  geom_text(data = dl, aes(label = name), size = 5/(14/5), color = "white") +
  scale_x_continuous(name = "Compound concentration (µM)", expand = c(0,0), breaks = seq(0,10,2.5)) +
  scale_y_continuous(labels = levels(d$Time_h), breaks = 1:4, name = "Time (h)") +
  scale_fill_manual(values = c("#999999", "#444444")) +
  facet_wrap(~sample_type, ncol = 1) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_markdown(size = 5, color = "black"), 
    axis.title = element_text(size = 6, color = "black"),
    axis.ticks.length = unit(0, "pt"),
    legend.position = "none",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    strip.text = element_markdown(size = 6, margin = margin(0,0,0,0,unit="pt"))
  )

ggsave("panels/fig4_ecoli_metabolomics.png", width = 6, height = 4, units = "cm", dpi = 1000, bg = "white")
