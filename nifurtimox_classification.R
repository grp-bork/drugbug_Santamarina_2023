#!/usr/bin/env Rscript

pacman::p_load(tidyverse, cowplot, glue, stringi, egg, ragg, ggrepel, ggtext) 

source("config.R")

mean_effects_ <- read_tsv("data/mean_effects.tsv") %>% 
  filter(flavor == "absolute", hit_def == "AUC") %>% 
  mutate( reduced_in_community = abundance_ratio < 1/2, overgrowth = abundance_ratio > 2 ) %>% 
  mutate(species_abbr = stri_replace_all_regex(species, "^(.)[^ ]+ ([^ ]+).*", "$1. $2"))

mean_effects_ %>% ungroup() %>% filter(reduced_in_community) %>% arrange(-abundance_ratio)

mean_effects <- mean_effects_ %>% filter(drug == "Nifurtimox", concentration == 20)

classification <- classification %>% 
  mutate(result = fct_rev(fct_inorder(result)))

percent_protected <- round(100*sum(!mean_effects$reduced_in_community & mean_effects$hit) / sum(mean_effects$hit), 0)
percent_sensitized <- round(100*sum(mean_effects$reduced_in_community & !mean_effects$hit) / sum(!mean_effects$hit), 0)

mean_effects %>% left_join(classification) %>% 
  filter(result == "sensitized in community") %>% 
  select(species, reduced_in_community, hit, AUC, abundance_ratio, result)

P_effect_y_min <- 0.08
P_effect_y_max <- 12

examples_black <- mean_effects %>% 
  filter(species %in% c("Bacteroides thetaiotaomicron", "Coprococcus comes", "Parabacteroides merdae", "Phocaeicola vulgatus", "Roseburia intestinalis")) 

examples_red <- mean_effects %>% 
  filter(species %in% c("Streptococcus parasanguinis", "Streptococcus salivarius", "Escherichia coli ED1a"))

examples_black <- examples_black %>% arrange(abundance_ratio) %>% 
  mutate(need_line = n()>1, y = 10**(log10(min(abundance_ratio)) + (row_number()-1)/8) + 0.1, x = 0.15)

  
H <- 4

p1 <- mean_effects %>% left_join(classification) %>% 
  ggplot(aes(AUC, abundance_ratio)) +
  annotate("rect", xmin = 0, xmax = 0.75, ymin = 0.5, ymax = Inf, fill=COLOR_PROTECTION, alpha = 1) +
  annotate("rect", xmin = 0, xmax = 0.75, ymin = 0, ymax = 0.5, fill=COLOR_EXP_REDUCED, alpha = 1) +
  annotate("rect", xmin = 0.75, xmax = Inf, ymin = 0, ymax = 0.5, fill=COLOR_SENSITIZATION, alpha = 1) +
  annotate("rect", xmin = 0.75, xmax = Inf, ymin = 0.5, ymax = Inf, fill=COLOR_EXP_NORMAL, alpha = 1) +
  geom_segment(data = examples_black %>% filter(need_line), aes(xend = x, yend = y), size = 0.2) +
  geom_point(shape = 16, size = 0.5) +
  geom_text(data = examples_black, aes(x = x, y = y, label = species_abbr), fontface = "italic", hjust = 0, vjust = 0.5, size = 2) +
  geom_text(data = examples_red %>% filter(species_abbr == "E. coli"), aes(label = species_abbr), fontface = "italic", hjust = 0, vjust = 0.5, size = 2, nudge_x = 0.02, color = "red") +
  geom_text(data = examples_red %>% filter(species_abbr != "E. coli"), aes(label = species_abbr), fontface = "italic", hjust = 0.25, vjust = 1.2, size = 2, color = "red") +
  scale_x_continuous(limits = c(-0.01, 1.3), name = "Ratio of growth curve AUCs in monoculture:<br><b>Nifurtimox (20 µM)</b> vs. control", breaks = c(0, 0.5, 1), labels = as.character, expand = c(0,0)) +
  scale_y_log10(name = "Ratio of abundances in community:<br><b>Nifurtimox (20 µM)</b> vs. control", labels = as.character, limits = c(P_effect_y_min, P_effect_y_max), expand = c(0,0)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        legend.direction = "vertical",
        legend.title = element_blank(), 
        axis.text = element_text(size = 5), 
        axis.title.x = element_markdown(size = 6, lineheight = 1.2),
        axis.title.y = element_markdown(size = 6, lineheight = 1.2),
        axis.ticks = element_line(color = "black", size = 0.2),
        legend.text = element_text(size = 5), 
        legend.spacing.x = unit(0, 'cm'), 
        legend.justification = "left",
        # plot.title.position = "plot", 
        plot.title = element_text(size = 6, face = "bold"),
  )

p1

pp1 <- grid.arrange(set_panel_size(p1, width = unit(H, "cm"), height = unit(4.5, "cm")))

ggsave("suppl_figures/nifurtimox_classification.png", pp1, height = 6, width = 5.7, units = "cm", dpi = 1000)
