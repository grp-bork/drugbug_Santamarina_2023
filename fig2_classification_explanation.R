#!/usr/bin/env Rscript

pacman::p_load(tidyverse, cowplot, glue, stringi, egg, ragg, ggtext) 

source("config.R")

mean_effects_ <- read_tsv("data/mean_effects.tsv") %>% 
  filter(flavor == "absolute", hit_def == "AUC") %>% 
  mutate( reduced_in_community = abundance_ratio < 1/2, overgrowth = abundance_ratio > 2 ) %>% 
  mutate(species_abbr = stri_replace_all_regex(species, "^(.)[^ ]+ ([^ ]+).*", "$1. $2"))

mean_effects_ %>% ungroup() %>% filter(reduced_in_community) %>% arrange(-abundance_ratio)

mean_effects <- mean_effects_ %>% filter(drug == "Methotrexate", concentration == 5)

classification <- classification %>% 
  mutate(result = fct_rev(fct_inorder(result)))

percent_protected <- round(100*sum(!mean_effects$reduced_in_community & mean_effects$hit) / sum(mean_effects$hit), 0)
percent_sensitized <- round(100*sum(mean_effects$reduced_in_community & !mean_effects$hit) / sum(!mean_effects$hit), 0)

mean_effects %>% left_join(classification) %>% 
  filter(result == "sensitized in community") %>% 
  select(species, reduced_in_community, hit, AUC, abundance_ratio, result)

P_effect_y_min <- 0.005
P_effect_y_max <- 3

examples <- mean_effects %>% filter(species %in% c("Erysipelatoclostridium ramosum", "Escherichia coli ED1a", "Fusobacterium nucleatum subsp. nucleatum", "Veillonella parvula")) 

H <- 4

p1 <- mean_effects %>% left_join(classification) %>% 
  ggplot(aes(AUC, abundance_ratio)) +
  annotate("rect", xmin = 0, xmax = 0.75, ymin = 0.5, ymax = Inf, fill=COLOR_PROTECTION, alpha = 1) +
  annotate("rect", xmin = 0, xmax = 0.75, ymin = 0, ymax = 0.5, fill=COLOR_EXP_REDUCED, alpha = 1) +
  annotate("rect", xmin = 0.75, xmax = Inf, ymin = 0, ymax = 0.5, fill=COLOR_SENSITIZATION, alpha = 1) +
  annotate("rect", xmin = 0.75, xmax = Inf, ymin = 0.5, ymax = Inf, fill=COLOR_EXP_NORMAL, alpha = 1) +
  geom_segment(data = examples %>% filter(species_abbr != "F. nucleatum"), aes(xend = AUC + 0.05, yend = abundance_ratio), linewidth = 0.2) +
  geom_segment(data = examples %>% filter(species_abbr == "F. nucleatum"), aes(xend = AUC, yend = abundance_ratio * 0.77), linewidth = 0.2) +
  geom_point(shape = 16, size = 0.5) +
  geom_text(data = examples %>% filter(species_abbr != "F. nucleatum"), aes(label = species_abbr), fontface = "italic", hjust = 0, vjust = 0.5, size = 2, nudge_x = 0.05) +
  geom_text(data = examples %>% filter(species_abbr == "F. nucleatum"), aes(label = species_abbr, y = abundance_ratio * 0.75), fontface = "italic", hjust = 0.5, vjust = 1, size = 2) +
  scale_x_continuous(limits = c(-0.01, 1.3), name = "Ratio of growth curve AUCs in monoculture:<br>**5 µM Methotrexate** vs. DMSO control", breaks = c(0, 0.5, 1), labels = as.character, expand = c(0,0)) +
  scale_y_log10(name = "Ratio of normalized abundances in community:<br>**5 µM Methotrexate** vs. DMSO control", labels = as.character, limits = c(P_effect_y_min, P_effect_y_max), expand = c(0,0)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        legend.direction = "vertical",
        legend.title = element_blank(), 
        axis.text = element_text(size = 5), 
        axis.title.x = element_markdown(size = 6, lineheight = 1.2),
        # fudge top margin so that the y-axis title fits
        axis.title.y = element_markdown(size = 6, lineheight = 1.2, margin = margin(t = 25)),
        axis.ticks = element_line(color = "black", size = 0.2),
        legend.text = element_text(size = 5), 
        legend.spacing.x = unit(0, 'cm'), 
        legend.justification = "left",
        plot.title = element_text(size = 6, face = "bold"),
  )

p1

pp1 <- grid.arrange(set_panel_size(p1, width = unit(H, "cm"), height = unit(3.2, "cm")))

ggsave("panels/fig2_classification_explanation_color_bg.pdf", pp1, height = 4.5, width = 6, units = "cm")
ggsave("panels/fig2_classification_explanation_color_bg.png", pp1, height = 4.5, width = 6, units = "cm", dpi = 1000)


P_comm_x_min <- 0.001/2 
P_comm_x_max <- max(mean_effects$control_abundance) / 0.9
P_comm_y_min <- min(mean_effects$abundance) * 0.9

ribbon <- tibble(x = c(P_comm_x_min, P_comm_x_max)) %>% mutate(ymin = x / 2, ymax = x * 2)

P_comm_y_max <- max(ribbon$ymax)

breaks <- 10**(-4:-1)
breaks_ <- as.character(breaks)
  
W2 <- H * log(P_comm_x_max/P_comm_x_min) / log(P_comm_y_max/P_comm_y_min) 

p2 <- mean_effects %>% 
  ggplot(aes(control_abundance, abundance)) + 
  geom_abline(color = "grey") +
  geom_point(shape = 16, size = 0.5) +
  geom_text(data = examples, aes(label = species_abbr), fontface = "italic", 
            hjust = c(1.05, 1.1, 0.75, 0.5), vjust = c(0.5, 0.5, 1.2, 1.2), size = 2) +
  scale_x_log10(name = "Normalized species abundance in<br>DMSO control", limits = c(P_comm_x_min, P_comm_x_max), expand = c(0,0), breaks = breaks, labels = breaks_) + 
  scale_y_log10(name = "Normalized species abundance in<br>**5 µM Methotrexate** treatment", limits = c(P_comm_y_min, P_comm_y_max), expand = c(0,0), breaks = breaks, labels = breaks_) +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), 
        axis.text = element_text(size = 5), 
        axis.title.x = element_markdown(size = 6, lineheight = 1.2),
        axis.title.y = element_markdown(size = 6, lineheight = 1.2), 
        plot.title = element_text(size = 6, face = "bold")
  ) + 
  ggtitle("Growth in community") +
  NULL


pp2 <- grid.arrange(set_panel_size(p2, width = unit(W2, "cm"), height = unit(H, "cm")))

pp2

ggsave("panels/fig1_drug_vs_control_example.pdf", pp2, height = 6, width = 6, units = "cm")
ggsave("panels/fig1_drug_vs_control_example.png", pp2, height = 6, width = 5, units = "cm", dpi = 1000)

