#!/usr/bin/env Rscript

pacman::p_load(ragg, tidyverse, stringi, egg, scales, ggtext)

GCsf <- read_tsv("data/single_species_growth_curves.tsv", col_types = c("strain" = "c"))

examples <- c("Erysipelatoclostridium ramosum", "Escherichia coli ED1a", "Fusobacterium nucleatum subsp. nucleatum", "Veillonella parvula")

treatment <- GCsf %>% filter(drug == "Methotrexate", conc == 5, !discard_conc, species_name %in% examples) 

controls <- GCsf %>% filter(control) %>% 
  semi_join(treatment, by = c("drug_plate", "tech_repl", "experiment_number", "biol_repl", "plate_number", "species_name"))

d <- bind_rows(treatment, controls) %>% filter(biol_repl == "biol1") %>% 
  mutate(species_abbr = stri_replace_all_regex(species_name, "^(.)[^ ]+ ([^ ]+).*", "$1. $2"))

d <- d %>%
  group_by(time, control, species_abbr) %>% summarise(ODc01 = median(ODc01)) %>%
  group_by(control, species_abbr) %>% mutate(ODc01 = ODc01 - ODc01[1])


p <- d %>% 
  ggplot(aes(time, ODc01, color = control, linetype = control)) + 
  geom_line(data = d %>% filter(!control)) + 
  geom_line(data = d %>% filter(control)) + 
  facet_wrap(~species_abbr, scales = "free_x", nrow = 2) +
  scale_color_manual(values = c("black", "darkgrey"), labels = c("**5 µM Methotrexate treatment**", "DMSO control")) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("**5 µM Methotrexate treatment**", "DMSO control")) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), name = "Normalized OD", expand = c(0.01,0)) +
  scale_x_continuous(breaks = breaks_pretty(2), name = "Time (h)", expand = c(0,0)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), 
        strip.text = element_text(size = 6, face = "italic"),
        axis.text = element_text(size = 6), axis.title = element_text(size = 6),
        legend.text = element_markdown(size = 6), 
        title = element_text(size = 7),
        legend.box.margin = margin(-10,0,0,0),
        legend.spacing.y = unit(0, 'cm'), 
        legend.direction = "vertical",
        legend.justification = "center",
        legend.key.height = unit(0.5, "line"),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 6, face = "bold")
  ) +
  ggtitle("Growth in monoculture") +
  NULL


pp <- grid.arrange(set_panel_size(p, width = unit(6/4, "cm"), height = unit(1, "cm")))

ggsave("panels/fig1_single_species_growth_curves.pdf", pp, height = 6, width = 4.5, units = "cm")
ggsave("panels/fig1_single_species_growth_curves.png", pp, height = 6, width = 4.5, units = "cm", dpi = 1000)

