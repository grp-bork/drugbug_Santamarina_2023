#!/usr/bin/env Rscript

pacman::p_load(ragg, tidyverse, egg, stringi, glue, ggtext)

counts <- read_tsv("data/counts.tsv")
counts <- counts %>% filter(flavor == "absolute") %>% 
  select(species, drug, concentration, abundance) %>% 
  mutate(species_abbr = stri_replace_all_regex(species, "^(.)[^ ]+ ([^ ]+).*", "$1. $2"))

example_drugs <- tribble(
  ~ drug, ~ concentration,
  "Chlorpromazine", 20,
  "Ciprofloxacin", 20,
  "Methotrexate", 5,
  "DMSO control", 0
) %>% mutate(
  label = ifelse(concentration > 0, glue("{concentration} µM {drug}"), drug ),
  label = ifelse(drug == "Methotrexate", glue("**{label}**"), label)
)

examples <- counts %>% semi_join(example_drugs, by = c("drug", "concentration"))

tp <- bind_rows(
  counts %>% filter(drug == "control") %>% mutate(drug = "DMSO control"),
  examples
) %>% group_by(drug, species_abbr) %>% summarise(abundance = mean(abundance))

species <- tp %>% group_by(species_abbr) %>% summarise(abundance = max(abundance)) %>% arrange(-abundance)
example_species <- species %>% filter(abundance > 0.1) %>% pull(species_abbr)

W <- 0.3

lbl_other <- glue("{ nrow(species) - length(example_species)} other species")

d <- tp %>% ungroup() %>% 
  mutate(species_abbr = fct_rev(fct_other(species_abbr, keep = example_species, other_level = lbl_other))) %>%
  group_by(species_abbr, drug) %>% summarise(abundance = sum(abundance)) %>% ungroup() %>% 
  mutate(drug = ordered(drug, levels = example_drugs$drug),
         y = as.integer(drug),
         ymin = y - W, ymax = y + W) %>%
  arrange(desc(species_abbr)) %>% 
  group_by(drug) %>% 
  mutate( xmax = cumsum(abundance), xmin = xmax - abundance)

ds <- d %>% arrange(drug) %>% group_by(species_abbr) %>% 
  mutate( x1 = xmax, x2 = lead(xmax), y1 = ymax, y2 = lead(ymin)) %>% 
  ungroup() %>% 
  select( x1, x2, y1, y2) %>% 
  filter(!is.na(x2))

dt <- 
  d %>% group_by(species_abbr) %>% filter(abundance > 0.5*max(abundance)) %>% 
    filter(drug == max(drug))

p <- d %>% 
  ggplot() +
  geom_segment(data = ds, aes(x = x1, xend = x2, y = y1, yend = y2), size = 0.2) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = species_abbr), color = "black", size = 0.2) +  
  geom_text(data = dt %>% filter(species_abbr != lbl_other), 
            aes(label = species_abbr, y = y, x = (xmin+xmax)/2), size = 2, fontface = "italic") +
  geom_text(data = dt %>% filter(species_abbr == lbl_other), 
            aes(label = species_abbr, y = y, x = (xmin+xmax)/2), size = 2) +
  scale_x_continuous(breaks = seq(0, 1, 0.5), name = "Normalised species abundance", limits = c(-0.01,1.05), expand = c(0,0)) +
  scale_y_continuous(breaks = 1:nrow(example_drugs), labels = example_drugs$label, name = "") +
  scale_fill_manual(values = c("#F2F2F2", scales::brewer_pal(type = "q", palette = "Pastel1", direction = -1)(7))) +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank(), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_markdown(size = 6, vjust = 0.6), # need to fudge to get baseline right
        axis.title = element_text(size = 6),
  ) 

pp <- grid.arrange(set_panel_size(p, width = unit(10, "cm"), height = unit(2, "cm")))

ggsave("panels/fig1_rel_abundance_barchart.pdf", pp, height = 15, width = 15, units = "cm")
ggsave("panels/fig1_rel_abundance_barchart.png", pp, height = 3.5, width = 13, units = "cm", dpi = 1000)
