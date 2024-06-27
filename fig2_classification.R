#!/usr/bin/env Rscript

pacman::p_load(tidyverse)
pacman::p_load(glue)
pacman::p_load(stringi)
pacman::p_load(ggtext)

source("config.R")

effects <- read_tsv("data/effects.tsv")
effect_counts_hit <- read_tsv("data/effect_counts_hit.tsv")
species_effect_counts_hit <- read_tsv("data/species_effect_counts_hit.tsv")

effect_stats <- effect_counts_hit %>% 
  ungroup() %>% 
  filter(flavor == "absolute") %>% 
  group_by(flavor, drug, concentration, category, result, hit, expected) %>% 
  summarise(fraction = sum(fraction)/n_repl[1], n = sum(n)/n_repl[1]) %>% 
  ungroup()

bind_rows(
  effect_stats,
  effect_stats %>% select(-category, -result) %>% filter(fraction == 1) %>% mutate(expected = !expected, fraction = 0) %>% left_join(classification)
) %>% group_by(result) %>% summarise(mean(fraction))

effect_stats %>% group_by(category, result, expected) %>% summarise(n = sum(n)) %>% group_by(category) %>% mutate(f = n/sum(n), total = sum(n) ) 
effect_stats %>% group_by(expected) %>% summarise(n = sum(n)) %>% ungroup() %>% mutate(f = n/sum(n), total = sum(n) )

effect_stats %>% filter(drug %in% c("Entacapone", "Niclosamide", "Nifurtimox")) %>% 
  group_by(drug) %>% filter(concentration == max(concentration)) %>% 
  group_by(category, result, expected) %>% summarise(n = sum(n)) %>% group_by(category) %>% mutate(f = n/sum(n), total = sum(n) ) 

effect_stats %>% filter(drug %in% c("Entacapone", "Niclosamide", "Nifurtimox")) %>% 
  group_by(drug) %>% filter(concentration == max(concentration)) %>% 
  group_by(hit) %>% summarise(n = sum(n)) %>% ungroup() %>% mutate(f = n/sum(n), total = sum(n) ) 
  
effect_stats %>% filter(drug %in% c("Entacapone", "Niclosamide", "Nifurtimox")) %>% 
  group_by(drug) %>% filter(concentration == max(concentration)) %>% 
  group_by(drug, hit) %>% summarise(n = sum(n)) %>% group_by(drug) %>% mutate(f = n/sum(n), total = sum(n) ) 

effects %>% select(species, drug, concentration, hit) %>% unique() %>% filter(is.na(hit))

d <- effect_stats %>% select(-hit, -category, -expected, -flavor) %>% pivot_wider(names_from = "result", values_from = c("fraction", "n"))

d <- d %>% mutate_if(is.double, function(x) replace_na(x, 0))

d <- d %>% mutate(`fraction_as expected` = (`n_expected (normal growth in both cases)` + `n_expected (normal growth in both cases)`) / 
                    (`n_expected (normal growth in both cases)` + `n_expected (normal growth in both cases)` + `n_protected in community` + `n_sensitized in community`))

d <- d %>% select(drug, concentration, `fraction_protected in community`, `fraction_sensitized in community`, `fraction_as expected`,
                  `n_protected in community`, `n_expected (reduced growth in both cases)`, `n_sensitized in community`, `n_expected (normal growth in both cases)`) 

d %>% write_tsv("data/table_sensitization_protection.tsv")


all_gut_conc <- read_tsv("data/closest_to_gut_concentration.tsv")

closest_conc <- all_gut_conc %>% semi_join(effect_stats) %>% arrange(abs_log_ratio) %>% 
  group_by(drug) %>% filter(row_number() == 1)

effect_stats %>% semi_join(closest_conc %>% filter(abs_log_ratio < log(3))) %>% 
  group_by(category, result) %>% 
  summarise(n = sum(n)) %>% 
  group_by(category) %>% 
  mutate(f = n / sum(n))

effect_stats %>% 
  group_by(category, result) %>% 
  summarise(n = sum(n)) %>% 
  group_by(category) %>% 
  mutate(f = n / sum(n), fp = round(100*f))

effect_stats %>% 
  group_by(expected) %>% 
  summarise(n = sum(n)) %>% 
  mutate(f = n / sum(n), fp = round(100*f))


# need to take care to include all classifications to get the average by drug right
crossing(effect_stats %>% select(drug, concentration) %>% unique(), classification) %>% 
  left_join(effect_stats %>% select(-flavor)) %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% 
  group_by(drug, concentration, category, result) %>% 
  summarise(n = sum(n)) %>% 
  group_by(drug, concentration, category) %>% 
  filter(sum(n) > 0) %>% 
  mutate(f = n / sum(n)) %>% 
  group_by(category, result) %>% 
  summarise(mean(f), median(f))
  

drug_levels <- effect_stats %>% semi_join(closest_conc) %>% 
  group_by(drug) %>% summarise(n_total = sum(hit*n)/sum(n), n_protected = sum(hit*n*!expected)/sum(n)) %>% 
  arrange(n_total, n_protected) %>% pull(drug) %>% as.character()

tp <- 
  effect_stats %>% semi_join(closest_conc) %>% 
  mutate(drug = factor(as.character(drug), levels = drug_levels), 
         y = as.integer(drug), 
         result = factor(result, levels = rev(classification$result))) %>% 
  arrange(result) %>% 
  group_by(drug) %>% 
  mutate(
    f = n / sum(n),
    x2 = cumsum(f) + 0.4/21*!hit,
    x1 = x2 - f
  ) 

tl <- bind_rows(
  tp %>% filter(!expected) %>% mutate(label = as.character(glue("{round(100*fraction)}%"))),
  tp %>% group_by(drug, category) %>% filter( n() == 1) %>% mutate(label = "0%")
)

to <- tp %>% group_by(y, category) %>% summarise(x1 = min(x1), x2 = max(x2))

B <- seq(0, 1, length.out = 3)


dl <- tl %>% ungroup() %>% select(drug, concentration, y) %>% unique() %>% arrange(drug) %>% 
    left_join(closest_conc) %>% 
    mutate(
      label = ifelse(drug == "Methotrexate",
                 paste0("<b>", drug, " (", concentration, " µM)</b>"),
                 paste0(drug, " (", concentration, " µM)")),
      label = 
           ifelse(!is.na(ratio) & abs_log_ratio < log(3), 
            paste0(label, "*"),
            paste0(label, "<span style='color:white'>*</i>"))
        )
            
p <- tp %>% 
  ggplot(aes(xmin = x1, xmax=x2, ymin = y - 0.4, ymax = y+0.4)) +
  geom_rect(aes(fill = result)) +
  geom_rect(data = to, fill=alpha("grey",0), color = "black", size = 0.2) +
  geom_text(data = tl %>% filter(hit), aes(x = -0.1/21, y = y, label = label), hjust = 1, inherit.aes = F, size = 5/14*5) +
  geom_text(data = tl %>% filter(!hit), aes(x = x2+4/21, y = y, label = label), hjust = 1, inherit.aes = F, size = 5/14*5) +
  annotate("text", x = -3.6/21, y = max(tp$y)+1, label = "Protection", hjust = 0, size = 2) +
  annotate("text", x = 25/21, y = max(tp$y)+1, label = "Sensitization", hjust = 1, size = 2) +
  scale_fill_manual(values = c(COLOR_PROTECTION, COLOR_EXP_REDUCED, COLOR_EXP_NORMAL, COLOR_SENSITIZATION)) +
  scale_y_continuous(labels = dl$label, breaks = seq_along(levels(tp$drug)), name = "", expand = c(0,0), limits = c(0.5, max(tp$y)+1.5)) +
  scale_x_continuous(limits = c(-3.7/21, 27/21), breaks = B, labels = paste0(B*100, "%"),  expand = c(0,0), name = "Percentage of species") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    legend.direction = "vertical",
    axis.text.x = element_text(size = 5), 
    axis.text.y = element_markdown(size = 5, vjust = 0.6), # need to fudge to get baseline right
    axis.title = element_text(size = 6),
    axis.ticks.length = unit(0, "pt")
  ) 

ggsave("panels/fig2_classification_closest_conc_barchart.png", p, width = 7, height = 8, dpi = 1000, bg = "white", units = "cm")
ggsave("panels/fig2_classification_closest_conc_barchart.pdf", p, width = 8, height = 8, units = "cm")

## all conditions

metab_conditions <- read_tsv("data/metabolomics_conditions.tsv")

tp <- effect_stats %>% 
  mutate(result = factor(result, levels = rev(classification$result))) %>% 
  arrange(result) %>% 
  group_by(drug) %>% 
  mutate(y = 3-as.integer(factor(concentration))) %>% 
  group_by(drug, concentration) %>% 
  mutate(
    f = n / sum(n),
    x2 = cumsum(f),
    x1 = x2 - f,
  ) 

tl <- bind_rows(
  tp %>% filter(!expected) %>% mutate(label = as.character(glue("{round(100*fraction)}%"))),
  tp %>% group_by(drug, concentration, category) %>% filter( n() == 1, fraction < 1) %>% mutate(label = "0%")
)

to <- tp %>% group_by(y, category) %>% summarise(x1 = min(x1), x2 = max(x2))

conc_labels <- tp %>% select(drug, concentration, y) %>% 
  unique() %>% left_join(all_gut_conc) %>% 
  mutate(within_factor_3 = ifelse(is.na(abs_log_ratio), FALSE, abs_log_ratio < log(3)),
         label = paste0(concentration, " µM"))

conc_labels %>% semi_join(metab_conditions)


tp %>% 
  ggplot(aes(xmin = x1, xmax=x2, ymin = y - 0.4, ymax = y+0.4)) +
  geom_rect(aes(fill = result)) +
  geom_rect(data = conc_labels %>% semi_join(metab_conditions),
            xmin = -0.25, xmax = 1, fill=alpha("grey",0), color = "black") +
  geom_text(data = conc_labels %>% filter(!within_factor_3),
            aes(x = -0.02, y = y, label = label), hjust = 1, inherit.aes = F, size = 2) +
  geom_text(data = conc_labels %>% filter(within_factor_3), 
            aes(x = -0.02, y = y, label = label), hjust = 1, inherit.aes = F, size = 2,
            fontface = "bold") +
  
  geom_text(data = tl %>% filter(hit), 
            aes(x = 0.02, y = y, label = label), hjust = 0, inherit.aes = F, size = 2) +
  geom_text(data = tl %>% filter(!hit), 
            aes(x = 0.98, y = y, label = label), hjust = 1, inherit.aes = F, size = 2) +
  scale_fill_manual(values = c(COLOR_PROTECTION, COLOR_EXP_REDUCED, COLOR_EXP_NORMAL, COLOR_SENSITIZATION)) +
  scale_x_continuous(
    limits = c(-0.26, 1.01), breaks = B, labels = paste0(B*100, "%"),  expand = c(0,0), name = "Fraction of 21 species across replicates") +
  facet_wrap(~drug, ncol = 4) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    strip.text = element_text(hjust = 0, size = 8),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
  ) 

ggsave("figures/figS2.png", width = 18, height = 25, units = "cm", dpi = 1000, bg="white")


## taxonomy-based analysis -- without significant results!

species_effect_stats <- species_effect_counts_hit %>% 
  ungroup() %>% 
  filter(flavor == "absolute") %>% 
  group_by(flavor, species, category, result, hit, expected) %>% 
  summarise(fraction = sum(fraction)/n_repl[1], n = sum(n)/n_repl[1]) %>% 
  ungroup() 


species_clades <- read_tsv("data/species_clades.tsv")
species_clades <- species_clades %>% semi_join(species_effect_stats, by = c(Species = "species"))

d <- species_effect_stats %>% left_join(species_clades, by = c(species = "Species")) %>% 
  filter(hit, !expected)

kruskal.test(fraction ~ Phylum, d)  
kruskal.test(fraction ~ Class, d)  
kruskal.test(fraction ~ Order, d)  
kruskal.test(fraction ~ Family, d)  
kruskal.test(fraction ~ Genus, d)  

## various correlations between species abundance and fraction of effects

dh <- effects %>%
  filter(flavor == "absolute") %>%
  group_by(species, drug, concentration, hit) %>% summarise(abundance_ratio = median(abundance_ratio)) %>%
  mutate( reduced_in_community = abundance_ratio < 1/2 ) %>%
  inner_join(classification, by = c("hit", "reduced_in_community")) %>%
  mutate(result = factor(result, levels = rev(classification$result)))





dh %>% ungroup() %>%
  mutate(g = paste(drug, concentration)) %>%
  mutate(g = fct_reorder(g, hit, mean), species = fct_reorder(species, hit, mean)) %>%
  group_by(g, category) %>% summarise(fraction_unexpected = 1-mean(expected)) %>%
  ggplot(aes(g, fraction_unexpected, color = category)) + geom_point()

dh %>% ungroup() %>%
  mutate(species = fct_reorder(species, hit, mean)) %>%
  group_by(species, category) %>% summarise(fraction_unexpected = 1-mean(expected)) %>%
  ggplot(aes(species, fraction_unexpected, color = category)) + geom_point()


ca <- effects %>%
  filter(flavor == "absolute") %>%
  select(species, control_abundance) %>% unique() %>% group_by(species) %>%
  summarise(control_abundance = median(control_abundance))

d1 <- dh %>% 
  group_by(species, category) %>% summarise(n_total = n(), n_unexpected = sum(!expected), f = n_unexpected / n_total) %>% 
  left_join(ca) %>% ungroup() %>% 
  mutate(species = fct_reorder(species, control_abundance)) 

pvs <- d1 %>% group_by(category) %>% 
  summarise(r_s = cor(f, control_abundance, method = "s"))

d1 %>% ggplot(aes(control_abundance, f, color = category)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "grey", se = F) +
  geom_point(aes(size = n_total), alpha = 0.5, shape = 16) +
  geom_text(data = pvs, aes(label = paste0("rs = ", round(r_s, 2))), x = Inf, y = Inf, color = "black", hjust = 1, vjust = 1, size = 3) +
  scale_x_log10("Relative species abundance in control") +
  scale_y_continuous(name = "Fraction of treatment conditions", expand = c(0.1,0)) +
  scale_size_area() +
  scale_color_manual(values = c(COLOR_PROTECTION, COLOR_SENSITIZATION)) +
  facet_wrap(~category, ncol = 1, scales = "free") + 
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        legend.direction = "vertical",
        legend.title = element_blank(), 
        axis.text = element_text(size = 5), axis.title = element_text(size = 6),
        legend.text = element_text(size = 5), 
        legend.spacing.x = unit(0, 'cm'), 
        legend.justification = "left",
        plot.title.position = "plot", plot.title = element_text(size = 6, face = "bold"),
  ) 

ggsave("panels/fig2_classification_species_vs_abundance.pdf", width = 4, height = 4)
ggsave("panels/fig2_classification_species_vs_abundance.png", width = 4, height = 4, dpi = 1000, bg = "white")



  
du <- effects %>% 
  filter(flavor == "absolute") %>%
  mutate( reduced_in_community = abundance_ratio < 0.85 ) %>% 
  inner_join(classification, by = c("hit", "reduced_in_community")) %>% 
  ungroup() %>% modify_if(is.character, as.factor) %>% 
  group_by(drug, concentration, expected, .drop = F) %>% count() %>% 
  semi_join(classification) %>% 
  group_by(drug, concentration) %>% 
  summarise( n_unexpected = n[!expected], n_total = sum(n), fraction = n_unexpected/n_total ) %>% 
  inner_join(
    effects %>% 
      filter(flavor == "absolute") %>%
      select(drug, concentration, species, hit) %>% unique() %>% 
      group_by(drug, concentration) %>% summarise(n_hit = sum(hit)) %>% 
      ungroup() %>% modify_if(is.character, as.factor),
    by = c("drug", "concentration")
  )

du %>% ggplot(aes(n_hit, fraction)) + 
  geom_line(aes(group = drug), alpha = 0.2) +
  geom_point(size = 1) +
  xlab("Number of affected species per monoculture treatment condition") +
  ylab("Fraction of unexpected outcomes") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        legend.direction = "vertical",
        legend.title = element_blank(), 
        axis.text = element_text(size = 5), axis.title = element_text(size = 6),
        legend.text = element_text(size = 5), 
        legend.spacing.x = unit(0, 'cm'), 
        legend.justification = "left",
        plot.title.position = "plot", plot.title = element_text(size = 6, face = "bold"),
  ) 

ggsave("panels/fig2_classification_fraction_unexpected.pdf", width = 4, height = 4)
ggsave("panels/fig2_classification_fraction_unexpected.png", width = 4, height = 4, dpi = 1000, bg = "white")

