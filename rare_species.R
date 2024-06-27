pacman::p_load(tidyverse, stringi, here, reshape2, scales, glue)

counts <- read_tsv("data/16S_counts.tsv")

control_counts <- counts %>% filter(control) 

detection_limits <- counts %>% group_by(sample_id, run, biol_repl, tech_repl, drug, concentration) %>% 
  summarise(detection_limit = 10 / sum(count)) 

counts <- counts %>% filter(!control, rare_species)

counts <- counts %>% semi_join( read_tsv("data/combined_monoculture_aucs.tsv") %>% select(drug, concentration = conc) %>% unique(), by = c("drug", "concentration"))


computePV_wilcox <- function(d) {
  tibble(
    p.value = wilcox.test(
      control_counts %>% semi_join(d, by = c("species")) %>% pull(absolute_abundance),
      d %>% pull(absolute_abundance),
      alternative = "l"
    )$p.value
  )
}

pvs <- counts %>% group_by(species, drug, concentration) %>% do(computePV_wilcox(.))

pvs %>% write_tsv("suppl_tables/rare_species_pvs.tsv")

median_pvs <- pvs %>% group_by(species, drug, concentration) %>% summarise(p.value = median(p.value)) %>% arrange(p.value)

dd <- counts %>% semi_join(median_pvs %>% ungroup() %>% filter(p.value < 0.001) %>% select(species, drug)) %>% 
  group_by(species)

median_pvs <- median_pvs %>% mutate(p.value_ = sprintf("%.1g", p.value))

plotSpecies <- function(d) {
  d <- bind_rows(
    d, 
    crossing(d %>% ungroup() %>% select(drug), control_counts %>% semi_join(d, by = c("species")) %>% select(-drug))
  ) %>% arrange(concentration) %>% mutate(concentration_ = factor(concentration)) 
  
  dl <- d %>% group_by(drug, concentration_) %>% filter(absolute_abundance == max(absolute_abundance)) %>% 
    inner_join(median_pvs, by = c("species", "drug", "concentration")) %>% 
    group_by(drug) %>% 
    mutate(absolute_abundance = 1.1 * max(absolute_abundance))
  
  p <- d %>% 
    ggplot(aes(concentration_, absolute_abundance, color = run)) + 
    geom_point() + 
    geom_text( data = dl %>% filter(p.value < 0.001), aes(label = p.value_), color="black") +
    facet_wrap(~drug, scales = "free") + 
    ggtitle(d %>% head(1) %>% pull(species)) +
    xlab("concentation (uM)") + ylab("absolute abundance") +
    theme(legend.position = "bottom") 
  
  print(p)

  data.frame()
}

d <- bind_rows(
  dd,
  full_join(
    dd %>% ungroup() %>% select(species, drug) %>% unique(),
    control_counts %>% semi_join(dd, by = c("species")) %>% select(-drug)
  )
) %>% ungroup() %>% arrange(concentration) %>% 
  mutate(concentration_ = fct_inorder(factor(concentration)), 
         species_abbr = stri_replace_all_regex(species, "([A-Z])[^ ]+ (.*)", "$1. $2")) 

dl <- d %>% group_by(species, drug, concentration_) %>% filter(absolute_abundance == max(absolute_abundance)) %>%
  inner_join(median_pvs, by = c("species", "drug", "concentration")) %>%
  group_by(drug, species) %>%
  mutate(absolute_abundance = 1.2 * max(absolute_abundance))

d %>% filter(absolute_abundance > 0) %>% pull(absolute_abundance) %>% min()

det_lims <- detection_limits %>% 
  right_join(d %>% select(run, tech_repl, biol_repl, species_abbr, drug, concentration, concentration_, finalOD) %>% unique(),
             by = join_by(run, biol_repl, tech_repl, drug, concentration))

det_lims <- bind_rows(
  det_lims,
  detection_limits %>% filter(concentration == 0) %>% ungroup() %>% select(-drug) %>% 
    inner_join(d %>% filter(concentration == 0) %>% 
                 select(run, tech_repl, biol_repl, species_abbr, drug, finalOD, concentration_) %>% unique()) 
) %>% 
  mutate(absolute_abundance = detection_limit * finalOD, g = paste(run, biol_repl, tech_repl))

det_lims <- det_lims %>% group_by(species_abbr, drug, concentration_) %>% summarise(absolute_abundance = median(absolute_abundance, na.rm = T))

d %>% 
  ggplot(aes(concentration_, absolute_abundance, color = absolute_abundance+1e-6)) + 
  geom_line(data = det_lims, aes(group=1), linewidth = 0.1, color = "black") +
  geom_point(size = 0.5) + 
  geom_text( data = dl %>% filter(p.value < 0.001), aes(label = p.value_), color="black", size = 2, vjust = 0.75) +
  xlab("Concentation (uM)") + 
  scale_y_continuous(name = "Normalized abundance", breaks = pretty_breaks(3)) +
  scale_color_viridis_c(option = "A", end = 0.8, trans = "log10") +
  facet_wrap(~glue("{species_abbr}\n{drug}"), scales = "free", ncol = 4) + 
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey", color = "lightgrey"), 
        panel.grid.minor = element_blank(), 
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(0,0,0,0),
        legend.position = "bottom"
  ) 

ggsave("suppl_figures/sfig_rare_species.pdf", width = 16, height = 20, units = "cm")

d %>% filter(species_abbr == "B. wadsworthia") %>% 
  ggplot(aes(concentration_, absolute_abundance, color = absolute_abundance+1e-6)) + 
  geom_line(data = det_lims %>% filter(species_abbr == "B. wadsworthia"), aes(group=1), linewidth = 0.1, color = "black") +
  geom_point(size = 0.5) + 
  geom_text( data = dl %>% filter(p.value < 0.001) %>% filter(species_abbr == "B. wadsworthia"), 
             aes(label = p.value_, y = Inf), color="black", size = 2, vjust = 1.25) +
  xlab("Concentation (uM)") + 
  scale_y_continuous(name = "Normalized abundance", breaks = pretty_breaks(3)) +
  scale_color_viridis_c(option = "A", end = 0.8, trans = "log10") +
  ggtitle("B. wadsworthia") +
  facet_wrap(~drug, scales = "free_x", ncol = 5) + 
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey", color = "lightgrey"), 
        panel.grid.minor = element_blank(), 
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(0,0,0,0),
        legend.position = "bottom"
  ) 

ggsave("suppl_figures/B_wadsworthia.pdf", width = 20, height = 10, units = "cm")



pv_filtered <- pvs %>% group_by(drug, species) %>% filter(any(p.value < 0.001))

d %>% semi_join(pv_filtered, by = c("drug", "species")) %>% 
  group_by(species, drug, concentration) %>% summarise(median_abundance = median(absolute_abundance)) %>% 
  left_join(pv_filtered, by = c("species", "drug", "concentration")) %>% 
  write_tsv("suppl_tables/rare_species_pvs_and_abundances.tsv")

