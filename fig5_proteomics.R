#!/usr/bin/env Rscript

pacman::p_load(tidyverse)
pacman::p_load(readxl)
pacman::p_load(ggrepel)
pacman::p_load(ggtext)
pacman::p_load(glue)

d <- bind_rows(
  read_xlsx("data/raw data_Bvulgatus_R.xlsx") %>% 
    select(PROTEIN, starts_with("log10")) %>% 
    pivot_longer(starts_with("log10")) %>% 
    separate(name, into = c(NA, "condition", "concentration")) %>% 
    mutate(species = "P. vulgatus"), 
  
  read_xlsx("data/Roseburia_raw_R_sorted by log 10 untreated.xlsx") %>% 
    select(PROTEIN, starts_with("log10")) %>% 
    pivot_longer(starts_with("log10")) %>% 
    separate(name, into = c(NA, "condition", "concentration")) %>% 
    mutate(species = "R. intestinalis") 
  
)


interesting_genes <-
  tibble(PROTEIN = c("ROSINTL182_06817", "ROSINTL182_07125", "BVU_2494", "BVU_0002", "BVU_1973", "BVU_0975", "BVU_2535", "BVU_1605", "BVU_2039"),
         gene = c("<i>Ri</i> C7GA87", "<i>Ri</i> C7GB41", "Pv2494", "Pv0002", "Pv1973", "Pv0975", "Pv2535", "Pv1605", "Pv2039")
  )

d <- d %>% pivot_wider(names_from = "concentration", names_prefix = "conc_")
d <- d %>% filter(condition != "pSG73")

# use gene name corresponding to condition
d <- d %>% mutate(condition = ifelse(condition == "pSG67", "Pv2039", condition))
d <- d %>% mutate(condition = ifelse(condition == "empty", "empty plasmid", condition))

# no plasmid insert in R. intestinalis
d <- d %>% mutate(condition = ifelse(species == "R. intestinalis", "", glue("<br>+ {condition}")))
d <- d %>% mutate(species_condition = glue("<i>{species}</i>{condition}"))

d <- d %>% filter(!is.na(conc_0uM), !is.na(conc_20uM))

d <- d %>% left_join(interesting_genes)
d <- d %>% mutate(interesting = !is.na(gene))

d <- d %>% arrange(desc(species)) %>% mutate(species_condition = fct_inorder(species_condition))

di <- d %>% filter(interesting)

di <- di %>% mutate(gene_family = "Nitroreductases")

d %>% filter(!interesting) %>% ggplot(aes(conc_0uM, conc_20uM)) + 
  geom_point(size = 0.75, alpha = 0.2, shape = 16) +
  geom_point(data = di, aes(color = gene_family), size = 0.75, shape = 16) +
  geom_richtext(data = di %>% filter(species == "R. intestinalis" | gene == "Pv2039"), aes(label = gene), size = 2, fill = NA, label.color = NA, hjust = c(1,0,1,1)) +
  scale_x_continuous(name = "Protein expression (log10), untreated", limits = c(6.1, 10.9), expand = c(0,0)) +
  scale_y_continuous(name = "Protein expression (log10),\n20 µM niclosamide", limits = c(6.1, 10.9), expand = c(0,0)) +
  facet_wrap(~species_condition) + 
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
      axis.text = element_text(size = 5), 
      axis.title = element_text(size = 6),
      legend.position = c(0,-0.175),
      legend.text = element_text(size = 6), 
      legend.title = element_blank(), 
      legend.key.width = unit(1, "pt"),
      strip.text = element_markdown(size = 6, lineheight = 1.2), 
      plot.title.position = "plot", 
      plot.title = element_text(size = 6, face = "bold")
  ) 
  
ggsave("panels/fig5_proteomics.png", width = 9, height = 4, units = "cm", dpi = 1000)


