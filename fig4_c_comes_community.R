#!/usr/bin/env Rscript

pacman::p_load(tidyverse)
pacman::p_load(stringi)
pacman::p_load(glue)
pacman::p_load(ggtext)


d <- read_tsv("data/c_comes_community.tsv")

d <- d %>% mutate(species = stri_replace_all_fixed(Condition, "_pellet", ""))

d <- d %>% 
  mutate(species = fct_relevel(species, "C. comes")) %>% arrange(species) %>% 
  mutate(
    Cl = ifelse(Coprococcus == "Cc", "<span  style='color:red'>with <i>C. comes</i></span>", "without <i>C. comes</i>"),
    sl = fct_inorder(glue("<i>{species}</i>"))
  ) 
d <- d %>% mutate(Cl = factor(Cl, level = rev(d$Cl %>% unique())))


p <- d %>%
  ggplot(aes(as.character(Concentration), rel_abundance, fill = sl)) + geom_col(width = 0.6) +
  scale_fill_manual(values=c("red", "black", "gray80", "gray66", "grey33"), name = "Species") +
  scale_x_discrete(name = "Niclosamide concentration (µM)") +
  scale_y_continuous(name = "Normalised relative\nspecies abundance") + 
  facet_grid(~Cl) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.text = element_markdown(size = 5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        strip.text = element_markdown(size = 6),
        axis.ticks.length = unit(0, "pt")
  )

ggsave("panels/fig4_c_comes_community.png", p, width = 7.5, height = 3.5, units = "cm", dpi = 1000)

