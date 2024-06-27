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
  scale_y_continuous(name = "Normalized relative\nspecies abundance") + 
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

format_percentage <- function(x) {
  ifelse(x >= 0, sprintf("  +%.0f%%  ", x), sprintf("  \u2212%.0f%%  ", abs(x)))
}

format_fold <- function(x) {
  ifelse(x >= 1, sprintf("  \u2191%g\u00D7  ", signif(x,2)), sprintf("  \u2193%g\u00D7  ", signif(1/x,2)))
}

dd <- d %>% 
  filter(!(species == "C. comes" & Coprococcus == "noCc")) %>% 
  select(Cl, Concentration, species, sl, rel_abundance) %>% 
  pivot_wider(names_from = Concentration, values_from = rel_abundance) %>% 
  mutate(delta = 100*(`10`/`0`-1), 
         delta_label = format_percentage(delta),
         fold = `10`/`0`,
         fold_label = format_fold(fold))

dd

p <- dd %>% 
  ggplot(aes(`0`, `10`, color = sl, shape = sl)) + 
  geom_abline(color = "gray50") +
  geom_point(size = 0.5) +
  geom_text(aes(label = fold_label, 
                hjust = ifelse(`10` > `0` | `0` > 0.75, 1, 0)), 
            color = "black", size = 1.5) +
  coord_fixed(clip = F) +
  scale_color_manual(values=c("red", "black", "gray80", "gray66", "grey33"), name = "Species") +
  scale_shape(name = "Species") +
  scale_x_log10(name = "Species abundance in untreated community", label = as.character) +
  scale_y_log10(name = "Species abundance under\nniclosamide treatment (10 µM)", label = as.character) + 
  facet_grid(~Cl) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.text = element_markdown(size = 5),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        strip.text = element_markdown(size = 6),
        axis.ticks.length = unit(0, "pt")
  )

ggsave("suppl_figures/S6F_c_comes_community.png", p, width = 8, height = 3.5, units = "cm", 
       dpi = 1000, bg = "white")

dd %>% filter(species != "C.comes", delta < 0) %>% 
  group_by(Cl) %>% summarise(mean(delta))


dd %>% filter(species != "C.comes") %>% mutate(abs_log_change = abs(log(`10`/`0`))) %>% 
  group_by(Cl) %>% summarise(abs_log_change = mean(abs_log_change), avg_fold_change = exp(abs_log_change))
