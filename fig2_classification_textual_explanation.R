#!/usr/bin/env Rscript

pacman::p_load(tidyverse, egg)

source("config.R")

d <- tribble(
  ~x, ~y, ~label, ~color,
  "reduced", "normal", "protected in\ncommunity", COLOR_PROTECTION,
  "normal", "reduced", "sensitized in\ncommunity", COLOR_SENSITIZATION,
  "reduced", "reduced", "expected", COLOR_EXP_REDUCED,
  "normal", "normal", "expected", COLOR_EXP_NORMAL
) %>% 
  mutate(x = fct_inorder(x), y = fct_rev(fct_inorder(y)), xy = paste(x, y)) %>% 
  arrange(xy)

p <- d %>% ggplot(aes(x, y)) +
  geom_tile(aes(fill = xy)) +
  geom_text(aes(label = label), size = 2, lineheight = 1) +
  scale_fill_manual(values = d$color) +
  theme_minimal() +
  scale_x_discrete(name = "Growth in monoculture\nupon drug treatment", expand = c(0,0)) +
  scale_y_discrete(name = "Growth in community\nupon drug treatment", expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),        
    axis.text = element_text(size = 5), 
    axis.title = element_text(size = 6),
    NULL
  )

pp <- grid.arrange(set_panel_size(p, width = unit(4, "cm"), height = unit(2, "cm")))

ggsave("panels/fig2_classification_explanation_text.pdf", pp, height = 3.3, width = 6, units = "cm")
ggsave("panels/fig2_classification_explanation_text.png", pp, height = 3.3, width = 6, units = "cm", dpi = 1000)

