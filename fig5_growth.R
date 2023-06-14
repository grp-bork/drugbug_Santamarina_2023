#!/usr/bin/env Rscript

pacman::p_load(tidyverse) 
pacman::p_load(readxl)
pacman::p_load(stringi)
pacman::p_load(ggtext)
pacman::p_load(tidyverse)

Bv_OE_supernatants <- read_tsv("data/Bv_OE_supernatants.tsv")
Bv_OE_supernatants <- Bv_OE_supernatants %>% filter (
  Supernatant %in% c("B. vulgatus", "B. vulgatus pempty", "B. vulgatus pSG62#2", "B. vulgatus pSG67#2", "B. vulgatus pSG70#2", "MGAM"))

Bv_OE_supernatants <- Bv_OE_supernatants %>% group_by (Sensitive_Strain, Supernatant, Concentration) %>% summarise (medianAUC = median(normAUC))

Bv_OE_supernatants <- Bv_OE_supernatants %>% full_join(tibble(
  Supernatant = c("B. vulgatus pSG62#2","B. vulgatus pSG67#2","B. vulgatus pSG70#2","B. vulgatus pempty","B. vulgatus","MGAM" ),
  supernatent_label = c("\\+ <span style='color:red'><i>Ri</i> C7GA87</span>", "\\+ Pv2039", "\\+ Pv1605", "\\+ empty plasmid", "<i>P. vulgatus</i>", "MGAM")
))

Bv_OE_supernatants$Supernatant <- factor(Bv_OE_supernatants$Supernatant, 
                                            levels = c("B. vulgatus pSG62#2","B. vulgatus pSG67#2","B. vulgatus pSG70#2","B. vulgatus pempty","B. vulgatus","MGAM" ))

Bv_OE_supernatants <- Bv_OE_supernatants %>% ungroup() %>% arrange(Supernatant) %>% mutate(supernatent_label = fct_inorder(supernatent_label))

Bv_OE_supernatants$Concentration<- factor(Bv_OE_supernatants$Concentration, levels = c("0","5", "10"),
                                             labels = c ("0 µM", "5 µM", "10 µM"))

Bv_OE_supernatants <- Bv_OE_supernatants %>% mutate(Sensitive_Strain = stri_replace_all_fixed(Sensitive_Strain, "B. vulgatus", "P. vulgatus"))

Bv_OE_supernatants %>% filter(Concentration %in% c("5 µM", "10 µM"), Sensitive_Strain != "P. copri") %>%
  ggplot(aes(x=medianAUC, y= supernatent_label)) +
  geom_point(aes(shape=Sensitive_Strain), size=1) +
  facet_wrap(.~Concentration) +
  scale_x_continuous(breaks = seq(0,1,0.5), name = "Relative growth (AUC) of sensitive species") +
  scale_y_discrete(name = "Incubated species") +
  scale_shape(name = "Sensitive species") +
  ggtitle("Initial niclosamide concentration") +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    legend.position = "right",
    axis.title = element_text(size = 6), 
    axis.text.x = element_text(size = 5), 
    axis.text.y = element_markdown(size = 5), 
    axis.ticks.length = unit(0, "pt"),
    legend.text = element_text(size = 5, face = "italic"),  
    legend.title = element_text(size = 6), 
    legend.box.margin = margin(0,0,0,-0.4, unit="cm"),
    legend.key.height = unit(0.5, "line"),
    plot.title = element_text(size = 6, hjust = 0.5, margin = margin(0,0,0,0)),
    strip.text = element_text(size = 6)
  )

ggsave("panels/fig5_growth.png", width = 8.2, height = 4, units = "cm", dpi = 1000)
