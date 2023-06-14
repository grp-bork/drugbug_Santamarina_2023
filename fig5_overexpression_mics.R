#!/usr/bin/env Rscript

pacman::p_load(tidyverse)
pacman::p_load(ggtext)
pacman::p_load(ggrepel)

MICs <- read_tsv("data/nitroreductase_overexpression.tsv") %>% filter(drug == "Niclosamide")

vectors <-  tribble(~strain, ~vector,
  "5001pww3864#1", "empty plasmid",
  "5001pSG72#1", "Pv0002",
  "5001pSG68#4", "Pv0975",
  "5001pSG69#1", "Pv1605",
  "5001pSG70#1", "Pv2535",
  "5001pSG71#2", "Pv2494",
  "5001pSG73#3", "Pv1973",
  "5001pSG67#1", "Pv2039",
  "5001pSG75#2", "<span style='color:red'><i>Ri</i> C7GB41</span>",
  "5001pSG62#4", "<span style='color:red'><i>Ri</i> C7GA87</span>")

#Figure 5a
p1 <- MICs %>% mutate(vector = factor(vector, levels = rev(vectors$vector))) %>%
  ggplot(aes(MIC90, vector)) + 
  geom_point(size = 0.5) + 
  scale_x_log10(breaks = 2.5*2**(-2:3), label = as.character, name = "Niclosamide IC<sub>90</sub> (µM) in<br>transformed <i>P. vulgatus</i> strain") +
  scale_y_discrete(name = "Overexpressed\nnitroreductase") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_markdown(size = 5, color = "black"), 
    axis.title.x = element_markdown(size = 6, color = "black"),
    axis.title.y = element_text(size = 6), 
    axis.text.x = element_text(size = 5), 
    axis.ticks.length = unit(0, "pt")
)

ggsave("panels/fig5_nitroreductase_overexpression.png", p1, width = 5.5, height = 4, units = "cm", dpi = 1000)

meanAUCs <- read_tsv("data/nitroreductase_overexpression_meanAUCs.tsv")

lbls <- meanAUCs %>% 
  group_by(drug, vector) %>% mutate(d = abs(log(AUC/0.5))) %>% filter(d == min(d))

meanAUCs %>% filter(conc > 0) %>% 
  ggplot(aes(conc, AUC, color = vector)) + geom_point() + geom_line() + 
  geom_label_repel(data = lbls, aes(label = vector), min.segment.length = 0, nudge_x = 0.5, size = 2) +
  facet_grid(~drug) + 
  scale_x_log10(name = "Drug concentration (µM)") +
  ylab("Normalized AUC") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_markdown(size = 5, color = "black"), 
    axis.title.x = element_markdown(size = 6, color = "black"),
    axis.title.y = element_text(size = 6), 
    axis.text.x = element_text(size = 5), 
    axis.ticks.length = unit(0, "pt"),
    legend.position = "none"
  )

ggsave("suppl_figures/nitroreductase_overexpression.png", width = 20, height = 7.5, units = "cm", dpi = 1000)


