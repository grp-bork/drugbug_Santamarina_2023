#!/usr/bin/env Rscript

pacman::p_load(ragg)
pacman::p_load(egg)
pacman::p_load(tidyverse)
pacman::p_load(patchwork)
pacman::p_load(ggstance)
pacman::p_load(stringi)
pacman::p_load(ggtext)
pacman::p_load(glue)

examples_red <- c("C. comes", "E. coli ED1a", "F. nucleatum", "R. intestinalis")

mb <- read_tsv("data/niclosamide_metabolomics/Metabolomics.tsv")

mbc_ <- mb %>% 
  group_by(TimePoint, Concentration, Strain, MeasuredDrug) %>% summarise( MS = mean(MSdata) ) %>% 
  rename(supernatant = Strain, measured_drug = MeasuredDrug, time = TimePoint, conc = Concentration) %>% 
  mutate(time = as.integer(stri_sub(time, 2))) %>% 
  filter(conc == 10, supernatant != "P. merdae", supernatant != "E. rectale")


gw <- read_tsv("data/niclosamide_metabolomics/AUCs_Supernatants_Matching_Metabolomics.tsv")

new_species_names <- read_tsv("new_species_names.tsv") %>% 
  mutate(
    new_name = stri_replace_all_regex(new_name, "^(.)[^ ]+ ([^ ]+).*", "$1. $2"),
    old_name = stri_replace_all_regex(old_name, "^(.)[^ ]+ ([^ ]+).*", "$1. $2"),
  )

gw <- bind_rows(
  gw %>% anti_join(new_species_names, by = c("strain" = "old_name")),
  gw %>% inner_join(new_species_names, by = c("strain" = "old_name")) %>% select(-strain, strain = new_name)
)

gws <- gw %>% filter(conc > 0) %>% group_by(strain, supernatant, conc, drug, experiment_number) %>% summarise(meanAUC = mean(AUC), sdAUC = sd(AUC) )
gwc <- gws %>% group_by(strain, supernatant, conc, drug) %>% summarise(medianAUC = median(meanAUC) )


mbc <- mbc_ %>% group_by(time, supernatant) %>% mutate(MS_niclo = MS[measured_drug == "niclo"]) %>% ungroup() %>% 
  mutate(supernatant = fct_rev(fct_reorder(supernatant, MS_niclo))) %>% 
  select(-MS_niclo) %>% 
  arrange(supernatant)

mbc <- mbc %>% ungroup() %>% mutate(color = ifelse(supernatant %in% examples_red, "red", "black"),
                                    name = ifelse(supernatant == "MGAM", "MGAM", glue("<i style='color:{color}'>{supernatant}</i>")),
                                    name = stri_replace_all_fixed(name, "ED1a</i>", "</i><span style='color:red'> ED1a</span>"),
                                    name = fct_inorder(name))

mbcl <- mbc %>% 
  semi_join(mbc %>% group_by(supernatant, measured_drug) %>% summarise(MS = mean(MS)) %>% group_by(measured_drug) %>% filter(MS == max(MS)), by=c("supernatant", "measured_drug")) %>% 
  group_by(measured_drug) %>% 
  filter(MS == max(MS)) %>%
  ungroup() %>%
  mutate(label = ifelse(measured_drug == "niclo", "Niclosamide", "Aminoniclosamide"), MS = MS / 2)


p1 <- 
  mbc %>% ggplot(aes(MS, name, fill = measured_drug)) + geom_colh() +
  geom_text(data = mbcl, aes(label = label), size = 5/(14/5), hjust = 0.5, color = "white") +
  scale_x_continuous(breaks = seq(0,10,5), name = paste0("Compound concentration of whole culture (µM)"), expand = c(0,0)) +
  scale_y_discrete(name = "\nIncubated species") +
  scale_fill_manual(values = c("#999999", "#444444")) +
  facet_grid(~paste(time, "h")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 5),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        axis.text.y = element_markdown(size = 5, vjust = 0.6), # need to fudge to get baseline right
        axis.ticks.length = unit(0, "pt"),
        legend.text = element_text(size = 5),  
        legend.title = element_text(size = 6), 
        legend.position = "none",
        strip.text = element_text(size = 6)
  )

ggsave("suppl_figures/suppl_fig_niclosamide_timecourse.pdf", p1, width = 9, height = 3, units = "cm")


mbc <- mbc_ %>% filter(time == 5)

mbc <- mbc %>% semi_join(gw)
mbc <- mbc %>% group_by(supernatant) %>% mutate(MS_niclo = MS[measured_drug == "niclo"]) %>% ungroup() %>% 
  mutate(supernatant = fct_rev(fct_reorder(supernatant, MS_niclo))) %>% 
  select(-MS_niclo) %>% 
  arrange(supernatant)

mbc <- mbc %>% ungroup() %>% mutate(color = ifelse(supernatant %in% examples_red, "red", "black"),
               name = ifelse(supernatant == "MGAM", "MGAM", glue("<i style='color:{color}'>{supernatant}</i>")),
               name = stri_replace_all_fixed(name, "ED1a</i>", "</i><span style='color:red'> ED1a</span>"),
               name = fct_inorder(name))

mbcl <- mbc %>% group_by(measured_drug) %>% 
  filter(MS == max(MS)) %>%
  ungroup() %>%
  mutate(label = ifelse(measured_drug == "niclo", "Niclosamide", "Aminoniclosamide"), MS = MS / 2)

  
p1 <- mbc %>% ggplot(aes(MS, name, fill = measured_drug)) + geom_colh() +
  geom_text(data = mbcl, aes(label = label), size = 5/(14/5), hjust = 0.5, color = "white") +
  scale_x_continuous(breaks = seq(0,10,5), name = paste0("Compound concentration\nof whole culture (µM)"), expand = c(0,0)) +
  scale_y_discrete(name = "\nIncubated species") +
  scale_fill_manual(values = c("#999999", "#444444")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 5),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        axis.text.y = element_markdown(size = 5, vjust = 0.6), # need to fudge to get baseline right
        axis.ticks.length = unit(0, "pt"),
        legend.text = element_text(size = 5),  
        legend.title = element_text(size = 6), 
        legend.position = "none"
  )

p2 <- gwc %>% filter(conc == 5) %>% 
  inner_join(mbc %>% select(supernatant, name) %>% unique(), by = "supernatant") %>% 
  ggplot(aes(medianAUC, name, shape = strain)) + 
  geom_point(size = 1) + 
  scale_x_continuous(breaks = seq(0,1,0.5), name = "Relative growth (AUC)\nof sensitive species") +
  scale_shape(name = "Sensitive species") +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = -10, unit = "pt"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = -5, unit = "pt"),
    legend.position = "right",
    axis.title.x = element_text(size = 6), 
    axis.title.y = element_blank(), 
    axis.text.x = element_text(size = 5), 
    axis.text.y = element_blank(), 
    axis.ticks.length = unit(0, "pt"),
    legend.text = element_text(size = 5, face = "italic"),  
    legend.title = element_text(size = 6), 
    legend.key.height = unit(0.5, "line")
  )

pp1 <- grid.arrange(set_panel_size(p1, width = unit(2.5, "cm"), height = unit(2.5, "cm")))
ggsave("panels/fig4_metabolomics_growth1.png", pp1, width = 4.9, height = 3.5, units = "cm", dpi = 1000)

pp2 <- grid.arrange(set_panel_size(p2, width = unit(2.5, "cm"), height = unit(2.5, "cm")))
ggsave("panels/fig4_metabolomics_growth2.png", pp2, width = 5.5, height = 3.5, units = "cm", dpi = 1000)



AUCs <- read_tsv("data/niclosamide_aminoniclosamide_AUCs.tsv") %>% filter(Sensitive_Strain != "P. merdae", Sensitive_Strain != "E. rectale")
AUCs <- AUCs %>% mutate(conc = as.integer(stri_sub(Concentration, 0, -3)))
AUCs <- AUCs %>% pivot_longer(cols = c("niclosamide", "aminoniclosamide"), names_to = "compound", values_to = "AUC")

AUCs <- AUCs %>% mutate(color = ifelse(Sensitive_Strain %in% examples_red, "red", "black"),
                        Sensitive_Strain = glue("<i style='color:{color}'>{Sensitive_Strain}</i>"))

AUCs %>% filter(conc > 0) %>% mutate(compound = stri_trans_totitle(compound)) %>% 
  ggplot(aes(AUC, Concentration, shape = Sensitive_Strain)) + 
  geom_point() + 
  scale_shape(name = "Species") +
  facet_wrap(~compound, ncol = 1) +
  scale_x_continuous(breaks = seq(0,1,0.5), name = "Relative growth (AUC) of species")+
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 6), 
    axis.text = element_text(size = 5), 
    axis.ticks.length = unit(0, "pt"),
    legend.text = element_markdown(size = 5),  
    legend.title = element_text(size = 6), 
    legend.key.height = unit(0.5, "line"),
    strip.text = element_markdown(size = 5),  
  )
  
ggsave("suppl_figures/suppl_fig_niclosamide_aminoniclosamide_growth.pdf", width = 9, height = 4, units = "cm")





