#!/usr/bin/env Rscript

pacman::p_load(ragg)
pacman::p_load(egg)
pacman::p_load(tidyverse)
pacman::p_load(patchwork)
pacman::p_load(ggstance)
pacman::p_load(stringi)
pacman::p_load(ggtext)
pacman::p_load(glue)

source("config.R")

mb <- read_tsv("data/Nifur_GM.tsv")

mb <- mb %>% mutate(Strain = stri_replace_all_fixed(Strain, "B. vulgatus", "P. vulgatus"))

mbc <- mb %>% filter(Concentration == 40) %>% 
  group_by(Time, Strain, MeasuredDrug) %>% summarise( MS = mean(`MSdata/IS`) ) %>% 
  rename(supernatant = Strain, measured_drug = MeasuredDrug)

mbc <- mbc %>% ungroup() %>% mutate(MS = MS / mean(MS[measured_drug == "oxidized" & Time == "t0"]))

mbc <- mbc %>% filter(Time != "t0") %>% mutate(Time = paste(stri_sub(Time, 2), "h"))

d_export <- mbc %>% ungroup() %>% rename(incubated_species = supernatant) %>% 
  mutate(measured_drug = ifelse(measured_drug == "oxidized", "Nifurtimox", "Aminonifurtimox"), Time = stri_replace_all_fixed(Time, " ", "")) %>%
  pivot_wider(names_from = Time, values_from = MS, names_prefix = "y_")

addOrUpdateWorksheet("Supplementary Table 3.xlsx", "Nifurtimox metabolomics", d_export, "Data underlying Fig. SE: Incubation of different species with 40 uM nifurtimox")


gw <- read_tsv("data/AUCs_Supernatants_Nifur_GM.tsv")

gw <- gw %>% mutate(strain = stri_replace_all_fixed(strain, "B. vulgatus", "P. vulgatus"))

gw <- gw %>% filter (supernatant %in% c("MGAM", "S. parasanguinis", "E. coli ED1a"), strain %in% c("P. vulgatus", "P. merdae", "C. comes", "B. thetaiotaomicron"))

gws <- gw %>% filter(conc > 0) %>% group_by(strain, supernatant, conc, drug, experiment_number) %>% summarise(meanAUC = mean(AUC), sdAUC = sd(AUC) )
gwc <- gws %>% group_by(strain, supernatant, conc, drug) %>% summarise(medianAUC = median(meanAUC) )

d <- gwc %>% ungroup() %>% 
  filter(conc == 20) %>% select(-conc) %>% 
  inner_join(mbc %>% filter(measured_drug == "reduced"))



mbc <- mbc %>% semi_join(gw)
v <- mbc %>% filter(measured_drug == "reduced") %>% group_by(supernatant) %>% summarise(MS = mean(MS)) %>% arrange(MS) %>% pull(supernatant)
mbc <- mbc %>% mutate(supernatant = factor(supernatant, levels = v)) %>% arrange(supernatant)

xt1 <- 0.0000000
xt2 <- 11
xb1 <- -0.2
xb2 <- 8.2
W <- 0.5*(xb2-xb1)/(xt2-xt1)/2

dt <- d %>% group_by(supernatant) %>% filter(medianAUC == max(medianAUC)) %>% 
  ungroup() %>% 
  mutate(
    supernatant = fct_rev(factor(supernatant, levels = v)),
    x = (as.integer(supernatant)-xt1)*(xb2-xb1)/(xt2-xt1)+xb1 )

dtp <- bind_rows(
  dt %>% mutate(x = MS, y = medianAUC),
  dt %>% mutate(x = x + W, y = Inf),
  dt %>% mutate(x = x - W, y = Inf)
)

examples_red <- c("E. coli ED1a", "S. parasanguinis")

mbc <- mbc %>% ungroup() %>% mutate(color = ifelse(supernatant %in% examples_red, "red", "black"),
                                    name = ifelse(supernatant == "MGAM", "mGAM", glue("<i style='color:{color}'>{supernatant}</i>")),
                                    name = stri_replace_all_fixed(name, "ED1a</i>", "</i><span style='color:red'> ED1a</span>"),
                                    name = fct_inorder(name))

mbc <- mbc %>% mutate(measured_drug = fct_rev(measured_drug))

mbcl <- mbc %>%
  group_by(measured_drug) %>% filter(MS == max(MS)) %>% 
  ungroup() %>% 
  mutate(label = ifelse(measured_drug == "oxidized", "Nifurtimox", "Aminonifurtimox"), MS = MS / 2)


p1 <- mbc %>% ggplot(aes(MS, name, fill = measured_drug)) + geom_colh() +
  geom_text(data = mbcl, aes(label = label), size = 5/(14/5), hjust = 0.5, color = "white") +
  scale_x_continuous(breaks = seq(0,1,0.5), name = paste0("Compound concentration of whole culture (AU)\n "), expand = c(0,0)) +
  scale_y_discrete(name = "\nIncubated species") +
  scale_fill_manual(values = c("#999999", "#444444")) +
  facet_grid(~Time) +
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

p2 <- gwc %>% mutate(Time = "5 h") %>% filter(conc == 20) %>% 
  inner_join(mbc %>% select(supernatant, name) %>% unique(), by = "supernatant") %>% 
  ggplot(aes(medianAUC, name, shape = strain)) + 
  geom_point(size = 1) + 
  scale_x_continuous(breaks = seq(0,1,0.5), name = "Relative growth (AUC)\nof sensitive species") +
  scale_shape(name = "Sensitive species") +
  facet_grid(~Time) +
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
    legend.key.height = unit(0.5, "line"),
    strip.text = element_text(size = 6)
  )

pp1 <- grid.arrange(set_panel_size(p1, width = unit(2.5, "cm"), height = unit(1.2, "cm")))
ggsave("suppl_figures/S6E_nifurtimox1.png", pp1, width = 10, height = 2.5, units = "cm", dpi = 1000, bg="white")

pp2 <- grid.arrange(set_panel_size(p2, width = unit(2.5, "cm"), height = unit(1.2, "cm")))
ggsave("suppl_figures/S6E_nifurtimox2.png", pp2, width = 5.5, height = 2.5, units = "cm", dpi = 1000, bg="white")
