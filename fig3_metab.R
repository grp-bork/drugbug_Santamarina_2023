#!/usr/bin/env Rscript

pacman::p_load(ragg)
pacman::p_load(egg)
pacman::p_load(tidyverse)
pacman::p_load(stringi)
pacman::p_load(glue)
pacman::p_load(openxlsx)

source("config.R")

dml <- read_tsv("data/metabolomics_treatment_normalised.tsv")
dmgam <- read_tsv("data/metabolomics_MGAM_normalised.tsv") %>% rename(MGAM_concentration = concentration)

dml <- dml %>% inner_join(DS,  by = "sample_type")

dms <- bind_rows(
  dml, dml %>% filter(time == 0) %>% mutate(sample_type = "SN")
) %>% 
  group_by(drug, concentration, sample_type, time) %>% 
  summarise(y = mean(y)) %>% 
  arrange(sample_type) %>% # gives order: SN, WC
  group_by(drug, concentration, time) %>% 
  mutate(ymax = y, ymin = c(0, y[1])) %>% 
  inner_join(DS,  by = "sample_type")

# pick closest MGAM concentration
dmgam <- dmgam %>% ungroup() %>% select(drug, MGAM_concentration) %>% unique() %>% 
  left_join(dml %>% ungroup() %>% select(drug, concentration) %>% unique()) %>% 
  group_by(drug, MGAM_concentration) %>% 
  mutate(r = abs(log(MGAM_concentration/concentration))) %>% 
  filter(r == min(r)) %>% select(-r) %>% 
  right_join(dmgam)

# correlation between technical replicates
bind_rows(
  dml, dml %>% filter(time == 0) %>% mutate(sample_type = "SN")
) %>% select(drug, concentration, community, sample_type, technical_replicate, time, y) %>% 
  pivot_wider(names_from = "technical_replicate", values_from = "y") %>% 
  group_by(drug, concentration, community, sample_type) %>%
  summarise(corr = cor(tech1, tech2, use = "c")) %>%
  group_by(sample_type) %>% summarise(corr = median(corr))


d_aucs <- bind_rows(
  dml, dml %>% filter(time == 0) %>% mutate(sample_type = "SN")
) %>% arrange(time) %>% 
  select(drug, concentration, community, sample_type, technical_replicate, time, y) %>% 
  group_by(drug, concentration, community, sample_type, technical_replicate) %>% 
  summarise(AUC = sum((time-lag(time))*(y+lag(y))/2, na.rm=T)/max(time))

d_aucs %>% group_by(drug, concentration, sample_type) %>% summarise(sd_AUC = sd(AUC)) %>% 
  ungroup() %>% summarise(sd_AUC = mean(sd_AUC))

# correlation between biological replicates
bind_rows(
  dml, dml %>% filter(time == 0) %>% mutate(sample_type = "SN")
) %>% group_by(drug, concentration, community, sample_type, time) %>% 
  summarise(y = mean(y)) %>% 
  pivot_wider(names_from = "community", values_from = "y") %>% 
  group_by(drug, concentration, sample_type) %>%
  summarise(corr = cor(C1, C2, use = "c")) %>%
  group_by(sample_type) %>% summarise(corr = median(corr))


examples <- tribble(
  ~drug, ~concentration,
  "Ciprofloxacin", 20,
  "Niclosamide", 20,
  "Mefloquine", 20,
  "Lansoprazole", 80
)

dmgam <- dmgam %>% mutate(sample_type = "MGAM") %>% inner_join(DS) %>% ungroup() %>% select(-MGAM_concentration)

dmlc <- dml %>% 
  group_by(time, drug, concentration, sample_info) %>% summarise(y_sd = sd(y), y = mean(y)) 

dmlc <- bind_rows(
  dmlc,
  dmlc %>% filter(time == 0) %>% mutate(sample_type = "SN") %>% select(-sample_info) %>% inner_join(DS,  by = "sample_type"),
)


d_export <- dmlc %>% select(-sample_type) %>% arrange(drug, concentration, time, sample_info) 

d_export <- bind_rows(
  d_export,
  d_export %>% group_by(time, drug, concentration) %>% 
    summarise(y = y[2]-y[1], y_sd = sqrt(y_sd[1]**2 + y_sd[2]**2)*(time[1] != 0) ) %>% 
    mutate(sample_info = "Bioaccumulation (inferred)")
) %>% 
  pivot_wider(names_from = time, values_from = c("y", "y_sd")) %>% 
  arrange(drug, concentration)

addOrUpdateWorksheet("Supplementary Table 3.xlsx", "Community metabolomics", d_export, "Data underlying Fig. 3: Community metabolomics")

dms <- bind_rows(
  dml, 
  dml %>% filter(time == 0) %>% mutate(sample_type = "SN") %>% select(-sample_info) %>% inner_join(DS,  by = "sample_type"),
  dmgam %>% select(-sample_info) %>% mutate(sample_type = "MGAM") %>% inner_join(DS,  by = "sample_type")
) %>% 
  ungroup() %>% 
  mutate(sample_type = fct_rev(fct_reorder(sample_type, y))) %>% 
  arrange(sample_type) %>% 
  group_by(drug, concentration, sample_type, sample_info, time) %>% 
  summarise(y = mean(y)) %>% 
  group_by(drug, concentration, time) %>% 
  mutate(ymax = y, ymin = lead(y, default = 0))


dx <- dms %>% select(-ymin, -ymax, -sample_info) %>% pivot_wider(names_from = "sample_type", values_from = "y")

dr <- bind_rows(
  dx %>% mutate(kind = "Biotransformation") %>% select(kind, drug, concentration, time, ymin = WC, ymax = MGAM, -SN),
  dx %>% mutate(kind = "Bioaccumulation") %>% select(kind, drug, concentration, time, ymin = SN, ymax = WC, -MGAM),
)

dmlc_ <- dmlc %>% semi_join(examples, by = c("drug", "concentration"))
dml_ <- dml %>% semi_join(examples, by = c("drug", "concentration"))
dms_ <- dms %>% semi_join(examples, by = c("drug", "concentration"))
dr_ <- dr %>% semi_join(examples, by = c("drug", "concentration"))
dmgam_ <- dmgam %>% semi_join(examples, by = c("drug"))


p <- dmlc_ %>% 
  ggplot(aes(time, y, color = sample_info)) + 
  geom_ribbon(data = dr_, aes(time, ymin = ymin, ymax = ymax, fill = kind), inherit.aes = F, alpha = 0.5) +
  geom_line(data = dmgam_, aes(time, y, color = sample_info), inherit.aes = F) +
  geom_errorbar(data = dmgam_, aes(time, ymin = y-y_sd, ymax = y+y_sd, color = sample_info), width = 0.25, inherit.aes = F) +
  geom_line() +
  geom_errorbar(aes(ymin = y-y_sd, ymax = y + y_sd), width = 0.25) + 
  scale_x_continuous(name = "Time (h)", labels = as.character) +
  scale_y_continuous(name = "Normalized MS value", breaks = c(0,1)) +
  scale_color_manual(name = "Measurement", values = c("#2171b5", "#6baed6","#bdd7e7"), breaks = dms_ %>% pull(sample_info) %>% unique(), 
                     labels = dms_ %>% mutate(sample_info = stri_replace_all_regex(sample_info, "\\(.*", "")) %>% pull(sample_info) %>% unique()) +
  scale_fill_manual(name = "Effect of community", values = c("grey","black"), breaks = c("Biotransformation", "Bioaccumulation")) +
  facet_wrap(~glue("{drug} ({concentration} µM)"), ncol = 1) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey", color = "lightgrey"), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 7),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        # legend.position = "none"
        legend.position = "bottom",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(-5,0,0,0),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.box.just = "left",
        axis.ticks.length = unit(0, "pt") 
  )


pp <- grid.arrange(set_panel_size(p, width = unit(3.2, "cm"), height = unit(1.5, "cm")))

ggsave("panels/fig3_metab_examples.png", plot = pp, width = 4, height = 12, units = "cm", dpi = 1000)
ggsave("panels/fig3_metab_examples.pdf", plot = pp, width = 6, height = 12, units = "cm")


p <- dmlc %>% 
  ggplot(aes(time, y, color = sample_info)) + 
  geom_ribbon(data = dr, aes(time, ymin = ymin, ymax = ymax, fill = kind), inherit.aes = F, alpha = 0.5) +
  geom_line(data = dmgam, aes(time, y, color = sample_info), inherit.aes = F) +
  geom_errorbar(data = dmgam, aes(time, ymin = y-y_sd, ymax = y+y_sd, color = sample_info), width = 0.25, inherit.aes = F) +
  geom_line() +
  geom_errorbar(aes(ymin = y-y_sd, ymax = y + y_sd), width = 0.25) + 
  scale_x_continuous(name = "Time (h)", labels = as.character) +
  scale_y_continuous(name = "Normalized MS value", breaks = c(0,1)) +
  scale_color_manual(name = "Measurement", values = c("#2171b5", "#6baed6","#bdd7e7"), breaks = dms %>% pull(sample_info) %>% unique(), 
                     labels = dms %>% mutate(sample_info = stri_replace_all_regex(sample_info, "\\(.*", "")) %>% pull(sample_info) %>% unique()) +
  scale_fill_manual(name = "Effect of community", values = c("grey","black"), breaks = c("Biotransformation", "Bioaccumulation")) +
  facet_wrap(~glue("{drug}\n({concentration} µM)"), ncol = 8) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey", color = "lightgrey"), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 7),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.position = "inside",
        legend.position.inside = c(0.81, 0.1),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(-5,0,0,0),
        legend.box = "horizontal",
        legend.direction = "vertical",
        legend.box.just = "left",
        axis.ticks.length = unit(0, "pt") 
  )

ggsave("suppl_figures/S4A_all_metab.png", plot = p, width = 18, height = 12, units = "cm", dpi = 1000, bg = "white")




L <- c("Biotransformation only", "Bioaccumulation only", "Biotransformation & bioaccumulation")

dcorr <- read_tsv("data/metab_effect_correlation.tsv") %>% 
  mutate(sample_type = factor(stri_replace_all_fixed(sample_type, "_", " "), levels = L))

de <- read_tsv("data/metab_effects.tsv") %>% 
  mutate(sample_type = factor(stri_replace_all_fixed(sample_type, "_", " "), levels = L))

dcorr %>% pivot_longer(starts_with("pvalue"), names_to = "correlation_type", values_to = "pvalue") %>%
  ggplot(aes(max_time, pvalue, color = sample_type)) + geom_line() + geom_point() +
  geom_hline(yintercept = 0.01) +
  geom_hline(yintercept = 0.05) +
  scale_y_log10(breaks = 0.05*5**(-2:2)) +
  facet_grid(correlation_type ~ result, scales = "free", space = "free") +
  theme(legend.position = "bottom")

asterisks <- function(x) ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", ""))

dcorr <- dcorr %>% mutate(
  label = asterisks(pvalue_s)
)

dcorr <- dcorr %>% arrange(-corr_s) %>% group_by(result, max_time) %>% 
  mutate(nudge = -1*(row_number() > 1 & lag(corr_s) - corr_s < 0.05))

p <- dcorr %>% 
  ggplot(aes(max_time, corr_s, color = sample_type)) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_text(aes(label = label), nudge_y = 0.02+0.09*dcorr$nudge, show.legend = F) +
  scale_y_continuous(name = "Spearman correlation between\nAUC and fraction of affected species") +
  scale_x_continuous(limits = c(0,10), labels = as.character, name = "Time cutoff for AUC calculation (h)") +
  scale_size_area(trans = "log10", breaks = c(0.002, 0.01, 0.05), max_size = 2, labels = as.character) +
  scale_color_grey(start = 0.8, end = 0, name = "") +
  facet_wrap(~ stri_trans_totitle(result, type='sentence'), ncol = 1) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 7),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        strip.background = element_rect(color = "lightgrey", fill = "lightgrey"),
        legend.margin = margin(0,0,0,0),
        legend.key.height = unit(0.5, "line"),
        legend.position = "bottom", 
        legend.direction = "vertical"
  ) 

ggsave("suppl_figures/S4B_metab_corr_vs_AUC_cutoff_time.png", plot = p, width = 6, height = 10, units = "cm", dpi = 1000, bg="white")


L <- "Other drugs"
de <- de %>% ungroup() %>% mutate(Drug = fct_other(drug, keep = examples %>% pull(drug), other_level = L))

plotOutcomes <- function(result_, max_time_) {
  de_ <- de %>% filter(result == result_, max_time == max_time_)
  p <- de_ %>% 
    ggplot(aes(AUC, fraction)) + 
    geom_point(data = de_ %>% anti_join(examples, by = c("drug", "concentration")), aes(color = Drug, size = n_total)) +
    geom_point(data = de_ %>% semi_join(examples, by = c("drug", "concentration")), aes(color = Drug, size = n_total)) +
    geom_text(
      data = dcorr %>% filter(result == result_, max_time == max_time_) %>% mutate(y = ifelse(result == "protected in community", 0.05, -0.05)),
      aes(x = 0, y = y, label = paste("r[s] == ", round(corr_s, 2),
                                      "~(p == ", sprintf("%.2g", pvalue_s), ")")),
      parse = TRUE, size = 2,
      vjust = 0, hjust = 0, color = "black"
    ) +
    scale_size_area(max_size = 1.5, name = "Number of affected species") +
    scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "black"), breaks = levels(de$Drug)) +
    xlab(glue("Decrease in drug concentration after {max_time_} h (AUC)")) +
    guides(colour = guide_legend(order = 1), 
           size = guide_legend(order = 2, override.aes = list(shape = 21, stroke = 0.25))) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 7),
          strip.text = element_text(size = 6),
          axis.text = element_text(size = 5), 
          axis.title = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          strip.background = element_rect(color = "lightgrey", fill = "lightgrey"),
          legend.margin = margin(0,0,0,-10),
          legend.key.height = unit(0.2, "line"),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.ticks.length = unit(0, "pt") 
    ) 
  
  pp <- p + facet_wrap(~sample_type, ncol = 1) + ylab(paste0("Fraction of species that are ", result_))

  ggsave(paste0("panels/fig3_outcome_", stri_split_fixed(result_, " ")[[1]][1], "_", max_time_, ".pdf"), pp, width = 5.5, height = 12, units = "cm")
  ggsave(paste0("panels/fig3_outcome_", stri_split_fixed(result_, " ")[[1]][1], "_", max_time_, ".png"), pp, width = 5.5, height = 12, units = "cm", dpi = 1000, bg = "white")

  pp <- p + facet_wrap(~sample_type, nrow = 1) + ylab(paste0("Fraction of species that are\n", result_))
  
  ggsave(paste0("panels/fig3_outcome_", stri_split_fixed(result_, " ")[[1]][1], "_", max_time_, "_row.pdf"), pp, width = 12, height = 6, units = "cm")
  ggsave(paste0("panels/fig3_outcome_", stri_split_fixed(result_, " ")[[1]][1], "_", max_time_, "_row.png"), pp, width = 12, height = 6, units = "cm", dpi = 1000, bg = "white")
  
}

plotOutcomes("protected in community", 1.25)
plotOutcomes("protected in community", 2.5)
plotOutcomes("protected in community", 10)
plotOutcomes("sensitized in community", 1.25)
plotOutcomes("sensitized in community", 2.5)
plotOutcomes("sensitized in community", 10)
