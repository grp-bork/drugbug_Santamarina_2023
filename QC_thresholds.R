#!/usr/bin/env Rscript

pacman::p_load(glue)
pacman::p_load(stringi)
pacman::p_load(tidyverse)
options(dplyr.summarise.inform = FALSE)

source("config.R")

lbl_AUC <- "Threshold for monoculture experiment (AUC)"
lbl_AR <- "Threshold for community experiment\n(ratio of abundances between treatment and control)"

lbl_absolute <- "Normalized species abundance"
lbl_relative <- "Standardised relative species abundance"

read_tsv("data/combined_monoculture_aucs.tsv") %>%
  rename(concentration = conc) %>% select(-species_name) %>% filter(drug == "MGAM")

effects <- read_tsv("data/effects.tsv") %>% filter(!rare_species, !is.na(AUC))

T_abundance_ratios <- tibble(t_abundance_ratio = seq(0.2, 0.9, 0.05))
T_AUCs <- tibble(t_AUC = seq(0.2, 0.9, 0.05))

d <- crossing(effects, T_abundance_ratios, T_AUCs) %>%
  mutate( reduced_in_community = abundance_ratio < t_abundance_ratio, hit = AUC < t_AUC ) %>%
  group_by(flavor, t_abundance_ratio, t_AUC, hit, reduced_in_community) %>% count() %>%
  left_join(classification, by = c("hit", "reduced_in_community"))

dd <- d %>% group_by(flavor, t_abundance_ratio, t_AUC, category) %>% mutate(fraction = n/sum(n)) %>% filter(!expected) %>% 
  # mutate(category = stri_trans_totitle(category, type = "sentence"))
  mutate(category = stri_replace_all_regex(stri_trans_totitle(category, type = "sentence"), " .*", ""))

dd <- dd %>% mutate(t_abundance_ratio_ = t_abundance_ratio,
                    t_AUC_ = t_AUC,
                    t_abundance_ratio = as.factor(t_abundance_ratio),
                    t_AUC = as.factor(t_AUC))


dd <- dd %>% mutate(
  flavor = ifelse(flavor == "absolute", lbl_absolute, lbl_relative),
  flavor = stri_replace_all_fixed(flavor, " ", "\n")
)

dlbl <- dd %>% filter(t_abundance_ratio_ + t_AUC_ == 1.1) %>%
  mutate(min_perc = floor(10*fraction)*10 ) %>%
  group_by(flavor, category, min_perc) %>%
  filter( n() > 1) %>%
  summarise( x = mean(as.integer(t_abundance_ratio)), y = mean(as.integer(t_AUC)), label = glue(">{min_perc[1]}%") )

dd %>% 
  ggplot(aes(t_abundance_ratio, t_AUC, z = fraction)) + geom_raster(aes(fill = fraction)) + geom_tile(aes(fill = fraction)) +
  geom_point(data = dd %>% filter(t_abundance_ratio_ == T_reduced_in_community, t_AUC_ == T_AUC), color = "white", size = 1) +
  geom_segment(
    data = dd %>% group_by(flavor, category, t_abundance_ratio) %>% filter( floor(10*fraction) != floor(10*lead(fraction)) ),
    aes(x = as.integer(t_abundance_ratio) - 0.5, y = as.integer(t_AUC) + 0.5, xend = as.integer(t_abundance_ratio) + 0.5, yend = as.integer(t_AUC) + 0.5),
    color = "white", lineend = "square"
  ) +
  geom_segment(
    data = dd %>% group_by(flavor, category, t_AUC) %>% filter( floor(10*fraction) != floor(10*lead(fraction)) ),
    aes(x = as.integer(t_abundance_ratio) + 0.5, y = as.integer(t_AUC) - 0.5, xend = as.integer(t_abundance_ratio) + 0.5, yend = as.integer(t_AUC) + 0.5),
    color = "white", lineend = "square"
  ) +
  geom_text(data = dlbl, aes(x = x, y = y, label = label), inherit.aes = F, color = "white", size = 1.5) +
  scale_fill_viridis_c(option = "B", limits = c(0,1), name = "Fraction") +
  scale_x_discrete(name = lbl_AR, breaks = levels(dd$t_abundance_ratio)[seq(1,20,2)]) +
  scale_y_discrete(name = lbl_AUC, breaks = levels(dd$t_AUC)[seq(1,20,2)]) +
  coord_fixed() +
  facet_grid(flavor~category) +
  theme_minimal() +
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0), legend.position = "bottom",
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.75, "line")
  )

ggsave("suppl_figures/S1D_QC_param_evaluation_protection_sensitization_fraction.png", width = 9, height = 9, units = "cm", dpi = 1000)



#### HILL CURVES FOR CONCENTRATION DEPENDENCE
lT_AUC <- 0.75
lT_reduced_in_community <- 0.5
lflavor <- "absolute"

f_pvs <- function(lT_AUC, lT_reduced_in_community, lflavor) {

  message(paste(lT_AUC, lT_reduced_in_community, lflavor))
  
  d <- effects %>%
    filter(flavor == lflavor) %>%
    mutate(hit_def = "AUC", hit = AUC < lT_AUC) %>%
    mutate( reduced_in_community = abundance_ratio < lT_reduced_in_community ) %>%
    inner_join(classification, by = c("hit", "reduced_in_community")) %>%
    ungroup()

  effect_counts_hit <- bind_rows(
    d %>% group_by(flavor, drug, concentration, run, tech_repl, hit, category, result, expected, .drop = F) %>% count() %>%
      semi_join(classification, by = c("hit", "category", "result", "expected")) %>%
      group_by(flavor, drug, concentration, run, tech_repl, hit) %>%
      mutate( n_total = sum(n), fraction = n/n_total ),
  )

  effect_counts_hit <- effect_counts_hit %>% group_by(drug, concentration) %>% mutate(n_repl = length(unique(paste(run, tech_repl))))

  d <- effect_counts_hit %>%
    group_by(flavor, drug, concentration, category, expected, hit) %>%
    summarise(fraction = sum(fraction)/n_repl[1], n = sum(n)/n_repl[1], n_total = sum(n_total)/n_repl[1]) %>%
    filter(!expected) %>%
    ungroup()

  dd <- d %>% filter(!expected, hit) %>% 
    mutate( lc = log(concentration)) %>% 
    group_by(drug) %>% filter(n()>1) %>% 
    mutate(drug = as.factor(as.character(drug)))
  
  # sigmoid function
  for (init_k in c(-1, -0.5, -0.25, -0.1, -0.05)) {
    b <- NULL
    try(
      b <- nls(fraction ~  1 / (1 + exp( -k* (lc - x0[drug]))), start=list(k=init_k, x0= dd %>% group_by(drug) %>% summarise(lc = mean(lc)) %>% pull(lc) ), data=dd),
      silent = TRUE
    )
    if (!is.null(b)) break;
  }
  
  # baseline: no concentration dependency
  b0 <- nls(fraction ~ const[drug], start= list(const = dd %>% group_by(drug) %>% summarise(fraction = mean(fraction)) %>% pull(fraction) ), data=dd)
  
  if (!is.null(b)) {
    pv_prot <- anova(b0, b)$`Pr(>F)`[2]
    k_prot <- coef(b)["k"]
  } else {
    pv_prot <- NA
    k_prot <- NA
  }

  
  # sensitization
  
  dd <- d %>% filter(!expected, !hit) %>%
    mutate( lc = log(concentration)) %>%
    group_by(drug) %>% filter(n()>1) %>%
    mutate(drug = as.factor(as.character(drug)))

  
  # sigmoid function
  for (init_k in c(1, 0.5, 0.25, 0.10, 0.05)) {
    b <- NULL
    try(
      b <- nls(fraction ~  1 / (1 + exp( -k* (lc - x0[drug]))), start=list(k=init_k, x0= dd %>% group_by(drug) %>% summarise(lc = mean(lc)) %>% pull(lc) ), data=dd),
      silent = TRUE
    )
    if (!is.null(b)) break;
  }
  b0 <- nls(fraction ~ const[drug], start= list(const = dd %>% group_by(drug) %>% summarise(fraction = mean(fraction)) %>% pull(fraction) ), data=dd)

  if (!is.null(b)) {
    pv_sens <- anova(b0, b)$`Pr(>F)`[2]
    k_sens <- coef(b)["k"]
  } else {
    pv_sens <- NA
    k_sens <- NA
  }

  tibble(pv_prot, k_prot, pv_sens, k_sens)
}

param_space <- crossing(T_abundance_ratios, T_AUCs, tibble(flavor = c("absolute", "relative"))) %>% 
  arrange(abs(t_abundance_ratio-lT_reduced_in_community)+abs(t_AUC-lT_AUC)) %>% 
  mutate(i = row_number())

pvs <- param_space %>%
  group_by(i, t_abundance_ratio, t_AUC, flavor) %>%
  group_modify(~ f_pvs(.y$t_AUC, .y$t_abundance_ratio, .y$flavor))

## Direction of effect always the same for significant outcomes -- these should not return anything!
pvs %>% filter(pv_prot < 0.01, k_prot > 0)
pvs %>% filter(pv_sens < 0.01, k_sens < 0)

pvs_ <- pvs %>% ungroup() %>% arrange(t_abundance_ratio, t_AUC) %>%
  pivot_longer(5:ncol(pvs), names_pattern = "(.*)_([^_]+)", names_to = c(".value", "kind")) %>%
  mutate(label = ifelse(pv < 0.05, "*", "")) %>%
  mutate(t_abundance_ratio_ = t_abundance_ratio,
         t_AUC_ = t_AUC,
         t_abundance_ratio = as.factor(t_abundance_ratio),
         t_AUC = as.factor(t_AUC))

pvs_ <- pvs_ %>% mutate(
  flavor = ifelse(flavor == "absolute", lbl_absolute, lbl_relative),
  flavor = stri_replace_all_fixed(flavor, " ", "\n"),
  kind = ifelse(kind == "prot", "Protection", "Sensitization")
)

pvs_ %>%
  ggplot(aes(t_abundance_ratio, t_AUC, fill = pv)) + geom_raster() + geom_tile() +
  geom_segment(
    data = pvs_ %>% group_by(flavor, kind, t_abundance_ratio) %>% filter( xor(pv < 0.05, lead(pv < 0.05)) ),
    aes(x = as.integer(t_abundance_ratio) - 0.5, y = as.integer(t_AUC) + 0.5, xend = as.integer(t_abundance_ratio) + 0.5, yend = as.integer(t_AUC) + 0.5),
    color = "white", lineend = "square"
  ) +
  geom_segment(
    data = pvs_ %>% group_by(flavor, kind, t_AUC) %>% filter( xor(pv < 0.05, lead(pv < 0.05)) ),
    aes(x = as.integer(t_abundance_ratio) + 0.5, y = as.integer(t_AUC) - 0.5, xend = as.integer(t_abundance_ratio) + 0.5, yend = as.integer(t_AUC) + 0.5),
    color = "white", lineend = "square"
  ) +
  geom_point(data = pvs_ %>% filter(t_abundance_ratio_ == T_reduced_in_community, t_AUC_ == T_AUC), color = "white", size = 1) +
  coord_fixed() +
  scale_fill_viridis_c(trans = "log10", direction =  -1, name = "p-value of concentration dependency") +
  scale_x_discrete(breaks = T_abundance_ratios$t_abundance_ratio[seq(1, 15, 2)], lbl_AR) +
  scale_y_discrete(breaks = T_AUCs$t_AUC[seq(1, 15, 2)], lbl_AUC) +
  facet_grid(flavor~kind) +
  theme_minimal() +
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0), legend.position = "bottom",
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.75, "line")
  )

ggsave("suppl_figures/S1E_QC_param_evaluation_pvalue_concentration_dependency.png", width = 9, height = 9, units = "cm", dpi = 1000)
