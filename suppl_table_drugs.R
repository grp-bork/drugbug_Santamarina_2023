#!/usr/bin/env Rscript

pacman::p_load(tidyverse)
pacman::p_load(stringi)
pacman::p_load(xlsx)

drugs <- read_tsv("data/drug_suppl_data.tsv")
metab_conditions <- read_tsv("data/metabolomics_conditions.tsv")
amplicon_drugs <- read_tsv("data/16S_counts.tsv") %>% select(drug, concentration) %>% unique() %>% filter(concentration > 0)
tested_conditions <- read_tsv("data/Plates_Layout.tsv") %>% rename_all(stri_trans_tolower)  %>% select(drug, concentration) %>% unique() %>% semi_join(amplicon_drugs, by = "drug")
auc_conditions <- read_tsv("data/combined_monoculture_aucs.tsv") %>% select(drug, concentration = conc) %>% unique() %>% semi_join(tested_conditions, by = c("drug", "concentration"))
n_hits <- read_tsv("data/Maier_et_al_hits.tsv")

drugs <- drugs %>% left_join(n_hits)

colon_conc_full <- read_tsv("data/closest_to_gut_concentration.tsv")

colon_conc <- colon_conc_full %>% select(drug, colon_conc_uM) %>% unique() 
colon_conc <- colon_conc %>% mutate(colon_conc_uM = signif(colon_conc_uM, 2))

dt <- drugs %>% 
  left_join(amplicon_drugs %>% group_by(drug) %>% arrange(concentration) %>% summarise(amplicon_conditions = paste(concentration, collapse = ", ")), by = "drug") %>%
  left_join(tested_conditions %>% anti_join(amplicon_drugs) %>% rename(condition_without_growth = concentration) ) %>% 
  left_join(metab_conditions %>% group_by(drug) %>% arrange(concentration) %>% summarise(metab_conditions = paste(concentration, collapse = ", ")), by = "drug") %>% 
  left_join(colon_conc) %>% 
  mutate(metab_conditions = replace_na(metab_conditions, ""), condition_without_growth = replace_na(as.character(condition_without_growth), ""))

wb <- createWorkbook(type="xlsx")

TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)

sheet <- createSheet(wb, sheetName = "Drugs")

column_names <- c("Drug", "Catalogue number", "Supplier", "Solvent", "Prestwick ID", "Target species", "Therapeutic effect", "Number of hits in Maier et al. (20 µM)", "Concentrations (µM) with community and monoculture data",
                  "Concentrations (µM) with no community growth", "Concentrations (µM) with metabolomics data", "Estimated colon concentration")

addDataFrame(dt %>% as.data.frame(), sheet, colnamesStyle = TABLE_COLNAMES_STYLE, row.names = FALSE)

cells <- getCells(getRows(sheet, 1))   # returns all non empty cells
mapply(setCellValue, cells, column_names)

saveWorkbook(wb, "suppl_tables/suppl_table_drugs.xlsx")


conditions <- read_tsv("data/table_sensitization_protection.tsv") %>% select(drug, concentration)

d <- bind_rows(
  conditions %>% mutate(kind = "tested concentration") %>% left_join(colon_conc_full),
  colon_conc_full %>% filter(!is.na(colon_conc_uM)) %>% select(drug, concentration = colon_conc_uM) %>% mutate(kind = "estimated colon concentration", is_closest = FALSE),
) %>% ungroup() %>% mutate(drug = fct_rev(factor(drug))) 

THRESHOLD <- 3

cd <- d %>% filter(kind == "tested concentration") %>% group_by(drug, kind) %>% 
  summarise(c1 = min(concentration)/THRESHOLD, c2 = max(concentration)*THRESHOLD)

dt <- drugs %>% select(drug, target_species) %>% mutate(drug = fct_rev(drug))

d %>% left_join(dt) %>% 
  ggplot(aes(concentration, drug, color = kind)) + 
  geom_segment(data = cd %>% left_join(dt), aes(x = c1, y = drug, xend = c2, yend = drug), size = 2, alpha = 0.25, show.legend = F, lineend = "round") +
  geom_point(aes(size = is_closest)) + 
  geom_text(data = drugs, aes(y = drug, label = paste(N_hit_among_32, "")), x = Inf, hjust = 1, inherit.aes = F, size = 3) +
  scale_x_log10(name = "Concentration (µM)", breaks = 10**(-5:5), label = as.character, expand = c(0,0), limits = c(0.15, 100000-1)) +
  scale_size_manual(breaks = c(TRUE, FALSE), values = c(4,2), labels = c("Shown in Figure 2 (closest to colon concentration, if available)", "other concentration"), name = "") +
  scale_color_discrete(name = "") +
  facet_grid(target_species~., scales = "free_y", space = "free") +
  theme_minimal() +
  ggtitle("", subtitle = "Number of hits\nin Maier et al.") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.background = element_rect(fill = "lightgrey", linewidth = 0),
    plot.subtitle = element_text(hjust = 1, size = 10)
  )

ggsave("suppl_figures/suppl_figure_conditions_colon_concentrations.pdf", width = 20, height = 20, units = "cm")
ggsave("suppl_figures/suppl_figure_conditions_colon_concentrations.png", width = 20, height = 20, units = "cm")

