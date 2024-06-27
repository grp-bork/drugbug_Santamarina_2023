#!/usr/bin/env Rscript

pacman::p_load(tidyverse)
pacman::p_load(stringi)
pacman::p_load(glue)

source("config.R")

effect_counts_hit <- read_tsv("data/effect_counts_hit.tsv")

d <- effect_counts_hit %>% 
  group_by(flavor, drug, concentration, category, expected, hit) %>% 
  summarise(fraction = sum(fraction)/n_repl[1], n = sum(n)/n_repl[1], n_total = sum(n_total)/n_repl[1]) %>% 
  filter(!expected) %>% 
  ungroup() 

## PROTECTION

dd <- d %>% filter(flavor == "absolute", !expected, hit) %>% 
  mutate( lc = log(concentration)) %>% 
  group_by(drug) %>% filter(n()>1) %>% 
  mutate(drug = as.factor(as.character(drug)))

nrow(dd)
dd %>% pull(drug) %>% unique() %>% length()

# sigmoid function
b <- nls(fraction ~  1 / (1 + exp( -k* (lc - x0[drug]))), start=list(k=-1, x0= dd %>% group_by(drug) %>% summarise(lc = mean(lc)) %>% pull(lc) ), data=dd)
# baseline: no concentration dependency
b0 <- nls(fraction ~ const[drug], start= list(const = dd %>% group_by(drug) %>% summarise(fraction = mean(fraction)) %>% pull(fraction) ), data=dd)

pv_prot <- anova(b0, b)$`Pr(>F)`[2]

BIC(b0)
BIC(b)

k <- coef(b)[1]
x0 <- tibble(x0 = coef(b)[2:length(coef(b))], drug = levels(dd$drug))

dp <- crossing( drug = unique(dd$drug), concentration = 5*2**seq(-4, 7, 0.1)) %>% mutate( lc = log(concentration))
dp <- dp %>% inner_join(x0) %>% mutate( fraction = 1 / (1 + exp( -k* (lc - x0))) )


d %>% filter(flavor == "absolute", !expected, hit) %>% 
  semi_join(dp, by = "drug") %>% 
  ggplot(aes(concentration, fraction)) +
  geom_line(data = dp, aes(group = drug), color = COLOR_PROTECTION) +
  geom_point(aes(size = n_total), alpha = 0.5) + 
  scale_y_continuous(name = "Fraction of species protected by the community") + 
  scale_x_log10(name = "Drug concentration (µM)") + 
  scale_size_area(name = "Number of species affected in monoculture experiments", max_size = 4) +
  facet_wrap(~drug) + 
  ggtitle(sprintf("Protection in community:\nLogistic function with common growth rate and drug-specific offset\nANOVA p-value compared to constant fraction per drug = %.2g", pv_prot)) +
  theme_minimal() +
  coord_cartesian(clip = F) +
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 7),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
  ) 

ggsave("suppl_figures/S3A_logistic_function_protection.png", width = 15, height = 12, units = "cm", dpi = 1000, bg = "white")

## SENSITIZATION

dd <- d %>% filter(flavor == "absolute", !expected, !hit) %>% 
  mutate( lc = log(concentration)) %>% 
  group_by(drug) %>% filter(n()>1) %>% 
  mutate(drug = as.factor(as.character(drug)))

nrow(dd)
dd %>% pull(drug) %>% unique() %>% length()

# sigmoid function
b <- nls(fraction ~  1 / (1 + exp( -k* (lc - x0[drug]))), start=list(k=1, x0= dd %>% group_by(drug) %>% summarise(lc = mean(lc)) %>% pull(lc) ), data=dd)
# baseline: no concentration dependency
b0 <- nls(fraction ~ const[drug], start= list(const = dd %>% group_by(drug) %>% summarise(fraction = mean(fraction)) %>% pull(fraction) ), data=dd)
  
pv_sens <- anova(b0, b)$`Pr(>F)`[2]

k <- coef(b)[1]
x0 <- tibble(x0 = coef(b)[2:length(coef(b))], drug = levels(dd$drug))

dp <- crossing( drug = unique(dd$drug), concentration = 5*2**seq(-4, 7, 0.1)) %>% mutate( lc = log(concentration))
dp <- dp %>% inner_join(x0) %>% mutate( fraction = 1 / (1 + exp( -k* (lc - x0))) )

d %>% filter(flavor == "absolute", !expected, !hit) %>% 
  semi_join(dp, by = "drug") %>% 
  ggplot(aes(concentration, fraction)) +
  geom_line(data = dp, aes(group = drug), color = COLOR_SENSITIZATION) +
  geom_point(aes(size = n_total), alpha = 0.5) + 
  # geom_text(data = colon_conc %>% , aes(x = X, y = 1, label = marker), color = "black") +
  scale_y_continuous(name = "Fraction of species sensitised by the community") + 
  scale_x_log10(name = "Drug concentration (µM)") + 
  scale_size_area(name = "Number of species not affected in monoculture experiments", max_size = 4) +
  coord_cartesian(clip=F) +
  facet_wrap(~drug) + 
  ggtitle(sprintf("Sensitization in community:\nLogistic function with common growth rate and drug-specific offset\nANOVA p-value compared to constant fraction per drug = %.2g", pv_sens)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 7),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
  ) 

ggsave("suppl_figures/S3B_logistic_function_sensitization.png", width = 15, height = 12, units = "cm", dpi = 1000, bg = "white")






db <- d %>% filter(flavor == "absolute") %>% 
  group_by(drug, category) %>% 
  mutate(i_concentration = factor(conc_names[rank(concentration)], levels = conc_names))
  

db %>% ungroup() %>% summarise(fraction = round(100*median(fraction)), .by = c(hit, i_concentration)) %>% arrange(hit)

fPV <- function(x, a) {
  y <- inner_join( x %>% filter(i_concentration == conc_names[a]), x %>% filter(i_concentration == conc_names[a+1]), by = c("drug") )
  tt <- wilcox.test(y$fraction.x, y$fraction.y, paired = T)
  print(tt)
  tibble(p.value = tt$p.value, category = as.character(x[1,"category"]), x = a+0.5)
}

pvs <- bind_rows(
  fPV(db %>% filter(hit), 1),
  fPV(db %>% filter(hit), 2),
  fPV(db %>% filter(!hit), 1),
  fPV(db %>% filter(!hit), 2)
)


plotBoxplot <- function(CATEGORY, TEXT, COLOR) {
  
  ddb <- db %>% 
    mutate(category = stri_trans_totitle(stri_extract_first_words(category))) %>% 
    filter(category == CATEGORY)

  levels(ddb$i_concentration) <- ddb %>% group_by(i_concentration) %>% summarise(lbl = paste0(as.character(i_concentration[1]), "\n(N=",n(),")")) %>% ungroup() %>%  arrange(i_concentration) %>% pull(lbl)
  
  p <- ddb %>% 
    ggplot(aes(i_concentration, fraction)) + 
    geom_boxplot(varwidth = T, size = 0.5, outlier.size = 0.5, color = COLOR, linewidth = 0.25) + 
    geom_text(data = pvs %>% mutate(category = stri_trans_totitle(stri_extract_first_words(category))) %>% filter(category == CATEGORY), 
              aes(x = x, y = 1, label = paste0("p=", round(p.value, 3))), vjust = 1.5, size = 1.5, color = "black") +
    xlab("Drug concentration step") +
    scale_y_continuous(name = glue("Percentage of {TEXT}\nspecies per drug"), breaks = c(0, 0.5, 1), label=function(v) paste0(v*100, "%")) +  
    theme_minimal() +
    theme(strip.text = element_text(size = 6),
          axis.text = element_text(size = 5), 
          axis.title = element_text(size = 6),
          axis.ticks.length = unit(0, "pt"),
          panel.grid.major.x = element_blank(),
          legend.position = "none"
    )  
  
  ggsave(glue("panels/fig2_concentration_boxplot_{CATEGORY}.png"), p, height = 4, width = 3.5, units = "cm", dpi = 1000)
  ggsave(glue("panels/fig2_concentration_boxplot_{CATEGORY}.pdf"), p, height = 4, width = 3.5, units = "cm")
  
  ddb %>% group_by(i_concentration) %>% summarise(signif(100*median(fraction), 2), n())
}

plotBoxplot("Protection", "protected", COLOR_PROTECTION)
plotBoxplot("Sensitization", "sensitized", COLOR_SENSITIZATION)
 




