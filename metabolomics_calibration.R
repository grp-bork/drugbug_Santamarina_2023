library(tidyverse)

d <- read_tsv("metabolomics_traces/Calibration_curves_data_Exp157.tsv")
d <- d %>% pivot_longer(cols = 3:4, names_to = "compound", values_to = "ion_count")

d <- d %>% mutate(compound = factor(compound, levels = c("Niclosamide", "Aminoniclosamide")))

d %>% ggplot(aes(`Concentration (uM)`, ion_count, color = compound)) + geom_point(size = 0.5) +
  xlab("Compound concentration (µM)") +
  ylab("Ion count (AU)") +
  facet_wrap(~compound, scales = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 6), 
    axis.title.y = element_text(size = 6), 
    axis.text.y = element_markdown(size = 5, vjust = 0.6), # need to fudge to get baseline right
    axis.ticks.length = unit(0, "pt"),
    legend.text = element_text(size = 5),  
    legend.title = element_text(size = 6), 
    strip.text = element_text(size = 6),
    legend.position = "none", 
    panel.grid.minor = element_blank(), 
    plot.title = element_markdown(hjust = 0.5, size = 6, margin = margin(b = 0)),
  )

ggsave("suppl_figures/S5A_metabolomics_calibration.png", width = 18, height = 4, units = "cm", dpi = 1000, bg = "white")


