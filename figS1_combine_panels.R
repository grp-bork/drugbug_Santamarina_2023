#!/usr/bin/env Rscript

pacman::p_load(magick)

panel_idx <- 0

next_letter <- function() {
  panel_idx <<- panel_idx + 1
  paste0(" ", LETTERS[panel_idx])
}

add_letter <- function(image, extend = FALSE) {
  if (extend) {
    image <- image_extent(image, 
                          geometry_area(image_info(image)["width"], 
                                        image_info(image)["height"]+80), 
                          gravity = "south", 
                          color = "none")
  }
  image_annotate(image, next_letter(), size = 1000*8/72, weight = 700)
}

p1 <- image_read("based_on_external_data/S1A_species_coverage.png") %>% add_letter(extend = TRUE)
p2 <- image_read("suppl_figures/S1B_conditions_colon_concentrations.png") %>% image_trim() %>% add_letter(extend = TRUE)
p3 <- image_read("suppl_figures/S1C_QC_corrs.png") %>% add_letter(extend = TRUE)
p4 <- image_read("suppl_figures/S1D_QC_param_evaluation_protection_sensitization_fraction.png") %>% 
  image_trim() %>% image_border(geometry = "40x10", color = "white") %>% add_letter(extend = FALSE)
p5 <- image_read("suppl_figures/S1E_QC_param_evaluation_pvalue_concentration_dependency.png") %>% 
  image_trim() %>% image_border(geometry = "0x10", color = "white") %>% add_letter(extend = FALSE)
p6 <- image_read("based_on_external_data/S1F_log_control_growth.png") %>% add_letter(extend = TRUE)


p13 <- image_append(c(p1, p3), stack = T)

p123 <- image_append(
  c(
    p13,
    image_extent(p2, geometry_area(image_info(p2)["width"], image_info(p13)["height"]), gravity = "north") 
  )
)

p45 <- image_append(c(p4, p5))
p6 <- image_extent(p6, geometry_area(image_info(p45)["width"], image_info(p6)["height"]), gravity = "west") 

whole <- image_append(c(p123, p45, p6), stack=T) %>% image_background("white")

whole %>% image_write("figures/figS1.png", density = 1000)
