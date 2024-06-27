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

p1 <- image_read("suppl_figures/S4A_all_metab.png") %>% add_letter()
p2 <- image_read("suppl_figures/S4B_metab_corr_vs_AUC_cutoff_time.png") %>% add_letter()
p3 <- image_append(
  c(image_read("panels/fig3_outcome_protected_10_row.png") %>% add_letter(), 
    image_read("panels/fig3_outcome_sensitized_2.5_row.png") %>% add_letter()),
  stack = T
  ) 

whole <- image_append(
  c(p1, 
    image_append(c(p2, p3))
  ), stack=T) %>% image_background("white")

whole %>% image_write("figures/figS4.png", density = 1000)
