#!/usr/bin/env Rscript

pacman::p_load(magick)

panel_idx <- 0

next_letter <- function() {
  panel_idx <<- panel_idx + 1
  paste0(" ", LETTERS[panel_idx])
}

add_letter <- function(image) {
  image_annotate(image, next_letter(), size = 1000*8/72, weight = 700)
}

p1 <- image_read("panels/fig3_metab_examples.png") %>% add_letter()
p2 <- image_read("panels/fig3_outcome_protected_2.5.png") %>% add_letter()

whole <- 
  image_append( c(p1, p2), stack = FALSE ) %>% 
  image_background("white")

whole %>% image_write("figures/fig3.png", density = 1000)
