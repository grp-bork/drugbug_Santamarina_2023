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

p1 <- image_read("panels/fig5_nitroreductase_overexpression.png") %>% add_letter()
p2 <- image_read("panels/fig5_proteomics.png") %>% add_letter()
p3 <- image_read("illustrations/fig5_schema.png") %>% add_letter()
p4 <- image_read("panels/fig5_growth.png") %>% add_letter()

p1_2 <- image_append( c(p1, p2))
p3_4 <- image_append( c(p3, p4))
                      
whole <- image_append( c(p1_2, p3_4), stack = TRUE) %>% 
  image_background("white")

whole %>% image_write("figures/fig5.png", density = 1000)
