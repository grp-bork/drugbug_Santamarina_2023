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

p1 <- image_read("suppl_figures/S3A_logistic_function_protection.png") %>% add_letter()
p2 <- image_read("suppl_figures/S3B_logistic_function_sensitization.png") %>% add_letter()

whole <- image_append(c(p1, p2), stack=T) %>% image_background("white")

whole %>% image_write("figures/figS3.png", density = 1000)
