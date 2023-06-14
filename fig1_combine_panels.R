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

p1 <- image_read("illustrations/fig1_drugbug_community_schema.png") %>% add_letter()

p2 <- image_read("panels/fig1_rel_abundance_barchart.png")
p2 <- image_extent(p2, geometry_area(image_info(p2)["width"]-300, image_info(p2)["height"]), gravity = "center", color = "none") %>% add_letter()

p3 <- image_read("panels/fig1_drug_vs_control_example.png") %>% image_trim() %>% image_border(geometry = "40x40", color = "white") %>% add_letter()
p4 <- image_read("illustrations/fig1_monoculture_schema.png") %>% add_letter()
p5 <- image_read("panels/fig1_single_species_growth_curves.png") %>% image_trim() %>% image_border(geometry = "40x40", color = "white") %>% add_letter()

# add padding around middle panel
w <- image_info(p2)["width"] - image_info(p3)["width"] - image_info(p5)["width"]
p4 <- image_extent(p4, geometry_area(w, image_info(p4)["height"]), gravity = "center", color = "none")

whole <- c(p1, p2, image_append(c(p3, p4, p5))) %>% 
  image_append(stack = T) %>% 
  image_background("white")

whole %>% image_write("figures/fig1.png", density = 1000)
