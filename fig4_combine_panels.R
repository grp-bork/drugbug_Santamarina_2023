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

p4 <- image_read("panels/fig4_ecoli_metabolomics.png") %>% image_trim() %>% image_border(geometry = "40x40", color = "white")

p1 <- image_read("panels/fig4_niclosamide_classification.png")
p1 <- image_extent(p1, geometry_area(image_info(p4)["width"], image_info(p1)["height"]), gravity = "east", color = "none") %>% add_letter()

p2 <- image_read("illustrations/fig4_schema.png") %>% add_letter()
  
p3 <- image_append( c(image_read("panels/fig4_metabolomics_growth1.png"), image_read("panels/fig4_metabolomics_growth2.png"))) %>% add_letter()

p4 <- p4 %>% add_letter()
p5 <- image_read("illustrations/fig4_schema_ccomes.png") %>% add_letter()
p6 <- image_read("panels/fig4_c_comes_community.png") %>% add_letter()


p1_4 <- image_append( c(p1, p4), stack = TRUE )
p2_3 <- image_append( c(p2, p3), stack = TRUE )
p2_3 <- image_extent(p2_3, geometry_area(image_info(p2_3)["width"], image_info(p1)["height"]), gravity = "None", color = "white")

whole <- 
  image_append(
    c(p1_4,
      image_append( c(p2_3, image_append(c(p5, p6))), stack = TRUE )
    )
  ) %>% 
  image_background("white")

whole %>% image_write("figures/fig4.png", density = 1000)

