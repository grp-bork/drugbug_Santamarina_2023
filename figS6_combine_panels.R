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

p1 <- image_read("suppl_figures/S6A_niclosamide_timecourse.png") %>% add_letter()
p2 <- image_read("suppl_figures/S6B_niclosamide_aminoniclosamide_growth.png") %>% add_letter()
p3 <- image_read("suppl_figures/S6C_niclosamide_MIC.png") %>% add_letter()
p4 <- image_read("suppl_figures/S6D_nifurtimox_classification.png") %>% add_letter()
p5 <- image_append(
  c(image_read("suppl_figures/S6E_nifurtimox1.png"), image_read("suppl_figures/S6E_nifurtimox2.png"))
) %>% image_border(color = "white", geometry = "100x0") %>% add_letter()
p6 <- image_read("suppl_figures/S6F_c_comes_community.png") %>% image_border(color = "white", geometry = "100x0") %>% add_letter()
p7 <- image_read("suppl_figures/S6G_nitroreductase_overexpression.png") %>% add_letter()

whole <- image_append(c(image_append(c(p1,p3)), image_append(c(p2,p4)), p5, p6, p7), stack=T)

whole %>% image_write("figures/figS6.png", density = 1000)
