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

p1 <- image_trim(image_read("panels/fig2_classification_explanation_text.png")) %>% 
  image_extent(geometry_area(2150, 1300), gravity = "southeast", color = "none") %>% 
  image_extent(geometry_area(2200, 1300), gravity = "west", color = "none") %>% 
  add_letter()

p2 <- image_trim(image_read("panels/fig2_classification_explanation_color_bg.png")) %>% 
  image_extent(geometry_area(2150, 1850), gravity = "east", color = "none") %>% 
  image_extent(geometry_area(2200, 1850), gravity = "west", color = "none") %>% 
  add_letter()

p3 <- image_read("panels/fig2_classification_closest_conc_barchart.png") %>% 
  image_trim()  %>% image_border(geometry = "30x30", color = "white") %>% add_letter()

p4 <- image_append(
  c(image_read("panels/fig2_concentration_boxplot_Protection.png") %>% add_letter(),
    image_read("panels/fig2_concentration_boxplot_Sensitization.png") %>% add_letter()
  ), stack = T
) 
    
p1_2 <- image_append( c(p1, p2) , stack = TRUE )

whole <- 
  image_append( c(p1_2, p3, p4) ) %>% 
  image_background("white")

whole %>% image_write("figures/fig2.png", density = 1000)
