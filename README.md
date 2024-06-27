# Emergence of community behaviors in the gut microbiota upon drug treatment

Sarela Garcia-Santamarina, Michael Kuhn, Saravanan Devendran, Lisa Maier, Marja Driessen, André Mateus, Eleonora Mastrorilli, Ana Rita Brochado, Mikhail M. Savitski, Kiran R Patil, Michael Zimmermann, Peer Bork, Athanasios Typas

[Preprint](https://www.biorxiv.org/content/10.1101/2023.06.13.544832v1.full)

## Figures and supplementary results

This repository contains all code to generate the main and supplementary figures from the pre-processed data. (Pre-processing includes: determining relative species abundances from 16S sequencing data, determining AUCs from growth curves, processing metabolomics data.) All preprocessed data is in the `data` folder. Statistical analyses and figures are all in R, and can be invoked on the shell using the `Makefile` (`make`).  [`pacman`](https://cran.r-project.org/web/packages/pacman/index.html) is used to load and, if necessary, install missing packages.

`make` will generate the figures in PNG format. `make TIFF` converts the PNGs to TIFFs needed for journal submission, using [ImageMagick](https://imagemagick.org/). 
