library(readxl)
library(tidyverse)
library(glue)
library(ggtext)
library(stringi)
library(patchwork)

file_list <- read_tsv("metabolomics_traces/files.tsv")

#### Read one file at the time, skip the first section, transform into long format ####
df = data.frame()

all_files_df <- vector("list", length(file_list$File))
names(all_files_df) <- file_list$File

for(x in file_list$File){
  
  f <- gsub(".d",".CSV",x)
  
  cat(paste("--------------------------------\n Reading file ",f,"\n",sep=""))
  df1 <- read_lines(file = file.path("metabolomics_traces/",f))
  
  df1_b = df1
  
  # find the breaks in the file
  int_rows = grep("^#.",df1,perl = T)
  
  # every file contains 3 profiles
  df1_b[int_rows]
  
  profiles <- vector("list", length(int_rows)/2)
  names(profiles) <- df1_b[int_rows[c(1,3,5)]]
  
  int_rows <- c(int_rows,length(df1_b)+1) # add artificial end for loop later
  
  for(i in 1:length(profiles)){
    profiles[[i]] <- sapply(df1_b[(int_rows[(i*2)] + 1):(int_rows[((i*2)+1)]-1)],
                            function(x) as.numeric(unlist(strsplit(x, ',',fixed = T))))
    profiles[[i]] <- as.data.frame(t(profiles[[i]]))
    names(profiles[[i]]) <- unlist(strsplit(df1_b[int_rows[2]], ',',fixed = T))
  }
  
  all_traces <- data.table::rbindlist(profiles,idcol = T)
  
  # adjust
  file_df <- all_traces %>% 
    filter(!str_detect(`.id`, pattern = "TIC")) %>% #we don't need the TIC
    # we need to extract the measured mass
    mutate(Measured_Mass = gsub('#"+/-ESI EIC(',
                                "",
                                gsub("\\) Scan Frag=.*$","",`.id`,perl = T), 
                                fixed = T)) 
  
  all_files_df[[x]] <- file_df
}

all_df <- data.table::rbindlist(all_files_df, idcol = 'File')

all_df <- left_join(all_df,file_list)

# I need to reorder factors if I want MGAM on top
all_df$Strain <- factor(all_df$Strain, levels = c("MGAM","C_comes","R_intestinalis"))
levels(all_df$Strain) <- c("mGAM", "C. comes", "R. intestinalis")

# I need to reorder factors if I want MGAM on top
all_df$Measured_Mass <- factor(all_df$Measured_Mass, levels = c("324.9782","297.0197"))

compound_names <- c(
  '324.9782'="Niclosamide",
  '297.0197'="Aminoniclosamide"
)

all_df$Compound <- factor(all_df$Compound, levels = c("Niclo", "Amino"))
levels(all_df$Compound) <- c("Niclosamide", "Aminoniclosamide")

d <- tibble(all_df) %>% filter(Time == "t0", Compound == "Niclosamide", Strain == "MGAM")

max_counts <- all_df %>% group_by(Measured_Mass) %>% summarise(`Y(Counts)` = max(`Y(Counts)`))

plotPanel <- function(Strain, Compound, Time, d) {
  Time <- stri_replace_all_regex(Time, "t(.)", "t = $1 h")
  if (Strain != "mGAM") {
    Strain <- glue("*{Strain}*")
  }
  title <- glue("{Compound} in {Strain}, {Time}")
  d <- d[[1]]
  p <- d %>% ggplot(aes(x = `X(Minutes)`, y = `Y(Counts)`, color = Measured_Mass))+
    geom_point(data = max_counts, x = Inf, color = "white", alpha = 0) +
    geom_line() +
    facet_wrap(. ~ Measured_Mass, scales = 'free', labeller = labeller(Measured_Mass=as_labeller(compound_names))) +
    xlab("Retention time (minutes)") +
    ylab("Ion count (AU)") +
    ggtitle(title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 5),
      plot.margin = margin(t = 10, r = 10, b = 0, l = 0, unit = "pt"),
      axis.title.x = element_text(size = 6), 
      axis.title.y = element_text(size = 6), 
      axis.text.y = element_markdown(size = 5, vjust = 0.6), # need to fudge to get baseline right
      axis.ticks.length = unit(0, "pt"),
      legend.text = element_text(size = 5),  
      legend.title = element_text(size = 6), 
      strip.text = element_text(size = 6),
      legend.position = "none", 
      panel.grid.minor = element_blank(), 
      plot.title = element_markdown(hjust = 0.5, size = 6, margin = margin(b = 0)),
    )
  
  list(p)
}

plots <- all_df %>% group_by(Strain, Compound, Time) %>% 
  nest() %>% 
  summarise(p = plotPanel(Strain, Compound, Time, data))

# plots[1,]$p

wrap_plots(plots$p, ncol = 2) + plot_layout(axis_titles = 'collect')

ggsave("suppl_figures/S5B_metabolomics.png", width = 18, height = 13, units = "cm", dpi = 1000, bg = "white")

