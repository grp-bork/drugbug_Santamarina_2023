library(readxl)
library(tidyverse)
library(patchwork)
library(ggtext)

file_list <- read_tsv("metabolomics_traces/files_ecoli.tsv")

# LOAD ----
####Read one file at the time, skip the first section, transform into long format####
df = data.frame()

# to change once I know which ones are the necessary ones
all_files_df <- vector("list", length(file_list$Integrated_Filename))
names(all_files_df) <- file_list$Integrated_Filename

for(x in file_list$Integrated_Filename){
  
  f <- paste(x,".CSV",sep = "")
  
  cat(paste("--------------------------------\n Reading file ",f,"\n",sep=""))
  df1 <- read_lines(file = file.path("metabolomics_traces/",f))
  
  df1_b = df1
  
  # find the breaks in the file
  int_rows = grep("^#.",df1,perl = T)
  
  # every file contains 5 profiles
  df1_b[int_rows]
  
  profiles <- vector("list", length(int_rows)/2)
  names(profiles) <- df1_b[int_rows[seq(1,length(int_rows),2)]]
  
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

all_df <- left_join(all_df,file_list, by = join_by(File == Integrated_Filename))

# Keep factor under control
all_df$Condition <- factor(all_df$Condition, levels = c("mGAM_ctrl","culture_and_supernatant","supernatant"))
all_df$Measured_Mass <- factor(all_df$Measured_Mass, levels = c("297.0197", "324.9782", "330.9991", "247.0638"))

mean_df <- all_df %>% 
  group_by(`X(Minutes)`, Measured_Mass, Condition) %>% 
  summarise(mean_Y_count = mean(`Y(Counts)`)) %>% 
  ungroup() %>% 
  group_by(Measured_Mass, Condition) %>% 
  mutate(rolling_mean_Y_count = (lead(mean_Y_count, n=2)+lead(mean_Y_count)+mean_Y_count+lag(mean_Y_count)+lag(mean_Y_count, n=2))/5) %>% 
  filter(!is.na(rolling_mean_Y_count)) %>% 
  select(`X(Minutes)`, Measured_Mass, Condition, mean_Y_count = rolling_mean_Y_count) %>% 
  group_by(`X(Minutes)`, Measured_Mass, Condition, mean_Y_count)

# # New facet label names for identifying the ions
# supp.labs <- c("Aminoniclosamide", "Niclosamide")
# names(supp.labs) <- c("297.0197", "324.9782")
# 
# all_df %>% filter(`X(Minutes)` < 3.72, `X(Minutes)` > 3.71, Measured_Mass == 247.0638)
# 
# mean_df %>% 
#   filter(Measured_Mass %in% c("297.0197", "324.9782")) %>% 
#   ggplot(aes(x = `X(Minutes)`, y = mean_Y_count, color = Measured_Mass))+
#   geom_line()+
#   theme_light()+
#   scale_color_manual(values = c("297.0197" = "darkred", "324.9782" = "lightblue"),
#                      labels = c("Aminoniclosamide","Niclosamide"))+
#   facet_grid(Measured_Mass ~ Condition, scales = 'free', labeller = labeller(Measured_Mass = supp.labs))+
#   ylab("Ion count [au]") +
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
# ggsave(filename = "overall_traces_horizontal.pdf", height = 6, width = 12, dpi = 300)
# 
# p2 <- mean_df %>% 
#   filter(Measured_Mass == "297.0197", Condition == "mGAM_ctrl", `X(Minutes)`>2.75, `X(Minutes)`<3.50) %>%  
#   ggplot(aes(x = `X(Minutes)`, y = mean_Y_count))+
#   geom_line(color = 'darkred')+
#   # geom_smooth(span = 0.03, se =F, color = 'darkred') +
#   theme_light() +
#   ggtitle("Aminoniclosamide standard in mGAM") +
#   ylab("Ion count [au]") +
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
# ggsave(p2, filename = "Aminoniclosamide_in_mGAM_smoothed_traces.pdf", height = 6, width = 6, dpi = 300)
# 
# p1 <- mean_df %>% 
#   filter(Measured_Mass == "324.9782", Condition == "mGAM_ctrl", `X(Minutes)`>4, `X(Minutes)`<4.75) %>% 
#   ggplot(aes(x = `X(Minutes)`, y = mean_Y_count))+
#   geom_line(color = 'lightblue')+
#   # geom_smooth(span = 0.03, se =F, color = 'lightblue') +
#   theme_light() +
#   ggtitle("Niclosamide standard in mGAM")+
#   ylab("Ion count [au]")+
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
# ggsave(p1, filename = "Niclosamide_in_mGAM_smoothed_traces.pdf", height = 6, width = 6, dpi = 300)
# 
# 
# tmp_p4 <- mean_df %>% 
#   filter(Measured_Mass == "297.0197", Condition %in% c("culture_and_supernatant","supernatant"), `X(Minutes)`>2.75, `X(Minutes)`<3.5) 
# tmp_p4_wide <- tmp_p4 %>% 
#   select(`X(Minutes)`, mean_Y_count, Condition) %>% 
#   unique() %>% 
#   pivot_wider(names_from = Condition, values_from = mean_Y_count) %>% 
#   as.data.frame()
# 
# # create plot object with loess regression lines
# g0 <- ggplot(tmp_p4_wide) + 
#   stat_smooth(aes(x = `X(Minutes)`, y = supernatant), colour = 'darkred', linetype = "dashed", method = "loess", se = FALSE, span = 0.03) +
#   stat_smooth(aes(x = `X(Minutes)`, y = culture_and_supernatant), colour = 'darkred', linetype = "solid", method = "loess", se = FALSE, span = 0.03) +
#   theme_light()
# g0
# 
# # build plot object for rendering 
# gg0 <- ggplot_build(g0)
# 
# # extract data for the loess lines from the 'data' slot
# df1 <- data.frame(x = gg0$data[[1]]$x,
#                   ymin = gg0$data[[1]]$y,
#                   ymax = gg0$data[[2]]$y) 
# 
# # use the loess data to add the 'ribbon' to plot 
# p4 <- g0 +
#   geom_ribbon(data = df1, aes(x = x, ymin = ymin, ymax = ymax),
#               fill = "brown4", alpha = 0.2) + 
#   theme_light() +
#   ggtitle("Aminoniclosamide standard in \nwhole cell culture and supernatant only")+
#   ylab("Ion count [au]")+
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
# ggsave(p4, filename = "Aminoniclosamide_in_SN_vs_mix_smoothed_traces.pdf", height = 6, width = 6, dpi = 300)
# 
# tmp_p3 <- mean_df %>% 
#   filter(Measured_Mass == "324.9782", Condition %in% c("culture_and_supernatant","supernatant"), `X(Minutes)`>4, `X(Minutes)`<4.75)  
# tmp_p3_wide <- tmp_p3 %>% 
#   select(`X(Minutes)`, mean_Y_count, Condition) %>% 
#   unique() %>% 
#   pivot_wider(names_from = Condition, values_from = mean_Y_count) %>% 
#   as.data.frame()
# 
# # create plot object with loess regression lines
# g1 <- ggplot(tmp_p3_wide) + 
#   stat_smooth(aes(x = `X(Minutes)`, y = supernatant), colour = 'lightblue', linetype = "dashed", method = "loess", se = FALSE, span = 0.03) +
#   stat_smooth(aes(x = `X(Minutes)`, y = culture_and_supernatant), colour = 'lightblue', linetype = "solid", method = "loess", se = FALSE, span = 0.03) +
#   theme_light()
# g1
# 
# # build plot object for rendering 
# gg1 <- ggplot_build(g1)
# 
# # extract data for the loess lines from the 'data' slot
# df2 <- data.frame(x = gg1$data[[1]]$x,
#                   ymin = gg1$data[[1]]$y,
#                   ymax = gg1$data[[2]]$y) 
# 
# # use the loess data to add the 'ribbon' to plot 
# p3 <- g1 +
#   geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
#               fill = "lightblue3", alpha = 0.2) + 
#   theme_light() +
#   ggtitle("Niclosamide standard in \nwhole cell culture and supernatant only")+
#   ylab("Ion count [au]")+
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
# ggsave(p3, filename = "Niclosamide_in_SN_vs_mix_smoothed_traces.pdf", height = 6, width = 6, dpi = 300)
# 
# 
# 
# p1 + p2 + p3 + p4 + plot_layout(nrow = 2, guides = 'collect')
# ggsave(filename = "Combined_traces.pdf", height = 10, width = 12, dpi = 300)



d <- tibble(mean_df) %>% 
  full_join(tibble(Measured_Mass = c("297.0197", "324.9782"), 
    Measured_Compound = c("Aminoniclosamide", "Niclosamide"),
    tmin = c(2.75, 4),
    tmax = c(3.75, 5)
  )) %>% ungroup() %>% 
  filter(tmin < `X(Minutes)`, tmax > `X(Minutes)`)

d <- d %>% mutate(Measured_Compound = fct_rev(Measured_Compound))

max_counts <- d %>% group_by(Measured_Compound) %>% summarise(mean_Y_count = max(mean_Y_count))

p1 <- d %>% filter(Condition == "mGAM_ctrl") %>% 
  ggplot(aes(x = `X(Minutes)`, y = mean_Y_count, color = Measured_Compound))+
  geom_point(data = max_counts, x = Inf, color = "white", alpha = 0) +
  geom_line() +
  facet_wrap(. ~ Measured_Compound, scales = 'free') +
  xlab("Retention time (minutes)") +
  ylab("Ion count (AU)") +
  ggtitle("Niclosamide in mGAM, t = 5 h") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5),
    plot.margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"),
    axis.title.x = element_text(size = 6), 
    axis.title.y = element_text(size = 6), 
    axis.text.y = element_markdown(size = 5, vjust = 0.6), # need to fudge to get baseline right
    axis.ticks.length = unit(0, "pt"),
    legend.text = element_text(size = 5),  
    legend.title = element_text(size = 6), 
    strip.text = element_text(size = 6),
    legend.position = "none", 
    panel.grid.minor = element_blank(), 
    plot.title = element_markdown(hjust = 0.5, size = 6, margin = margin(b = 0))
  )

dd <- d %>% filter(Condition != "mGAM_ctrl") %>% 
  mutate(`X(Minutes)` = round(`X(Minutes)`, )) %>%
  group_by(Measured_Compound, Condition, `X(Minutes)`) %>% summarise(mean_Y_count = mean(mean_Y_count)) %>% 
  pivot_wider(names_from = Condition, values_from = mean_Y_count)

dp <- d %>% filter(Condition != "mGAM_ctrl") %>% arrange(Condition, (Condition == "supernatant")*(-1)*`X(Minutes)`)
  
p2 <- d %>% filter(Condition != "mGAM_ctrl") %>% 
  ggplot(aes(x = `X(Minutes)`, y = mean_Y_count, color = Measured_Compound))+
  # geom_ribbon(data = dd, aes(x = `X(Minutes)`, ymin = supernatant, ymax = culture_and_supernatant, fill=Measured_Compound), 
  #             inherit.aes = F, alpha = 0.25) +
  geom_polygon(data = dp, aes(fill=Measured_Compound), alpha = 0.25, linewidth = 0) +
  geom_point(data = max_counts, x = Inf, color = "white", alpha = 0) +
  geom_line(aes(linetype=Condition)) +
  facet_wrap(. ~ Measured_Compound, scales = 'free') +
  xlab("Retention time (minutes)") +
  ylab("Ion count (AU)") +
  ggtitle("Niclosamide in *E. coli*, t = 5 h") +
  scale_color_discrete(guide="none")+
  scale_fill_discrete(labels = c("/", "Bioaccumulated compound"), name = "" )+
  scale_linetype(labels = c("Whole culture", "Supernatant"), name = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5),
    plot.margin = margin(t = 20, r = 10, b = 0, l = 0, unit = "pt"),
    axis.title.x = element_text(size = 6), 
    axis.title.y = element_text(size = 6), 
    axis.text.y = element_markdown(size = 5, vjust = 0.6), # need to fudge to get baseline right
    axis.ticks.length = unit(0, "pt"),
    legend.text = element_text(size = 5),  
    legend.title = element_text(size = 6), 
    strip.text = element_text(size = 6),
    legend.position = "bottom", 
    panel.grid.minor = element_blank(), 
    plot.title = element_markdown(hjust = 0.5, size = 6, margin = margin(b = 0)),
  )

p1/p2+plot_layout(nrow = 2, axis_titles = 'collect')
ggsave("suppl_figures/S5C_metabolomics_ecoli.png", width = 18, height = 8, units = "cm", dpi = 1000, bg = "white")
