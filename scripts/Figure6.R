#!/usr/bin/env Rscript
#Generate plots for Hawaii 2017 collection
library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

args <- commandArgs(trailingOnly = TRUE)

load("data/FST/Genomewide_Fst_complete.rda") 

# Fst
fst_plot <- all_subsample_df_summarized %>%
  ggplot() +
  aes(x = position/1e6, y = Fst)+
  geom_point(alpha = 0.5, color = "black", fill = "#999999", shape = 21)+
  geom_smooth(span = 0.1, se = F, method = "loess")+
  facet_grid(.~CHROM, scales = "free")+
  theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Hudson's Fst")

# Pi
pi_plot <- all_subsample_df_summarized %>%
  dplyr::select(CHROM, position, Diversity_a, Diversity_b) %>%
  tidyr::gather(Pi, value, -CHROM, -position) %>%
  ggplot() +
  aes(x = position/1e6, y = value, 
      color = factor(Pi, levels = c("Diversity_b", "Diversity_a"), labels = c("World", "Hawaii")),
      fill = factor(Pi, levels = c("Diversity_b", "Diversity_a"), labels = c("World", "Hawaii")))+
  geom_point(alpha = 0.5, shape = 21, color = "black")+
  facet_grid(.~CHROM, scales = "free")+
  geom_smooth(span = 0.1, se = F, method = "loess")+
  scale_color_manual(values = c("hotpink3", "cadetblue3"))+
  scale_fill_manual(values = c("hotpink3", "cadetblue3"))+
  theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = expression(pi), color = "Population", fill = "Population")

# Tajima's D
for(chrom in 1:length(grep(pattern = "Neutr", list.files("data/FST/"), value = T))){
  fch <- strsplit(grep(pattern = "Neutr", list.files("data/FST/"), value = T)[chrom], split = "_")[[1]][1]
  if(chrom == 1){
    load(grep(pattern = "Neutr", list.files("data/FST/", full.names = T), value = T)[chrom])
    td_df <- neutrality_df %>%
      dplyr::mutate(CHROM = fch)
  } else {
    load(grep(pattern = "Neutr", list.files("data/FST/", full.names = T), value = T)[chrom])
    td_df <- neutrality_df %>%
      dplyr::mutate(CHROM = fch) %>%
      dplyr::bind_rows(., td_df)
  }
}


td_plot <- td_df %>% 
  na.omit() %>%
  dplyr::distinct(WindowPosition, Population, value, .keep_all = T) %>%
  dplyr::select(CHROM, position = WindowPosition, Population, statistic, value) %>%
  dplyr::filter(statistic == "Tajima.D") %>%
  ggplot() +
  aes(x = position/1e6, y = value, 
      color = factor(Population, levels = c("A", "B"), labels = c("Hawaii", "World")),
      fill = factor(Population, levels = c("A", "B"), labels = c("Hawaii", "World")))+
  geom_point(alpha = 0.5, shape = 21, color = "black")+
  geom_smooth(span = 0.1, se = F, method = "loess")+
  facet_grid(.~CHROM, scales = "free")+
  scale_color_manual(values = c( "cadetblue3","hotpink3"))+
  scale_fill_manual(values = c( "cadetblue3","hotpink3"))+
  theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Tajima's D", color = "Population", fill = "Population")

final_plot <- cowplot::plot_grid(pi_plot, 
                   td_plot,
                   fst_plot,
                   ncol = 1,
                   align = "vh", labels = c("A","B","C"),
                   rel_heights = c(1, 0.8, 0.9))

ggsave(final_plot, filename = "plots/Figure6.pdf", height = 8, width = 12)
ggsave(final_plot, filename = "Plots/Figure6.png", height = 8, width = 12, dpi = 300)

#Calculate geome wide averages and fold change in pi
averages_df <- all_subsample_df_summarized %>%
  dplyr::mutate(genome_avg_pi_word = mean(Diversity_b),
                genome_avg_pi_hi = mean(Diversity_a),
                genome_avg_FST = mean(Fst)) %>%
  dplyr::distinct(genome_avg_pi_word, genome_avg_pi_hi, genome_avg_FST) %>%
  dplyr::mutate(fold_pi = genome_avg_pi_hi/genome_avg_pi_word)
