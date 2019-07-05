library(tidyverse)
library(ggmap)
library(memoise)
library(lubridate)
library(cowplot)
library(pals)
library(grid)
library(gridExtra)
library(DT)
library(FSA)
library(scales)
library(ggrepel)

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# assign Hawaii isotypes
hi_only_samples <- read.csv(file = "data/fulcrum/hawaii_isotypes.csv") 

#load admixture proportions for K=11 LD8
admix <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K11_Processed_Ancestry.tsv",header = T) %>%
  dplyr::rename(isotype = samples) %>%
  tidyr::gather(pop, frac_pop, - isotype) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(max_pop_frac = max(frac_pop)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(pop, max_pop_frac) %>%
  dplyr::mutate(isotype = factor(isotype)) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(pop_assignment = ifelse(max_pop_frac == frac_pop, pop, NA)) %>%
  dplyr::arrange(isotype, pop_assignment) %>%
  tidyr::fill(pop_assignment) %>%
  tidyr::spread(pop, frac_pop) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE"))

#####################################################################
# Confirm most common haplotype in admix pop B is the global sweep  #
#####################################################################

# Load processed haps
load("data/HAPLOTYPE_LD8/haplotype_plot_df.Rda")

# join admixture info
hap_admix_df <- dplyr::left_join(plot_df, admix)

# Plot % sharing of max haplotype between global admix population D and Hawaiian poulations C, F, G. No filtering of sweep.
admix_sharing_all <- hap_admix_df  %>% 
  distinct(isotype, chromosome, .keep_all= TRUE) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(genome_swept_hap_length = sum(isotype_swept_haplotype_length),
                genome_max_length = sum(max_swept_haplotype_length),
                genome_frac_swept = genome_swept_hap_length/genome_max_length) %>%
  dplyr::ungroup()

admix_sharing <- hap_admix_df  %>%
  dplyr::filter(pop_assignment %in% c("C","E", "A", "H", "F")) %>%
  distinct(isotype, chromosome, .keep_all= TRUE) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(genome_swept_hap_length = sum(isotype_swept_haplotype_length),
                genome_max_length = sum(max_swept_haplotype_length),
                genome_frac_swept = genome_swept_hap_length/genome_max_length) %>%
  dplyr::ungroup()

max_hap_sharing_plot_chrom_all_pops <- ggplot(admix_sharing_all) +
  aes(x = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")), y = max_haplotype_shared, fill = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 1.5, shape = ifelse(admix_sharing_all$Hawaiian == T, 21, 25)) +
  facet_wrap(~chromosome, scales = "free") +
  labs(y = "Fraction most common global haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  geom_text(aes(label=ifelse(max_haplotype_shared > 0.15 & pop_assignment %in% c("E", "A", "H", "F"),
                             as.character(isotype),'')),hjust=1,vjust=.2, size = 2)
max_hap_sharing_plot_chrom_all_pops

ggsave(paste("plots/max_hap_sharing_chr_K=11_LD8.png"), width = 7.5, height = 7.5)
ggsave(paste("plots/max_hap_sharing_chr_K=11_LD8.pdf"), width = 7.5, height = 7.5)

max_hap_sharing_plot_chrom <- ggplot(admix_sharing) +
  aes(x = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")), y = max_haplotype_shared, fill = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 1.5, shape = ifelse(admix_sharing$Hawaiian == T, 21, 25)) +
  facet_wrap(~chromosome, scales = "free") +
  labs(y = "Fraction most common global haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  geom_text_repel(aes(label=ifelse(max_haplotype_shared > 0.3 & pop_assignment != "D" & Hawaiian == T,
                             as.character(isotype),'')),hjust=-2,vjust=.2, size = 2)
max_hap_sharing_plot_chrom

ggsave(paste("plots/max_hap_sharing_chr_K=11_LD8_HI.png"), width = 7.5, height = 7.5)
ggsave(paste("plots/max_hap_sharing_chr_K=11_LD8_HI.pdf"), width = 7.5, height = 7.5)


max_hap_sharing_plot_genome <- ggplot(admix_sharing %>% dplyr::distinct(isotype, .keep_all=T)) +
  aes(x = factor(pop_assignment, levels = c("A", "C", "E", "F", "H")), y = genome_frac_swept, fill = factor(pop_assignment, levels = c("A", "C", "E", "F", "H"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, shape = ifelse(admix_sharing %>% dplyr::distinct(isotype, .keep_all=T) %>% .$Hawaiian == T, 21, 25)) +
  labs(y = "Fraction of genome swept haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  geom_text_repel(aes(label=ifelse(genome_frac_swept > 0.2 & pop_assignment != "C" & Hawaiian == T,as.character(isotype),'')),hjust=-1,vjust=.2, size = 3)
max_hap_sharing_plot_genome

ggsave(paste("plots/max_hap_sharing_genome_K=11_LD8_HI.png"), width = 3.5, height = 5)
ggsave(paste("plots/max_hap_sharing_genome_K=11_LD8_HI.pdf"), width = 3.5, height = 5)


# Analyze differences among populations 
options(scipen=999)
admix_sharing_anlaysis <- admix_sharing %>%
  distinct(isotype, .keep_all= TRUE) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")))

D_test <- dunnTest(data =admix_sharing_anlaysis, genome_frac_swept~pop_assignment, method = "bonferroni") 
D_test[[2]]

# # plot fraction hawaiian for each admixture group
# fraction_admix_population_hawaiian <- admix_sharing_all %>%
#   dplyr::distinct(isotype, pop_assignment, .keep_all = T) %>%
#   dplyr::group_by(pop_assignment) %>%
#   dplyr::mutate(num_isotypes_in_pop = n(),
#                 num_hi_isotypes_in_pop = sum(Hawaiian == "TRUE"),
#                 fraction_hi = num_hi_isotypes_in_pop/num_isotypes_in_pop,
#                 fraction_world = 1-fraction_hi) %>%
#   dplyr::select(pop_assignment, fraction_hi, fraction_world, num_isotypes_in_pop, num_hi_isotypes_in_pop) %>%
#   dplyr::distinct(pop_assignment, .keep_all=T) %>%
#   tidyr::gather(group, value, -pop_assignment, -num_isotypes_in_pop, -num_hi_isotypes_in_pop)
# 
# 
# admix_fractions_pie_chart <- ggplot(fraction_admix_population_hawaiian, aes(x="", y=value, fill=group))+
#   geom_bar(width = 1, stat = "identity") +
#   facet_wrap(~pop_assignment) +
#   coord_polar("y", start=0) +
#   labs(y = "", x = "", fill = "") +
#   theme_minimal()+
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_blank(),
#     panel.grid=element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.x=element_blank()) +
#   geom_text(aes(x=1, y = value, label=ifelse(group == "fraction_hi", num_hi_isotypes_in_pop, num_isotypes_in_pop-num_hi_isotypes_in_pop)))
#   admix_fractions_pie_chart
# 
# #ggsave(paste("plots/admix_hi_fractions.png"), width = 3.75, height = 3.75)
# #ggsave(paste("plots/admix_hi_fractions.pdf"), width = 3.75, height = 3.75)

###### plot hawaii haplotypes and pop D haplotypes by admix population using filtered sweep ratio
mcolor_grp <- plot_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
mcolor <- mcolor_grp$color
names(mcolor) <- mcolor_grp$haplotype

hap_admix_df_ordered <- hap_admix_df %>%
  dplyr::filter(Hawaiian == "TRUE" | pop_assignment == "C") %>%
  dplyr::arrange(desc(pop_assignment), plotpoint) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("C", "E",  "A", "F", "H")))

plotpoints <- hap_admix_df_ordered %>%
  dplyr::distinct(pop_assignment, isotype) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("C", "E",  "A", "F", "H"))) %>%
  dplyr::arrange(pop_assignment) %>%
  dplyr::mutate(plotpoint_hi = row_number())

hap_admix_df_ordered <- dplyr::left_join(hap_admix_df_ordered, plotpoints) #%>%
  # dplyr::mutate(filter = ifelse(plotpoint_hi %in% seq(10, 99, by=1), TRUE, FALSE)) %>%
  # dplyr::filter(filter != TRUE) %>%
  # dplyr::mutate(plotpoint_hi = ifelse(plotpoint_hi > 10, plotpoint_hi-89, plotpoint_hi))
  
mcolor_grp_hi <- hap_admix_df_ordered %>%
  dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "C") %>%
  dplyr::select(haplotype, color) %>% dplyr::distinct()

mcolor_hi <- hap_admix_df_ordered$color
names(mcolor_hi) <- mcolor_grp_hi$haplotype

strain_labels_hi <- hap_admix_df_ordered %>% 
  dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "C") %>%
  dplyr::select(isotype, plotpoint_hi)

ggplot(hap_admix_df_ordered %>% dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "C"),
       aes(xmin = start/1E6, xmax = stop/1E6,
           ymin = plotpoint_hi - 0.5, ymax = plotpoint_hi + 0.5,
           fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = unique(strain_labels_hi$plotpoint_hi),
                     labels = unique(strain_labels_hi$isotype),
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(factor(pop_assignment, levels=c("H", "F", "A", "E", "C"))~chromosome, scales="free", space="free") +
  theme(legend.position="none")

ggsave("plots/haplotype_admix_order_K=11_LD8.png", height = 10, width = 10)
ggsave("plots/haplotype_admix_order_K=11_LD8.pdf", height = 10, width = 10)

# # gray out other haplotypes
# mcolor_grp_hi <- hap_admix_df_ordered %>% 
#   dplyr::filter(Hawaiian == "TRUE") %>%
#   dplyr::select(haplotype, color) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(color = ifelse(color == "#CC0000", color, "grey"))
# 
# mcolor_hi <- mcolor_grp_hi$color
# names(mcolor_hi) <- mcolor_grp_hi$haplotype
# 
# ggplot(hap_admix_df_ordered %>% dplyr::filter(Hawaiian == "TRUE"),
#        aes(xmin = start/1E6, xmax = stop/1E6,
#            ymin = plotpoint_hi - 0.5, ymax = plotpoint_hi + 0.5,
#            fill = haplotype)) +
#   geom_rect() +
#   scale_fill_manual(values = mcolor_hi) +
#   scale_y_continuous(breaks = unique(strain_labels_hi$plotpoint_hi),
#                      labels = unique(strain_labels_hi$isotype),
#                      expand = c(0, 0)) +
#   xlab("Position (Mb)") +
#   theme_bw() +
#   facet_grid(.~chromosome, scales="free", space="free") +
#   theme(legend.position="none")
# 
# test <- hap_admix_df %>%
#   dplyr::filter(pop_assignment == "D") %>%
#   dplyr::distinct(isotype, .keep_all = T)
# 
# ggsave("plots/haplotype_admix_order_grey_K=11.png", height = 5.5, width = 7.5)
# ggsave("plots/haplotype_admix_order_grey_K=11.pdf", height = 5.5, width = 7.5)

### Normal haplotype plots
# # Normal haplotype plot with all strains and haplotype colors
# mcolor_grp <- hap_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
# mcolor <- mcolor_grp$color
# names(mcolor) <- mcolor_grp$haplotype
# 
# strain_labels <- hap_df %>%
#   dplyr::select(isotype, plotpoint)
# 
# ggplot(hap_df,
#        aes(xmin = start/1E6, xmax = stop/1E6,
#            ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
#            fill = haplotype)) +
#   geom_rect() +
#   scale_fill_manual(values = mcolor) +
#   scale_y_continuous(breaks = strain_labels$plotpoint,
#                      labels = strain_labels$isotype,
#                      expand = c(0, 0)) +
#   xlab("Position (Mb)") +
#   theme_bw() +
#   facet_grid(.~chromosome, scales="free", space="free") +
#   theme(legend.position="none")
# 
# #ggsave("haplotype.png", height = 10, width = 10)
# 
# # plot occurance of max haplotype across genome for full set
# strain_labels <- hap_df %>%
#   dplyr::ungroup() %>%
#   dplyr::select(plotpoint, isotype) %>%
#   dplyr::distinct() %>%
#   dplyr::arrange(plotpoint)
# 
# ggplot(hap_df,
#        aes(xmin = start/1E6, xmax = stop/1E6,
#            ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
#            fill = swept_haplotype)) +
#   geom_rect() +
#   scale_fill_manual(values = c("Gray", "Red")) +
#   scale_y_continuous(breaks = strain_labels$plotpoint,
#                      labels = strain_labels$isotype,
#                      expand = c(0, 0)) +
#   xlab("Position (Mb)") +
#   theme_bw() +
#   facet_grid(.~chromosome, scales="free", space="free") +
#   theme(legend.position="none")
# 
# #ggsave(paste("max_haplotype_genome_wide.png"), width = 32, height = 28)

