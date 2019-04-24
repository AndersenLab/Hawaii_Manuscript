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

ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# assign Hawaii isotypes
hi_only_samples <- read.csv(file = "data/fulcrum/hawaii_isotypes.csv") 

#load admixture proportions for LD8
admix <- data.table::fread("data/ADMIXTURE_LD8/BEST_K/K7_Processed_Ancestry.tsv",header = T) %>%
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
#load haplotypesfrom 276 strain set
#load("~/Dropbox/AndersenLab/Hawaii_manuscript/data/HAPLOTYPE/haplotype_plot_df.Rda")
# load("data/HAPLOTYPE/processed_haps.Rda")
# 
# # code for coloring haplotypes
# color_plotpoint <- processed_haps[[5]] %>%
#   dplyr::mutate(cvalue = row_number()) %>%
#   dplyr::rename(color = value)
# 
# hap_df <-
#   processed_haps[[3]] %>%
#   dplyr::rename(isotype=haplotype,
#                 haplotype=value) %>%
#   dplyr::group_by(chromosome, isotype) %>%
#   dplyr::arrange(chromosome, isotype, start, stop) %>%
#   # Condense data structure
#   mutate(segment = data.table::rleid(haplotype)) %>%
#   dplyr::group_by(chromosome, isotype, segment) %>%
#   dplyr::mutate(start = min(start), stop = max(stop)) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(cvalue = as.integer(haplotype),
#                 hap_length = stop - start,
#                 haplotype=as.character(haplotype)) %>%
#   dplyr::select(chromosome, start, stop, haplotype, isotype, dplyr::everything()) %>%
#   dplyr::left_join(color_plotpoint, by = c("cvalue")) %>%
#   # Filter empty haplotypes
#   dplyr::filter(!is.na(haplotype)) %>%
#   # Determine the swept haplotype based on max(summation of lengths)
#   dplyr::group_by(chromosome, haplotype) %>%
#   dplyr::mutate(chrom_haplotype_sum = sum(hap_length)) %>%
#   dplyr::group_by(chromosome) %>%
#   dplyr::mutate(swept_haplotype = max(chrom_haplotype_sum) == chrom_haplotype_sum) %>%
#   dplyr::mutate(swept_haplotype_name = ifelse(swept_haplotype, haplotype, NA),
#                 swept_haplotype_name = base::Filter(Negate(is.na), unique(swept_haplotype_name)),
#                 isotype_has_swept_haplotype = sum(swept_haplotype_name == haplotype) > 0,
#                 isotypes_w_haplotype = length(unique(isotype))) %>%
#   # Determine max length of swept haplotype and % shared
#   dplyr::group_by(chromosome, haplotype, isotype) %>%
#   dplyr::mutate(isotype_swept_haplotype_length = sum(ifelse(swept_haplotype, hap_length, 0))) %>%
#   dplyr::group_by(chromosome, isotype) %>%
#   dplyr::mutate(isotype_swept_haplotype_length = max(isotype_swept_haplotype_length)) %>%
#   dplyr::group_by(chromosome) %>%
#   dplyr::mutate(max_swept_haplotype_length = max(isotype_swept_haplotype_length)) %>%
#   dplyr::group_by(chromosome, isotype) %>%
#   dplyr::mutate(max_haplotype_shared = isotype_swept_haplotype_length / max_swept_haplotype_length) %>%
#   dplyr::mutate(filtered_swept_haplotype_len = ifelse(
#     (
#       (hap_length > 1E6)
#       &
#         (max_haplotype_shared > 0.03)
#       &
#         (swept_haplotype == TRUE)
#     ), hap_length, 0)
#   ) %>%
#   dplyr::mutate(filtered_sweep_len = sum(filtered_swept_haplotype_len), 
#                 filtered_sweep_ratio =  (sum(filtered_swept_haplotype_len) / max_swept_haplotype_length),
#                 is_swept = (sum(filtered_swept_haplotype_len) / max_swept_haplotype_length) > 0.03) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(swept_haplotype = ifelse(is_swept == F, F, swept_haplotype)) %>%
#   dplyr::ungroup()

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
  dplyr::filter(pop_assignment %in% c("D","C", "F", "G")) %>%
  distinct(isotype, chromosome, .keep_all= TRUE) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(genome_swept_hap_length = sum(isotype_swept_haplotype_length),
                genome_max_length = sum(max_swept_haplotype_length),
                genome_frac_swept = genome_swept_hap_length/genome_max_length) %>%
  dplyr::ungroup()

max_hap_sharing_plot_chrom_all_pops <- ggplot(admix_sharing_all) +
  aes(x = factor(pop_assignment, levels = c("A", "B", "D", "E", "C", "F", "G")), y = max_haplotype_shared, fill = factor(pop_assignment, levels = c("A", "B", "D", "E", "C", "F", "G"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 1.5, shape = ifelse(admix_sharing_all$Hawaiian == T, 21, 25)) +
  facet_wrap(~chromosome, scales = "free") +
  labs(y = "Fraction most common global haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  geom_text(aes(label=ifelse(max_haplotype_shared > 0.15 & pop_assignment %in% c("C", "F", "G"),
                             as.character(isotype),'')),hjust=1,vjust=.2, size = 2)
max_hap_sharing_plot_chrom_all_pops

ggsave(paste("plots/max_hap_sharing_chr_K=7_LD8.png"), width = 7.5, height = 7.5)
ggsave(paste("plots/max_hap_sharing_chr_K=7_LD8.pdf"), width = 7.5, height = 7.5)

max_hap_sharing_plot_chrom <- ggplot(admix_sharing) +
  aes(x = factor(pop_assignment, levels = c("D", "C", "F", "G")), y = max_haplotype_shared, fill = factor(pop_assignment, levels = c("D", "C", "F", "G"))) +
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

ggsave(paste("plots/max_hap_sharing_chr_K=7_LD8_HI.png"), width = 7.5, height = 7.5)
ggsave(paste("plots/max_hap_sharing_chr_K=7_LD8_HI.pdf"), width = 7.5, height = 7.5)


max_hap_sharing_plot_genome <- ggplot(admix_sharing %>% dplyr::distinct(isotype, .keep_all=T)) +
  aes(x = factor(pop_assignment, levels = c("D", "C", "F", "G")), y = genome_frac_swept, fill = factor(pop_assignment, levels = c("D", "C", "F", "G"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, shape = ifelse(admix_sharing %>% dplyr::distinct(isotype, .keep_all=T) %>% .$Hawaiian == T, 21, 25)) +
  labs(y = "Fraction of genome swept haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  geom_text_repel(aes(label=ifelse(genome_frac_swept > 0.15 & pop_assignment != "D" & Hawaiian == T,as.character(isotype),'')),hjust=-2.75,vjust=.2, size = 3)
max_hap_sharing_plot_genome

ggsave(paste("plots/max_hap_sharing_genome_K=7_LD8_HI.png"), width = 3.75, height = 3.75)
ggsave(paste("plots/max_hap_sharing_genome_K=7_LD8_HI.pdf"), width = 3.75, height = 3.75)


# Analyze differences among populations 
options(scipen=999)
admix_sharing_anlaysis <- admix_sharing %>%
  distinct(isotype, .keep_all= TRUE) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("A", "B", "D", "E", "C", "F", "G")))

D_test <- dunnTest(data =admix_sharing_anlaysis, genome_frac_swept~pop_assignment, method = "bonferroni") 
D_test[[2]]

# plot fraction hawaiian for each admixture group
fraction_admix_population_hawaiian <- admix_sharing_all %>%
  dplyr::distinct(isotype, pop_assignment, .keep_all = T) %>%
  dplyr::group_by(pop_assignment) %>%
  dplyr::mutate(num_isotypes_in_pop = n(),
                num_hi_isotypes_in_pop = sum(Hawaiian == "TRUE"),
                fraction_hi = num_hi_isotypes_in_pop/num_isotypes_in_pop,
                fraction_world = 1-fraction_hi) %>%
  dplyr::select(pop_assignment, fraction_hi, fraction_world, num_isotypes_in_pop, num_hi_isotypes_in_pop) %>%
  dplyr::distinct(pop_assignment, .keep_all=T) %>%
  tidyr::gather(group, value, -pop_assignment, -num_isotypes_in_pop, -num_hi_isotypes_in_pop)


admix_fractions_pie_chart <- ggplot(fraction_admix_population_hawaiian, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  facet_wrap(~pop_assignment) +
  coord_polar("y", start=0) +
  labs(y = "", x = "", fill = "") +
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank()) +
  geom_text(aes(x=1, y = value, label=ifelse(group == "fraction_hi", num_hi_isotypes_in_pop, num_isotypes_in_pop-num_hi_isotypes_in_pop)))
  admix_fractions_pie_chart

ggsave(paste("plots/admix_hi_fractions.png"), width = 3.75, height = 3.75)
ggsave(paste("plots/admix_hi_fractions.pdf"), width = 3.75, height = 3.75)

###### plot hawaii haplotypes and pop B haplotypes by admix population using filtered sweep ratio
mcolor_grp <- hap_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
mcolor <- mcolor_grp$color
names(mcolor) <- mcolor_grp$haplotype

hap_admix_df_ordered <- hap_admix_df %>%
  dplyr::filter(Hawaiian == "TRUE" | pop_assignment == "B") %>%
  dplyr::arrange(desc(pop_assignment), plotpoint) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("B", "F",  "E", "D")))

plotpoints <- hap_admix_df_ordered %>%
  dplyr::distinct(pop_assignment, isotype) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("B", "F",  "E", "D"))) %>%
  dplyr::arrange(pop_assignment) %>%
  dplyr::mutate(plotpoint_hi = row_number())

hap_admix_df_ordered <- dplyr::left_join(hap_admix_df_ordered, plotpoints) %>%
  dplyr::mutate(filter = ifelse(plotpoint_hi %in% seq(10, 105, by=1), TRUE, FALSE)) %>%
  dplyr::filter(filter != TRUE) %>%
  dplyr::mutate(plotpoint_hi = ifelse(plotpoint_hi > 10, plotpoint_hi-95, plotpoint_hi))
  
mcolor_grp_hi <- hap_admix_df_ordered %>% 
  dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "B") %>%
  dplyr::select(haplotype, color) %>% dplyr::distinct()

mcolor_hi <- hap_admix_df_ordered$color
names(mcolor_hi) <- mcolor_grp_hi$haplotype

strain_labels_hi <- hap_admix_df_ordered %>% 
  dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "B") %>%
  dplyr::select(isotype, plotpoint_hi)

ggplot(hap_admix_df_ordered %>% dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "B"),
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
  facet_grid(factor(pop_assignment, levels=c("E","D","F", "B"))~chromosome, scales="free", space="free") +
  theme(legend.position="none")

ggsave("plots/haplotype_admix_order_K=6.png", height = 10, width = 10)
ggsave("plots/haplotype_admix_order_K=6.pdf", height = 10, width = 10)

# gray out other haplotypes
mcolor_grp_hi <- hap_admix_df_ordered %>% 
  dplyr::filter(Hawaiian == "TRUE") %>%
  dplyr::select(haplotype, color) %>%
  dplyr::distinct() %>%
  dplyr::mutate(color = ifelse(color == "#CC0000", color, "grey"))

mcolor_hi <- mcolor_grp_hi$color
names(mcolor_hi) <- mcolor_grp_hi$haplotype

ggplot(hap_admix_df_ordered %>% dplyr::filter(Hawaiian == "TRUE"),
       aes(xmin = start/1E6, xmax = stop/1E6,
           ymin = plotpoint_hi - 0.5, ymax = plotpoint_hi + 0.5,
           fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor_hi) +
  scale_y_continuous(breaks = unique(strain_labels_hi$plotpoint_hi),
                     labels = unique(strain_labels_hi$isotype),
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme(legend.position="none")

ggsave("plots/haplotype_admix_order_grey.png", height = 5.5, width = 7.5)
ggsave("plots/haplotype_admix_order_grey.pdf", height = 5.5, width = 7.5)

### Normal haplotype plots
# Normal haplotype plot with all strains and haplotype colors
mcolor_grp <- hap_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
mcolor <- mcolor_grp$color
names(mcolor) <- mcolor_grp$haplotype

strain_labels <- hap_df %>%
  dplyr::select(isotype, plotpoint)

ggplot(hap_df,
       aes(xmin = start/1E6, xmax = stop/1E6,
           ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
           fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = strain_labels$plotpoint,
                     labels = strain_labels$isotype,
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme(legend.position="none")

#ggsave("haplotype.png", height = 10, width = 10)

# plot occurance of max haplotype across genome for full set
strain_labels <- hap_df %>%
  dplyr::ungroup() %>%
  dplyr::select(plotpoint, isotype) %>%
  dplyr::distinct() %>%
  dplyr::arrange(plotpoint)

ggplot(hap_df,
       aes(xmin = start/1E6, xmax = stop/1E6,
           ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
           fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("Gray", "Red")) +
  scale_y_continuous(breaks = strain_labels$plotpoint,
                     labels = strain_labels$isotype,
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme(legend.position="none")

#ggsave(paste("max_haplotype_genome_wide.png"), width = 32, height = 28)
