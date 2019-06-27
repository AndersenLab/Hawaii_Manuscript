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
library(rcompanion)

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

# Load processed haps
load("data/HAPLOTYPE_LD8/haplotype_plot_df.Rda")

# join admixture info
hap_admix_df <- dplyr::left_join(plot_df, admix)

# Get sweep fractions without filtering of sweep.
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

swept <- admix_sharing_all %>%
  dplyr::filter(!chromosome %in% c("II", "III")) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(swept_chrs = sum(is_swept)) %>%
  dplyr::distinct(isotype, .keep_all = T) %>%
  dplyr::select(isotype, swept_chrs)

# read in climate data (Katie 2018) (these data are for isotype ref strains only right now)
clim <- data.table::fread("~/Hawaii_Manuscript/data/20190606_allstrain_12mo_noaa.csv") %>%
  dplyr::select(-V1) %>%
  dplyr::rename(strain = isotype) %>%
  dplyr::filter(strain %in% admix$isotype)

#look at distribution of distance of station to sample and elevation difference
ggplot(clim %>% dplyr::distinct(strain, .keep_all = T)) +
  aes(x = station_distance) +
  geom_vline(xintercept = 100, color = "red") +
  geom_histogram(binwidth = 50) +
  labs(x = "Station Distance (m)")

# get wild isolate info for all 276 isotypes.
isotype_data <- data.table::fread("~/Hawaii_Manuscript/data/WI_strain_list.csv") %>%
  dplyr::filter(isotype %in% c(as.character(admix$isotype))) %>%
  dplyr::select(isotype, strain, latitude, longitude, isolation_date)
  
# combine admix data with climate data (filter out strains with station distance greter than 100km)
pop_clim <- full_join(admix, clim  %>% dplyr::rename(isotype = strain)) %>% # %>% dplyr::filter(station_distance <= 100)
  tidyr::spread(trait, value) %>% 
  #add in swept chr count
  full_join(., swept)
  

# plot sd for populations 
ggplot(pop_clim) +
  aes(x = factor(pop_assignment, levels = c("A", "B", "D", "E", "C", "F", "G")), y = as.double(sd.temp), fill = pop_assignment) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .25, shape = ifelse(pop_clim$Hawaiian == TRUE, 21, 24)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(y = "12-month temperature SD (°C)", x = "Ancestral population") +
  theme_bw()

# save climate data plot
ggsave('plots/clim_admix_SD_temp.pdf', width = 12, height = 7.5)

ggplot(pop_clim) +
  aes(x = as.character(swept_chrs), y = as.double(sd.temp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = pop_assignment), width = .25, shape = ifelse(pop_clim$Hawaiian == TRUE, 21, 24)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(y = "12-month temperature SD (°C)", x = "Number of swept chromosomes") +
  theme_bw()

# save climate data plot
ggsave('plots/clim_num_chr_swept_SD_temp.pdf', width = 12, height = 7.5)



# # look at climate parameters for admix populations. 
# ggplot(pop_clim) +
#   aes(x = pop_assignment, y = value) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = .25, shape = ifelse(pop_clim$Hawaiian == TRUE, 21, 24)) +
#   facet_wrap(~trait, scales = "free")


