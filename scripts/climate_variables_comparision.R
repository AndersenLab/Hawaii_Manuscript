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

# read in climate data (Katie manuscript 2017)
clim <- data.table::fread("~/Hawaii_Manuscript/data/strain_climate_data.csv")

# combine sweep data with climate data
sweep_clim <- full_join(admix_sharing_all, clim) %>%
  dplyr::filter(isotype %in% clim$isotype) %>%
  tidyr::gather(env_trait, value, ceil_hgt_min:ws_var) %>%
  dplyr::mutate(env_trait = paste0(time_period, "_", env_trait))
  
# plot climate variables vs fraction sweep
test <- ggplot(sweep_clim %>% dplyr::filter(chromosome == "V")) +
  aes(x = value, y = max_haplotype_shared, fill = factor(pop_assignment, levels = c("A", "B", "D", "E", "C", "F", "G"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_point(shape = 21) +
  #geom_jitter(width = 0.25, size = 1.5, shape = ifelse(admix_sharing_all$Hawaiian == T, 21, 25)) +
  facet_wrap(~env_trait, scales = "free") +
  #labs(y = "Fraction most common global haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none")



temp_traits <- c("3yr_temp_min", "3mo_temp_min", "1yr_temp_min","3yr_temp_avg",
                 "3mo_temp_avg","1yr_temp_avg", "3yr_temp_max", "3mo_temp_max",
                 "1yr_temp_max", "3yr_temp_var", "3mo_temp_var", "1yr_temp_var")
temp_sweep_clim <- sweep_clim %>%
  dplyr::filter(env_trait %in% temp_traits) %>%
  dplyr::distinct(isotype, env_trait, .keep_all = T)

test <- ggplot(temp_sweep_clim) +
  aes(x = pop_assignment, y = value, fill = fill = factor(pop_assignment, levels = c("A", "B", "D", "E", "C", "F", "G"))) +
  #scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 1.5, shape = ifelse(temp_sweep_clim$Hawaiian == T, 21, 25)) +
  facet_wrap(~env_trait, scales = "free") +
  scale_fill_manual(values=c(ancestry.colours)) +
  theme_bw() +
  theme(legend.position="none") +
  labs(y = "temperature", x = "admixture group")
test


# checking correlation of fraction admixed with average sequencing depth and number of strains in isotype
check_index <- admix %>%
  dplyr::mutate(frac_mixed = 1-max_pop_frac)

# load depth summary for isotypes in 330 set
depth_summary <- data.table::fread("~/Hawaii_Manuscript/data/depth_summary.tsv") %>%
  dplyr::rename(isotype = INDV)

depth_summary_cendr <- data.table::fread("~/Hawaii_Manuscript/data/depth_summary_CeNDR_hard-filtered.tsv") %>%
  dplyr::rename(isotype = INDV, MEAN_DEPTH_HF = MEAN_DEPTH, N_SITES_HF = N_SITES)
  

depth_summary <- left_join(depth_summary, depth_summary_cendr)

wi_data <- data.table::fread("~/Hawaii_Manuscript/data/WI_strain_list.csv") %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(strain_num = n()) %>%
  dplyr::distinct(isotype, .keep_all = T) %>%
  dplyr::select(isotype, strain_num)

join1 <- full_join(depth_summary, wi_data)
join2 <- full_join(check_index, join1) %>%
  dplyr::mutate(strain_num_group = case_when(
                                            strain_num == 1 ~ "1",
                                            strain_num <= 3 & strain_num > 1 ~ "2-3",
                                            strain_num > 3 ~ "3+")) %>%
  dplyr::mutate(admixed = ifelse(max_pop_frac < 0.95, "admix", "no_admix"))

# Plot correlations
strain_num <- ggplot(join2) +
  aes(x = strain_num_group, y = frac_mixed) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 1.5, shape = 21) +
  
  #labs(y = "Fraction most common global haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  labs(y="Fraction admixed", x="number of strains in isotype") 
strain_num

# Run non-parametric kruskal-wallis test with Dunn's post-hoc test. Kruskal–Wallis test does not assume a normal distribution of the residuals
dtest <- dunnTest(join2$frac_mixed ~ join2$strain_num_group, method = "bonferroni")

# Run enrichment test on admixed strains
# shape data for analysis at Caenorhabditis level
enrich_test <- join2 %>%
  dplyr::filter(!is.na(admixed)) %>%
  dplyr::group_by(strain_num_group, admixed) %>%
  dplyr::mutate(admix_count = n()) %>%
  dplyr::distinct(strain_num_group, admix_count, admixed) %>%
  dplyr::ungroup() %>%
  tidyr::spread(admixed, admix_count) %>%
  column_to_rownames(var = "strain_num_group") %>%
  as.matrix(.)

# show p-values without scientific notation
options(scipen=999)

# full chi squared test to see if there are differernces 
chisq.test(enrich_test, simulate.p.value = TRUE)

# Post-hoc pairwise Fisher tests with pairwise.table
pairwiseNominalIndependence(enrich_test,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)

# plot with coverage depth
cov_depth <- ggplot(join2) +
  aes(x = MEAN_DEPTH, y = frac_mixed) +
  #geom_boxplot(outlier.shape = NA) +
  geom_point(width = 0.25, size = 1.5, shape = 21) +
  
  #labs(y = "Fraction most common global haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  labs(y="Fraction admixed", x="mean depth")  +
  geom_smooth()
cov_depth

cov_strain <- ggplot(join2) +
  aes(x = strain_num, y = MEAN_DEPTH) +
  #geom_boxplot(outlier.shape = NA) +
  geom_point(width = 0.25, size = 1.5, shape = 21) +
  
  #labs(y = "Fraction most common global haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  labs(y="Mean depth", x="number of strains in isotype") +
  geom_smooth()
cov_strain



