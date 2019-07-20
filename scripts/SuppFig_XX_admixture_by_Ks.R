#!/usr/bin/env Rscript
library(pophelper)
library(tidyverse)
library(ggthemes)
library(stats)
library(ggrepel)
library(scatterpie)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#set VCF path
VCF_PATH="data/ANNOTATE_VCF/LD_0.8.vcf.gz"

# hawaii strains
strain_islands <- c("XZ1514" = "#E69F00", "XZ1516" = "#E69F00","XZ1513" = "#E69F00","ECA372" = "#E69F00","ECA701" = "#E69F00","XZ1515" = "#E69F00",
                    "CB4856" = "#56B4E9",
                    "ECA369" = "#009E73","ECA738" = "#009E73",
                    "QX1792" = "#0072B2", "QX1794" = "#0072B2", "QX1793" = "#0072B2", "QX1791" = "#0072B2", "ECA740" = "#0072B2", "ECA741" = "#0072B2", "ECA363" = "#0072B2", "ECA743" = "#0072B2", "ECA742" = "#0072B2",
                    "ECA760" = "#CC79A7","ECA768" = "#CC79A7","ECA777" = "#CC79A7","ECA706" = "#CC79A7","ECA705" = "#CC79A7","ECA703" = "#CC79A7","ECA807" = "#CC79A7","ECA778" = "#CC79A7",
                    "ECA812" = "#CC79A7","ECA710" = "#CC79A7","ECA744" = "#CC79A7","ECA745" = "#CC79A7","ECA732" = "#CC79A7","ECA733" = "#CC79A7","ECA746" = "#CC79A7","DL238" = "#CC79A7",
                    "ECA347" = "#CC79A7","ECA730" = "#CC79A7","ECA724" = "#CC79A7","ECA722" = "#CC79A7","ECA189" = "#CC79A7","ECA191" = "#CC79A7","ECA723" = "#CC79A7","ECA712" = "#CC79A7",
                    "ECA396" = "#CC79A7")

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# pie_chart define a color pallette with maximal contrast
ancestry.colours.pie <- c("perc_a"="gold2", "perc_b"="plum4","perc_c"= "darkorange1", 
                          "perc_d"="lightskyblue2", "perc_e"="firebrick","perc_f"= "burlywood3",
                          "perc_g"="gray51",
                          "perc_h" = "springgreen4",
                          "perc_i" = "lightpink2",
                          "perc_j" = "deepskyblue4",
                          "perc_k" = "mediumpurple4",
                          "perc_l" = "orange",
                          "perc_m"= "maroon","perc_n"= "yellow3","perc_o"= "brown4", 
                          "perc_p"="yellow4", "perc_q"="sienna4", "perc_r"="chocolate", "perc_s"="gray19")

# generate sample list
# vcf name says 330, but it is a lie
#system("bcftools query -l data/ANNOTATE_VCF/Ce330_annotated.vcf.gz > data/ANNOTATE_VCF/samples.txt")

# get Hawaii strains
hi_only_samples <- read.csv(file = "data/fulcrum/hawaii_isotypes.csv") 

# get processed K summarys
k7 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K7_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 7)

k8 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K8_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 8)

k9 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K9_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 9)

k10 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K10_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 10)

k11 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K11_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 11)

k12 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K12_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 12)

k13 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K13_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 13)

k14 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K14_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 14)

k15 <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K15_Processed_Ancestry.tsv", header = T) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE")) %>%
  dplyr::mutate(K_num = 15)

# define hawaiian clusters for diferent Ks
k7_hipops <- c("C", "F", "G")
k8_hipops <- c("C", "F", "D")
k9_hipops <- c("C", "F", "H", "I") 
k10_hipops<- c("I", "E", "F", "G", "B")
k11_hipops<- c("A", "E", "H", "F")
k12_hipops<- c("C", "K", "E", "B")
k13_hipops<- c("J", "I", "H", "B", "E")
k14_hipops<- c("A", "E", "D", "L", "M")
k15_hipops<- c("B", "O", "M", "F")

# Join Ks into one dataframe
k_full <- list(k7, k8, k9, k10, k11, k12, k13, k14, k15) %>%
  reduce(full_join) %>%
  dplyr::mutate(admixed_10 = ifelse(max_pop_frac > 0.9, "TRUE", "FALSE"),
                admixed_5 = ifelse(max_pop_frac > 0.95, "TRUE", "FALSE")) %>%
  dplyr::filter(Hawaiian == "TRUE") %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(max_pop_sum = sum(max_pop_frac)) %>% # measures degree of admixture across all Ks; lower value is more admixed 
  dplyr::group_by(isotype, K_num) %>%
  dplyr::mutate(non_hi_admix_frac = ifelse(K_num == 7, 1-sum(C, F, G),
                                           ifelse(K_num == 8, 1-sum(C, F, D),
                                                  ifelse(K_num == 9, 1-sum(C, F, H, I),
                                                         ifelse(K_num == 10, 1-sum(I, E, F, G, B),
                                                                ifelse(K_num == 11, 1-sum(A, E, H, F),
                                                                       ifelse(K_num == 12, 1-sum(C, K, E, B),
                                                                              ifelse(K_num == 13, 1-sum(J, I, H, B, E),
                                                                                     ifelse(K_num == 14, 1-sum(A, E, D, L, M), 1-sum(B, O, M,F)))))))))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(max_pop_sum) %>%
  tidyr::gather(cluster, frac_cluster, -isotype, -max_pop_frac, -pop_assignment, -Hawaiian, -K_num, -admixed_10, -admixed_5, -max_pop_sum, -non_hi_admix_frac) %>%
  dplyr::mutate(frac_cluster = ifelse(is.na(frac_cluster), 0, frac_cluster),
                perc_cluster = frac_cluster*100)

# get admix order from above
isotype_admix_order <- k_full %>%
  dplyr::distinct(isotype)

# plot pie charts in grid
pop_assignment_pie_chart <- ggplot(k_full, aes(x="", y=perc_cluster, fill=cluster))+
  geom_bar(width = 1, stat = "identity") +
  facet_grid(factor(isotype, levels = isotype_admix_order$isotype)~K_num) +
  coord_polar("y", start=0) +
  labs(y = "", x = "", fill = "") +
  scale_fill_manual(values = ancestry.colours) +
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank())
#geom_text(aes(x=1, y = value, label=ifelse(group == "fraction_hi", num_hi_isotypes_in_pop, num_isotypes_in_pop-num_hi_isotypes_in_pop)))
pop_assignment_pie_chart

ggsave(paste("plots/ADMIXTURE_pie_charts_LD8_K7-15.png"), width = 10, height = 30)
ggsave(paste("plots/ADMIXTURE_pie_charts_LD8_K7-15.pdf"), width = 10, height = 30)

# get non-hi admix order from above
isotype_non_hi_admix_order <- k_full %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(strain_sum_non_hi_admix_frac = sum(non_hi_admix_frac)) %>%
  dplyr::mutate(average_nonHW_admix_frac = mean(non_hi_admix_frac)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(strain_sum_non_hi_admix_frac) %>%
  dplyr::distinct(isotype, .keep_all=T) %>%
  dplyr::select(isotype, average_nonHW_admix_frac) 

# plot heat map for degree admixed
p <- ggplot(k_full %>% dplyr::mutate(isotype = factor(isotype, levels = isotype_non_hi_admix_order$isotype),
                                     non_hi_admix_frac),
            aes(K_num, isotype)) + 
  geom_tile(aes(fill = non_hi_admix_frac), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue")
p

ggsave(paste("plots/ADMIXTURE_heatmap_LD8_K7-15.png"), width = 10, height = 10)
ggsave(paste("plots/ADMIXTURE_heatmap_LD8_K7-15.pdf"), width = 10, height = 10)

####################
# plot heatmaps for non-HW admixture for K=11 and K=13
k_full_hz <- k_full %>%
  dplyr::filter(K_num %in% c(11, 13) & non_hi_admix_frac>0.05) %>%
  dplyr::mutate(global_pop = ifelse(K_num == 11 & cluster %in% c("A", "E", "H", "F"), "FALSE",
                                ifelse(K_num == 11 & !(cluster %in% c("A", "E", "H", "F")), "TRUE",
                                       ifelse(K_num == 13 & cluster %in% c("J", "I", "H", "B", "E"), "FALSE", "TRUE"))))

# plot pie charts in grid for K=11 and K=13 for just the non-hawaiian admixture values
p2 <- ggplot(k_full_hz, aes(x="", y=perc_cluster, fill=cluster, alpha=global_pop))+
  geom_bar(width = 1, stat = "identity") +
  facet_grid(K_num~factor(isotype, levels = rev(isotype_non_hi_admix_order$isotype))) +
  coord_polar("y", start=0) +
  labs(y = "", x = "", fill = "") +
  scale_fill_manual(values = ancestry.colours) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank())
#geom_text(aes(x=1, y = value, label=ifelse(group == "fraction_hi", num_hi_isotypes_in_pop, num_isotypes_in_pop-num_hi_isotypes_in_pop)))
p2

ggsave(paste("plots/ADMIXTURE_nonHW_ADMIXED_LD8_K11_13.png"), width = 10, height = 4)
ggsave(paste("plots/ADMIXTURE_nonHW_ADMIXED_LD8_K11_13.pdf"), width = 10, height = 4)

##################
# plot hybrid zone maps for K=11 and K=13

# load in wi strain info for lat long coordinates for map
df <- data.table::fread("data/WI_strain_list.csv")

# load strain names, generated by Figure7.R
sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples.txt", header = F) %>% dplyr::pull(V1))

#Set K=11 as admixture qfile
k = 11
qfile <- as.data.frame(pophelper::readQ("data/ADMIXTURE_LD8_K12/BEST_K/LD_0.8_MAF_0.004.11.Q"))
colnames(qfile) <- LETTERS[1:k]
qfile <- qfile %>%
  dplyr::mutate(samples = as.vector(sample_names))

# make long and determin order of plotting
long_admix_pops <- qfile %>%
  dplyr::mutate(samples = sample_names) %>%
  dplyr::filter(samples %in% names(strain_islands)) %>%
  tidyr::gather(cluster, frac_cluster, -samples) %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(max_frac = max(frac_cluster)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))

# load admix info for full 276_set for K = 12
admix <- data.table::fread("data/ADMIXTURE_LD8_K12/BEST_K/K11_Processed_Ancestry.tsv") %>%
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

#5 plot HI map file supplement  
strain_pops <- dplyr::filter(long_admix_pops, frac_cluster == max_frac) %>%
  dplyr::rename(strain=samples)

isolation_info <- df %>%
  #dplyr::filter(reference_strain == 1)%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude, state, country)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(strain_pops,.,by="strain")%>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::distinct(strain, long, lat, .keep_all = TRUE)

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)

dev.off()
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

HW_plot <- ggplot()+ geom_map(data=world, map=world,
                                 aes(x=long, y=lat, map_id=region),
                                 color="black", fill="#E6E6E6", size=0.25)+
  scale_fill_manual(values = ancestry.colours,name = "Population")+
  theme_map()+
  # geom_label_repel(aes(long, lat, label = strain, fill = cluster),
  #                  data=dplyr::filter(isolation_info, state == "Hawaii"),
  #                  fontface = 'bold', color = 'white',
  #                  segment.color = 'cyan')+
  geom_point(aes(long, lat, fill = pop_assignment),
                        data=full_isolation_info %>% dplyr::filter(!(isotype %in% c("XZ1515", "ECA372", "ECA738", "ECA189", "XZ1513", "QX1794", "ECA710"))) %>%
                          dplyr::mutate(pop_assignment = factor(pop_assignment, level = c("A", "D", "F", "G", "H", "I",
                                                                                          "J", "B", "C", "E", "K"))) %>%
                          dplyr::arrange(pop_assignment),
                        shape = 21, size = 3) +
  theme(legend.position = "none") +
  coord_quickmap(xlim = c(-160,-155), ylim = c(19,22.3))  # coord_quickmap keeps correct aspect ratio
  HW_plot + geom_scatterpie(aes(x=long, y=lat),
                                    data=dplyr::filter(full_isolation_info, state == "Hawaii" & isotype %in% c("XZ1515", "ECA372", "ECA738", "ECA189", "XZ1513", "QX1794", "ECA710")), cols=LETTERS[1:11],  alpha=.8)

#pull global C poplist
k11_popC <- full_isolation_info %>%
  dplyr::filter(state %in% c("california", "oregon", "washington") & C> 0.5)

  
dev.off()
states <- map_data("state")
states <- states %>%
  filter(region %in% c("california", "oregon", "washington"))   # intercourse other states

 CAL_plot <- ggplot()+ geom_map(data=states, map=states,
                                aes(x=long, y=lat, map_id=region),
                                color="black", fill="#E6E6E6", size=0.25)+
    scale_fill_manual(values = ancestry.colours,name = "Population")+
    theme_map()+
    # geom_label_repel(aes(long, lat, label = strain, fill = cluster),
    #                  data=dplyr::filter(isolation_info, state == "Hawaii"),
    #                  fontface = 'bold', color = 'white',
    #                  segment.color = 'cyan')+
    geom_point(aes(long, lat, fill = pop_assignment),
               data=full_isolation_info %>% dplyr::filter(state %in% c("California", "Oregon", "Washington")) %>%
                 dplyr::mutate(pop_assignment = factor(pop_assignment, level = c("A", "D", "F", "G", "H", "I",
                                                                                 "J", "B", "C", "E", "K"))) %>%
                 dplyr::arrange(pop_assignment),
               shape = 21, size = 3) +
    theme(legend.position = "none") +
    coord_quickmap()  # coord_quickmap keeps correct aspect ratio
 CAL_plot + geom_scatterpie(aes(x=long, y=lat),
                            data=dplyr::filter(full_isolation_info, state == "California" &isotype %in% c("XZ1515", "ECA372", "ECA738", "ECA189", "XZ1513", "QX1794", "ECA710")), cols=LETTERS[1:11],  alpha=.8)  
  
# ggsave(HW_plot, filename = "plots/K11_HW_MAP_Admix_pie_admixed.pdf",
#        height = 10,
#        width = 10)

# get admix
full_isolation_info <- df %>%
  #dplyr::filter(reference_strain == 1) %>%
  dplyr::select(isotype, long = longitude, lat = latitude, state, country) %>%
  dplyr::filter(lat != "None") %>%
  dplyr::left_join(.,admix,by="isotype") %>%
  dplyr::distinct(isotype, long, lat, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(pop_assignment))

full_world_plot <- ggplot()+ geom_map(data=world, map=world,
                                      aes(x=long, y=lat, map_id=region),
                                      color="black", fill="#E6E6E6", size=0.25)+
  scale_fill_manual(values = ancestry.colours,name = "Population")+
  theme_map()+
  # geom_point(aes(long, lat, fill = pop_assignment),
  #            data=full_isolation_info %>%
  #              dplyr::mutate(pop_assignment = factor(pop_assignment, level = c("A", "D", "F", "G", "H", "I",
  #                                                                              "J", "L", "B", "C", "E", "K"))) %>%
  #              dplyr::arrange(pop_assignment),
  #            shape = 21, size = 3) +
  theme(legend.position = "none") +
  coord_quickmap() 
full_world_plot + geom_scatterpie(aes(x=long, y=lat),
                    data=full_isolation_info, cols=LETTERS[1:11], color=NA, alpha=.8) 
# ggsave(full_world_plot, filename = "plots/K11_world_MAP_Admix.pdf", useDingbats = F,
#        height = 7.5,
#        width = 7.5)

#####


