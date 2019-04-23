#Generate plots for Hawaii 2017 collection
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

####################################################
#  A: Define Functions                             #
####################################################
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

####################################################
#  B: Setup data frames and add color palettes     #
####################################################
# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# color palettes from pals package: https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html
species_palette <- c("C. elegans" = "#BE0032", #7
                     "C. sp. 53" =  "#875692", #4 
                     "C. tropicalis" = "#F38400", #5 
                     "Panagrolaimus sp." = "#C2B280", #8
                     "Oscheius sp." = "#F3C300", #3
                     "C. briggsae" = "#A1CAF1", #6 
                     "Other PCR +" = "#008856", #10
                     "PCR -" = "#848482", #9
                     "Not genotyped" = "#F2F3F4", #1
                     "No Worm" = "#222222") #2
                     
substrate_palette <- c("Leaf litter" = "#E68FAC",
                       "Fruit/nut/veg" = "#0067A5",
                       "Flower" = "#DCD300",
                       "Fungus" = "#604E97",
                       "Compost" = "#F6A600",
                       "Other" = "#B3446C")

island_palette <- c("Kauai" = "#E69F00",
                    "Oahu" = "#56B4E9",
                    "Molokai" = "#009E73",
                    "Maui" = "#F0E442",
                    "Big Island" = "#D55E00")

ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")


####################################################
#  Figure 3: Distribution of envi. parameters      #
####################################################
# Kruskal-Wallis non-parametric tests on distinct collections, i.e. if a c-label had three c. elegans and 2 C. briggsae isolated from it we are counting it once for C. elegaans and once for C. briggsae.
stat_df <- cso %>%
  dplyr::filter(pcr_rhpositive == 1) %>%
  dplyr::mutate(species_family = ifelse(spp_id %in% c("Chabertia ovina", 
                                                      "Choriorhabditis cristata",
                                                      "Choriorhabditis sp.",
                                                      "Heterhabditis zealandica",
                                                      "Mesorhabditis sp.",
                                                      "no match",
                                                      "C. kamaaina",
                                                      "Rhabditis terricola", 
                                                      "Rhanditis tericola",
                                                      "Teratorhabditis sp.",
                                                      "Unknown",
                                                      "unknown",
                                                      NA), 
                                        "PCR +", spp_id)) %>%
  dplyr::distinct(c_label, spp_id, .keep_all = T) %>%
  dplyr::select(c_label, s_label, species_family, substrate_temperature, substrate_moisture, ambient_humidity, 
                ambient_temperature, altitude) %>%
  dplyr::mutate(species_family = forcats::as_factor(species_family),
                species_family = forcats::fct_relevel(species_family,
                                                      "C. elegans",
                                                      "C. briggsae",
                                                      "Oscheius sp.",
                                                      "C. tropicalis",
                                                      "C. sp. 53",
                                                      "Panagrolaimus sp.",
                                                      "PCR +")) %>%
  tidyr::gather(env_par, value, substrate_temperature, substrate_moisture, ambient_humidity, ambient_temperature, altitude)%>%
  dplyr::group_by(env_par) %>%
  dplyr::mutate(KM_pvalue = kruskal.test(value ~ species_family)[[3]])

# perform multiple comparisions test using Dunn's test with pvalues adjusted with Bonferroni method
Dunn_list <- list()

for (e in 1:length(unique(stat_df$env_par))){
  KM_df <- stat_df %>%
    dplyr::filter(env_par == (unique(stat_df$env_par)[e]))
  
  D_test <- dunnTest(KM_df$value ~ KM_df$species_family, method = "bonferroni") 
  Dunn_list[[unique(stat_df$env_par)[e]]] <- D_test
}

# Plot data as box plots for each environmental parameter
env_par_palette <- c("#ff0000", "#ff8000", "#a6aeff", "#ffff00", "#0080ff", "#ff69b4", "#33CC00")
names(env_par_palette) = c("C. elegans", "C. tropicalis", "Panagrolaimus sp.", "Oscheius sp.", "C. briggsae", "C. sp. 53", "PCR +")

plot_atemp <- ggplot(data = stat_df %>% dplyr::filter(env_par == "ambient_temperature")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(env_par_palette)) +
  geom_jitter(aes(x = species_family, y = value), size = .5, width = .25) +  
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black")) + 
  labs(fill = "Species", x = "", y = "Ambient temperature (°C)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_stemp <- ggplot(data = stat_df %>% dplyr::filter(env_par == "substrate_temperature")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(env_par_palette)) +
  geom_jitter(aes(x = species_family, y = value), size = .5, width = .25) +  
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black")) + 
  labs(fill = "Species", x = "", y = "Substrate temperature (°C)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_ahum <- ggplot(data = stat_df %>% dplyr::filter(env_par == "ambient_humidity")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(env_par_palette)) +
  geom_jitter(aes(x = species_family, y = value), size = .5, width = .25) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black")) + 
  labs(fill = "Species", x = "", y = "Ambient humidity (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_smoist <- ggplot(data = stat_df %>% dplyr::filter(env_par == "substrate_moisture")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(env_par_palette)) +
  geom_jitter(aes(x = species_family, y = value), size = .5, width = .25) +  
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black")) + 
  labs(fill = "Species", x = "", y = "Substrate moisture (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_elev <- ggplot(data = stat_df %>% dplyr::filter(env_par == "altitude")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(env_par_palette)) +
  geom_jitter(aes(x = species_family, y = value), size = .5, width = .25) +  
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black")) + 
  labs(fill = "Species", x = "", y = "Elevation (m)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

legend <- get_legend(plot_elev + theme(legend.position="right"))

plots <- cowplot::plot_grid(plot_atemp, plot_stemp, plot_ahum, plot_smoist, plot_elev, legend, labels = c("A", "B", "C", "D", "E", ""), ncol = 2, nrow = 3, align = "hv", axis = "l")  

plots

#ggsave('figure/20180807_env_pars.pdf', width = 12, height = 8.6)
#ggsave('figure/20180807_env_pars.png', width = 12, height = 8.6)
####################################################
#  Figure 9: ADMIX matrix plot                     #
####################################################
#  ~ ~ ~ ~ command to extract CV value from the ADMIXTURE log files ~ ~ ~ ~ #

# grep -h CV log*.out > Ksummary.txt

base.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/Manuscripts/Hawaii_Manuscript/"
setwd(base.dir)

# generate sample list
system("bcftools query -l Data/VCF/hard-filtered-used-for-trees/249-hawaii.genome.vcf.gz > Data/ADMIXTURE/249_Hawaii/samples.txt")

# strain names
# get sample information
samples <- sort(data.table::fread(paste0(base.dir, "Data/ADMIXTURE/249_Hawaii/samples.txt"),header = F) %>% dplyr::pull(V1))
k_summary <- data.table::fread(paste0(base.dir, "Data/ADMIXTURE/249_Hawaii/CV_Summary/admix_replicates_CV.tsv"),header = T) 

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ PLOT ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

setwd(base.dir)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Extract minimal information for minimal K ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

plot.dir <- glue::glue("{base.dir}Plots/ADMIXTURE/249_Hawaii")
dir.create(plot.dir)

ggplot(k_summary)+
  aes(x = factor(K), y = CV)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = .1)+
  theme_bw()+
  labs(x = "K")

# ggsave(glue::glue("{plot.dir}/K_v_CV.pdf", height = 6, width = 12))

# plot admixture
base.admixture.plot.dir <- glue::glue("{base.dir}Plots/ADMIXTURE/249_Hawaii/")
base.admixture.data.dir <- glue::glue("{base.dir}Processed_Data/ADMIXTURE/249_Hawaii/")
dir.create(base.admixture.plot.dir)
dir.create(base.admixture.data.dir)

K <- 6

setwd(glue::glue("{base.dir}Data/ADMIXTURE/249_Hawaii/BEST_K/"))

# load Q files
qfile_name <- grep(pattern = paste0("Q$"), value = T, x = list.files())
qfile <- pophelper::readQ(files = qfile_name)

# look at qfil
# label Q file rownames and colnames
rownames(qfile[[1]]) <- samples
colnames(qfile[[1]]) <- LETTERS[1:K]

# make Q file df
hi_only_samples <- read.csv(file = "~/Dropbox/AndersenLab/Hawaii_manuscript/data/fulcrum/hawaii_isotypes.csv") %>%
  dplyr::arrange(strain) %>%
  dplyr::filter(strain %in% samples) %>%
  dplyr::select(samples = strain)


hi_only_qfile_df <- qfile[[1]] %>%
  dplyr::filter(row.names(.) %in% hi_only_samples$samples) %>%
  dplyr::mutate(samples = hi_only_samples$samples) %>%
  tidyr::gather(cluster, frac_cluster, - samples) %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(max_frac = max(frac_cluster)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, max_frac) %>%
  dplyr::mutate(samples = factor(samples))

plot_order <- hi_only_qfile_df %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = samples))


admix_plot <- ggplot(hi_only_qfile_df) +
  geom_bar(stat = "identity", aes(x = factor(samples, levels =plot_order$samples), y = frac_cluster, fill = cluster)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(fill = "", x = "", y = "ADMIXTURE") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.x=element_blank(),    
        axis.text.y=element_blank(),
        axis.title.y = element_text(angle = 0, vjust = .5),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

admix_plot_no_legend <- admix_plot +
  theme(legend.position="none",
  axis.line=element_blank(),
  plot.margin = unit(c(0,0,0,0), units = "cm"))

grid.newpage()
admix_legend <- cowplot::get_legend(admix_plot)
grid.draw(admix_legend) +
  theme(plot.margin = unit(c(0,0,0,0), units = "cm"))

# Filter cso to isolates assigned to an isotype
hm_hi1 <- cso %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::select(isotype,
                strain,
                c_label,
                s_label,
                substrate,
                island,
                altitude,
                latitude,
                longitude,
                ambient_temperature,
                substrate_temperature,
                ambient_humidity,
                substrate_moisture)

# load hawaii isolates found before hawaii 2017 trip
old_hi_isolates <- googlesheets::gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  googlesheets::gs_read("WI C. elegans") %>%
  dplyr::filter(isotype %in% hi_only_samples$samples & !(isotype %in% hm_hi1$isotype)) %>%
  dplyr::mutate(altitude = as.numeric(NA)) %>%
  dplyr::mutate(longitude = as.numeric(longitude),
                substrate_temperature = substrate_temp,
                ambient_temperature = ambient_temp)

# join old and new hawaiian isolates
hm_hi2 <- full_join(hm_hi1, old_hi_isolates)

#set substrate to manuscript categories
substrate_merge <- c("Fruit",
                     "Rotting fruit",
                     "Nut",
                     "Rotting nut",
                     "Rotting vegetable")

hm_hi3 <- hm_hi2 %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(substrate = ifelse(substrate %in% substrate_merge, "Fruit/nut/vegetable", substrate)) %>%
  dplyr::mutate(substrate = ifelse(substrate %in% c("Rotting fungus"), "Fungus", substrate)) %>%
  dplyr::mutate(substrate = ifelse(is.na(substrate), "Other", substrate)) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Rotting bark",
                                                                               "Soil",
                                                                               "Grass"), 
                                                              "Other",substrate))))) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(main_substrate = names(sort(summary(as.factor(fixed_substrate)), decreasing=T))[1]) %>% # find most common factor in group
  dplyr::ungroup()

#set altitude from gps except for unknown collection locations
options(geonamesUsername="katiesevans")
options(geonamesHost="api.geonames.org")
hm_hi4 <- hm_hi3 %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(altitude = ifelse(is.na(altitude),
                                  geonames::GNsrtm3(latitude, longitude)$srtm3,
                                  altitude)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(altitude = ifelse(isotype %in% c("CB4856", "DL238"), NA, altitude))

#set islands from gps
hm_hi4$island <- "?"
hm_hi4[filter_box(hm_hi4$longitude, hm_hi4$latitude, c(-158.3617,21.1968,-157.5117,21.7931)), "island"] <- "Oahu"
hm_hi4[filter_box(hm_hi4$longitude, hm_hi4$latitude, c(-159.9362, 21.6523, -159.1782, 22.472)), "island"] <- "Kauai"
hm_hi4[filter_box(hm_hi4$longitude, hm_hi4$latitude, c(-157.327, 21.0328, -156.685, 21.2574)), "island"] <- "Molokai"
hm_hi4[filter_box(hm_hi4$longitude, hm_hi4$latitude, c(-156.7061, 20.4712, -155.9289, 21.0743)), "island"] <- "Maui"
hm_hi4[filter_box(hm_hi4$longitude, hm_hi4$latitude, c(-156.1346, 18.6619, -154.6985, 20.4492)), "island"] <- "Big Island"

# Identify most common island for isotypes
hm_hi5 <- hm_hi4 %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(main_island = names(sort(summary(as.factor(island)), decreasing=T))[1]) %>%
  dplyr::ungroup()

# Identify mean values for continous variables for isotypes
hm_hi6 <- hm_hi5 %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise_if(is.numeric, funs(mean), na.rm = TRUE) %>%
  dplyr::select(isotype,
                altitude,
                ambient_temperature,
                substrate_temperature,
                ambient_humidity,
                substrate_moisture) %>%
  dplyr::ungroup()

# Combine isotype data for plotting continous data
hm_hi_proc <- full_join(hm_hi5 %>% dplyr::select(isotype, strain, main_island, main_substrate), hm_hi6) %>%
  dplyr::distinct(isotype, .keep_all = TRUE) %>%
  dplyr::select(-strain) %>%
  dplyr::mutate(altitude_s = scale_this(altitude),
                ambient_temperature_s = scale_this(ambient_temperature),
                substrate_temperature_s = scale_this(substrate_temperature),
                ambient_humidity_s = scale_this(ambient_humidity),
                substrate_moisture_s = scale_this(substrate_moisture)) %>%
  tidyr::gather(trait, value, -isotype)

# plot continous variables on heat map
hm_plot <- ggplot(hm_hi_proc %>% 
                    dplyr::filter(trait %in% c("altitude_s",
                                               "ambient_temperature_s",
                                               "substrate_temperature_s",
                                               "ambient_humidity_s",
                                               "substrate_moisture_s"))) +
  aes(x = factor(isotype, levels =plot_order$samples), y = trait, fill = as.numeric(value)) +
  geom_tile(aes(x = factor(isotype, levels =plot_order$samples), y = trait, fill = as.numeric(value))) +
  theme_bw() + 
  scale_fill_viridis(name="") +
  coord_equal() +
  theme_bw() +
  theme(axis.line=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

hm_plot_no_legend <- hm_plot +
  theme(legend.position="none",
        plot.margin = unit(c(0,0,0,0), units = "cm"))

hm_plot_legend <- cowplot::get_legend(hm_plot)

# 2 
# plot discrete variable, island, substrate
sub_plot <- ggplot(hm_hi_proc %>% 
                     dplyr::filter(trait == "main_substrate") %>%
                     dplyr::mutate(value = factor(value, levels=names(substrate_palette)))) +
  geom_tile(aes(x = factor(isotype, levels =plot_order$samples), y = trait, fill = value)) + #factor(value, levels=names(substrate_palette))
  theme_bw() + 
  coord_equal() +
  labs(fill="") +
  scale_fill_manual(values = substrate_palette) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

sub_plot_no_legend <- sub_plot +
  theme(legend.position="none",
  plot.margin = unit(c(0,0,0,0), units = "cm"))

sub_plot_legend <- cowplot::get_legend(sub_plot)

#3
# plot islands
isl_plot <- ggplot(hm_hi_proc %>% 
                     dplyr::filter(trait == "main_island") %>%
                     dplyr::mutate(value = factor(value, levels=names(island_palette)))) +
  geom_tile(aes(x = factor(isotype, levels =plot_order$samples), y = trait, fill = value)) +
  theme_bw() + 
  coord_equal() +
  labs(fill="") +
  scale_fill_manual(values=island_palette) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

isl_plot_no_legend <- isl_plot +
  theme(legend.position="none",
  plot.margin = unit(c(0,0,0,0), units = "cm"))

isl_plot_legend <- cowplot::get_legend(isl_plot)

#4
# add all together
admixture_panel <- cowplot::plot_grid(admix_plot_no_legend,
                                              admix_legend,
                                              isl_plot_no_legend,
                                              isl_plot_legend,
                                              sub_plot_no_legend,
                                              sub_plot_legend,
                                              hm_plot_no_legend,
                                              hm_plot_legend,
                                              ncol = 2, nrow = 4,
                                              rel_heights = c(.8, .5, .5, 1),
                                              rel_widths = c(1, .2),
                                              align = "vh",
                                              axis = 'l')
# reassign working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

ggsave('plots/admixture_panel.pdf', width = 20, height = 10)

###############################
# multiple isotypes histogram #
###############################

# number of isotypes on a c_plate (2 collections are excluded because no isotypes were recovered from them do to extinction before sequencing)
num_isotypes_per_sample <- cso %>%
  dplyr::filter(!is.na(s_label)) %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(c_label, isotype, .keep_all = TRUE) %>%
  dplyr::select(isotype, spp_id, c_label, s_label, substrate, latitude, longitude) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(distinct_isotypes = n()) %>%
  dplyr::distinct(c_label, .keep_all=TRUE)


# histogram of isotypes per C plate (All but one c-plate with mutliple isotypes come from Fig tree)
isotype_hist <- ggplot(num_isotypes_per_sample) +
  aes(x = distinct_isotypes) +
  geom_histogram(binwidth = 1) +
  labs(y = "Count", x = "Distinct isotypes per C-plate")
isotype_hist

ggsave('plots/multiple_isotypes_per_cplate.pdf', width = 85, height = 85, units = "mm")

################################################################
# find admixture assignment population  for all strains in set #
################################################################

admix_pops <- qfile[[1]] %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(isotype = rowname) %>%
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
  tidyr::spread(pop, frac_pop)


