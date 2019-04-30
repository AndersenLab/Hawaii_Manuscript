#library(imager)
library(cowplot)
#library(jpeg)
#library(grid)
#library(gridExtra)
library(scales)
library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# Color palettes
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

# All Caenorhabditis collections ---------------------------------------------------------
pop_size_caen_by_sub <- cso %>%
  dplyr::filter(spp_id %in% c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae")) %>%
  dplyr::mutate(species_family = spp_id) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass",
                                                                               "Isopod",
                                                                               "Slug",
                                                                               "Millipede"), 
                                                              "Other",substrate))))) %>%
  dplyr::mutate(fixed_approximate_number_of_worms = ifelse(approximate_number_of_worms == "Very Few (1-3)", "1-3",
                                                           ifelse(approximate_number_of_worms == "Few (4-10)", "4-10",
                                                                  ifelse(approximate_number_of_worms == "Some (11-25)", "11-25",
                                                                         ifelse(approximate_number_of_worms == "Proliferating (25+)", "25+",
                                                                                approximate_number_of_worms))))) %>%
  dplyr::distinct(c_label, .keep_all=T) %>%
  dplyr::group_by(fixed_substrate, fixed_approximate_number_of_worms) %>%
  dplyr::mutate(worm_per_popsize = n()) %>%
  dplyr::distinct(fixed_substrate, fixed_approximate_number_of_worms, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_approximate_number_of_worms) %>%
  dplyr::mutate(total_popsize = sum(worm_per_popsize), perc_worm_popsize = worm_per_popsize / total_popsize * 100) %>%
  dplyr::arrange(total_popsize) %>%
  dplyr::select(fixed_substrate, fixed_approximate_number_of_worms, worm_per_popsize, total_popsize, perc_worm_popsize, species_family) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fixed_approximate_number_of_worms = factor(fixed_approximate_number_of_worms, levels = c("1-3", "4-10", "11-25", "25+")),
                fixed_substrate = factor(fixed_substrate, levels = c("Leaf litter", "Fruit/nut/veg", "Flower","Fungus", "Other","Isopod", "Slug", "Millipede")))

pop_size_caen_by_sp  <- cso %>%
  dplyr::filter(spp_id %in% c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae")) %>%
  dplyr::mutate(species_family = spp_id) %>%
  dplyr::mutate(fixed_approximate_number_of_worms = ifelse(approximate_number_of_worms == "Very Few (1-3)", "1-3",
                                                           ifelse(approximate_number_of_worms == "Few (4-10)", "4-10",
                                                                  ifelse(approximate_number_of_worms == "Some (11-25)", "11-25",
                                                                         ifelse(approximate_number_of_worms == "Proliferating (25+)", "25+",
                                                                                approximate_number_of_worms))))) %>%
  dplyr::mutate(fixed_approximate_number_of_worms = factor(fixed_approximate_number_of_worms, levels = c("1-3", "4-10", "11-25", "25+"))) %>%
  dplyr::distinct(c_label, species_family, .keep_all=T) %>%
  dplyr::select(species_family, fixed_approximate_number_of_worms, c_label) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(plot_class = ifelse(n() > 1, "multiple", species_family)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::group_by(plot_class, fixed_approximate_number_of_worms) %>%
  dplyr::mutate(count_per_size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_approximate_number_of_worms) %>%
  dplyr::mutate(total_size = n(), perc_size = count_per_size / total_size * 100) %>%
  dplyr::arrange(total_size) %>%
  #dplyr::select(fixed_substrate, plot_class, total_substrates, perc_worm_sub, worm_per_substrate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(plot_class, fixed_approximate_number_of_worms, .keep_all = T) %>%
  dplyr::mutate(plot_class = factor(plot_class, levels = c("multiple", "Other PCR +", "C. briggsae", "C. tropicalis", "C. sp. 53", "C. elegans")))

# Plotting estimated population size of Caenorhabditis collections by species and by substrate type
caen_species_pop_size_a <- ggplot(data = pop_size_caen_by_sp) +
  geom_bar(stat = "identity", aes(x = factor(fixed_approximate_number_of_worms), y = perc_size, fill = plot_class), colour = "black") + 
  scale_fill_manual(values =c(species_palette)) +
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black")) +
  labs(fill = "", x = "Approx. # of nematodes", y = "Percentage of Caenorhabditis\ncollections")+
  geom_text(aes(x=fixed_approximate_number_of_worms, y= 104, label=paste0("n=",total_size)), 
            position = position_dodge(width=1), size = 2.5, color = "black") +
  guides(fill = guide_legend(reverse = T)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 104))

caen_species_pop_size_b <- ggplot(data = pop_size_caen_by_sub) +
  geom_bar(stat = "identity", aes(x = factor(fixed_approximate_number_of_worms), y = perc_worm_popsize, fill = fixed_substrate), colour = "black") + 
  scale_fill_manual(values =c(substrate_palette)) +
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black")) +
  labs(fill = "", x = "Approx. # of nematodes", y = "Percentage of Caenorhabditis\ncollections")+
  geom_text(aes(x=fixed_approximate_number_of_worms, y= 104, label=paste0("n=",total_popsize)), 
            position = position_dodge(width=1), size = 2.5, color = "black") +
  guides(fill = guide_legend(reverse = F)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 104))

####### Figure 4 C and D ----------------------------------------------
# Define po regions
po_region = list(
  "Hawaii" = c(-160.7605,17.3977,-153.7671,24.3923),
  "NU" = c(-89.0033,41.5047,-86.9942,42.5936))
  
# label evanston po team
evanston_collectors <- c("joost.vanderzwaag@wur.nl",
                           "steffen.hahnel@northwestern.edu",
                           "robyn.tanny@northwestern.edu",
                           "tcrombie@northwestern.edu",
                           NA)

# Filter cso to c-labels with appropriate dauer info
dauer_cso <-  cso %>%
  dplyr::mutate(isolation_location = ifelse(po_created_by %in% evanston_collectors, "Evanston", "Hawaii")) %>%
  dplyr::filter(isolation_location == "Hawaii")

dauer_fractions <- dauer_cso %>%
  dplyr::filter(spp_id %in% c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae")) %>%
  dplyr::mutate(species_family = spp_id) %>%
  dplyr::mutate(dauer = ifelse(dauers_on_sample == "yes", 1, 0)) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass",
                                                                               "Isopod",
                                                                               "Slug",
                                                                               "Millipede"), 
                                                              "Other",substrate))))) %>%
  dplyr::distinct(c_label, species_family, .keep_all=T) %>%
  dplyr::group_by(fixed_substrate) %>%
  dplyr::mutate(total_collections_per_sub = n(),
                dauers_per_sub = sum(dauer)) %>%
  dplyr::group_by(species_family) %>%
  dplyr::mutate(total_collections_per_sp = n(),
                dauers_per_sp = sum(dauer)) %>%
  dplyr::select(spp_id, c_label,  species_family, dauer, fixed_substrate,
                total_collections_per_sub, dauers_per_sub,  total_collections_per_sp, dauers_per_sp) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dauer_fraction_sub = dauers_per_sub/total_collections_per_sub * 100) %>%
  dplyr::mutate(dauer_fraction_sp = dauers_per_sp/total_collections_per_sp * 100) %>%
  dplyr::distinct(fixed_substrate, species_family, .keep_all = T) %>%
  dplyr::mutate(species_family = factor(species_family, levels = c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae")),
                fixed_substrate = factor(fixed_substrate, levels = c("Leaf litter", "Fruit/nut/veg", "Flower", "Fungus", "Other")))

caen_species_dauer_subs <- ggplot(data = dauer_fractions %>% dplyr::distinct(fixed_substrate, dauer_fraction_sub, .keep_all =T)) +
  geom_bar(stat = "identity", aes(x = factor(fixed_substrate), y = dauer_fraction_sub), colour = "black", fill = "light grey") + 
  #scale_fill_manual(values =c(substrate_palette)) +
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "", x = "", y = "Percentage of Caenorhabditis\ncollections with dauers")+
  geom_text(aes(x=fixed_substrate, y= 104, label=paste0("n=",total_collections_per_sub)), 
            position = position_dodge(width=1), size = 2.5, color = "black") +
  guides(fill = guide_legend(reverse = F)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 104))

caen_species_dauer_sp <- ggplot(data = dauer_fractions %>% dplyr::distinct(species_family, dauer_fraction_sp, .keep_all =T)) +
  geom_bar(stat = "identity", aes(x = factor(species_family), y = dauer_fraction_sp), colour = "black", fill = "light grey") + 
  #scale_fill_manual(values =c(substrate_palette)) +
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "", x = "", y = "Percentage of collections\nwith dauers")+
  geom_text(aes(x=species_family, y= 104, label=paste0("n=",total_collections_per_sp)), 
            position = position_dodge(width=1), size = 2.5, color = "black") +
  guides(fill = guide_legend(reverse = F)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 104))


full <- cowplot::plot_grid(caen_species_pop_size_b,
                           caen_species_dauer_subs,
                           caen_species_pop_size_a, 
                           caen_species_dauer_sp,
                           labels = c("A", "C", "B", "D"),
                           ncol = 2, 
                           align = "hv",
                           axis = "tb",
                           rel_widths = c(1.5,1))
full

ggsave('plots/Figure4.pdf', width = 7.5, height = 5)

####### Figure 4B without multiples ----------------------------------------
pop_size_caen_by_sp_no_multiples  <- cso %>%
  dplyr::filter(spp_id %in% c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae")) %>%
  dplyr::mutate(species_family = spp_id) %>%
  dplyr::mutate(species_family = factor(species_family, levels = c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae"))) %>%
  dplyr::mutate(fixed_approximate_number_of_worms = ifelse(approximate_number_of_worms == "Very Few (1-3)", "1-3",
                                                           ifelse(approximate_number_of_worms == "Few (4-10)", "4-10",
                                                                  ifelse(approximate_number_of_worms == "Some (11-25)", "11-25",
                                                                         ifelse(approximate_number_of_worms == "Proliferating (25+)", "25+",
                                                                                approximate_number_of_worms))))) %>%
  dplyr::mutate(fixed_approximate_number_of_worms = factor(fixed_approximate_number_of_worms, levels = c("1-3", "4-10", "11-25", "25+"))) %>%
  dplyr::distinct(c_label, .keep_all=T) %>%
  dplyr::select(species_family, fixed_approximate_number_of_worms, c_label) %>%
  dplyr::group_by(species_family, fixed_approximate_number_of_worms) %>%
  dplyr::mutate(count_per_size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_approximate_number_of_worms) %>%
  dplyr::mutate(total_size = n(), perc_size = count_per_size / total_size * 100) %>%
  dplyr::arrange(total_size) %>%
  #dplyr::select(fixed_substrate, plot_class, total_substrates, perc_worm_sub, worm_per_substrate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(species_family, fixed_approximate_number_of_worms, .keep_all = T) %>%
  dplyr::mutate(species_family = factor(species_family, levels = c("C. briggsae", "C. tropicalis", "C. sp. 53", "C. elegans")))


caen_species_pop_size_b_no_multiples <- ggplot(data = pop_size_caen_by_sp_no_multiples) +
  geom_bar(stat = "identity", aes(x = factor(fixed_approximate_number_of_worms), y = perc_size, fill = species_family), colour = "black") + 
  scale_fill_manual(values =c(species_palette)) +
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black")) +
  labs(fill = "", x = "Approx. # of nematodes", y = "Percentage of Caenorhabditis\ncollections")+
  geom_text(aes(x=fixed_approximate_number_of_worms, y= 104, label=paste0("n=",total_size)), 
            position = position_dodge(width=1), size = 2.5, color = "black") +
  guides(fill = guide_legend(reverse = T)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 104))

full_no_multiples <- cowplot::plot_grid(caen_species_pop_size_b,
                           caen_species_dauer_subs,
                           caen_species_pop_size_b_no_multiples, 
                           caen_species_dauer_sp,
                           labels = c("A", "C", "B", "D"),
                           ncol = 2, 
                           align = "hv",
                           axis = "tb",
                           rel_widths = c(1.5,1))
full_no_multiples

ggsave('plots/Figure4_no_multiples.pdf', width = 7.5, height = 5)


#### Spearman's rho statistic is used to estimate a rank-based measure of association between substrate and population size
cor_test_df <- cso %>%
  dplyr::filter(spp_id %in% c("C. elegans", "C. tropicalis", "C. briggsae",  "C. sp. 53")) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass",
                                                                               "Isopod",
                                                                               "Slug",
                                                                               "Millipede"), 
                                                              "Other",substrate))))) %>%
  dplyr::distinct(c_label, spp_id, .keep_all = T) %>%
  dplyr::select(spp_id, fixed_substrate, approximate_number_of_worms) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = c("Leaf litter", "Flower", "Fruit/nut/veg", "Fungus",  "Other")),
                fixed_substrate1 = as.numeric(fixed_substrate),
                approx_pop = factor(approximate_number_of_worms, levels = c("Very Few (1-3)", "Few (4-10)", "Some (11-25)", "Proliferating (25+)")),
                approx_pop1 = as.numeric(approx_pop),
                approx_pop2 = ifelse(approximate_number_of_worms == "Proliferating (25+)", 2,1))

# setup test levels for cb and ct and do test inside of dataframe
cor_test_df_cb_ct_test <- cor_test_df %>%
  dplyr::filter(fixed_substrate %in% c("Flower", "Leaf litter")) %>%
  dplyr::group_by(spp_id) %>%
  dplyr::do(broom::tidy(cor.test(.$fixed_substrate1, .$approx_pop2, alternative = "g")))

# test levels for cb and ct for flower vs leaf litter plot
cor_test_df_cb_ct <- cor_test_df %>%
  dplyr::filter(fixed_substrate %in% c("Flower", "Leaf litter")) %>%
  dplyr::group_by(spp_id)

# setup test levels for co and do test inside of dataframe
cor_test_df_co_test <- cor_test_df %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = c("Fruit/nut/veg", "Flower", "Leaf litter", "Fungus",  "Other")), # reorder for cor test between fruit and flower 
                fixed_substrate1 = as.numeric(fixed_substrate)) %>%
  dplyr::filter(fixed_substrate %in% c("Fruit/nut/veg", "Flower"))%>%
  dplyr::group_by(spp_id) %>%
  dplyr::do(broom::tidy(cor.test(.$fixed_substrate1, .$approx_pop2, alternative = "g")))

# test levels for co for flower vs Fruit plot
  cor_test_df_co <- cor_test_df %>%
    dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = c("Fruit/nut/veg", "Flower", "Leaf litter", "Fungus",  "Other")), # reorder for cor test between fruit and flower 
                  fixed_substrate1 = as.numeric(fixed_substrate)) %>%
    dplyr::filter(fixed_substrate %in% c("Fruit/nut/veg", "Flower"))
  
# plot population size vs substrate scatter 
plot_co_test <- ggplot(cor_test_df_co %>% dplyr::filter(spp_id == "C. sp. 53")) +
  aes(x = fixed_substrate1, y = approx_pop2, fill = fixed_substrate) +
  geom_jitter(width = 0.025, height = 0.025, shape = 21, size = 3) +
  geom_smooth(method = lm, inherit.aes = F, aes(x = fixed_substrate1, y = approx_pop2))
plot_co_test


######################
## count multiples
multiples <- cso %>%
  dplyr::filter(spp_id %in% c("C. elegans", "C. tropicalis", "C. briggsae",  "C. sp. 53")) %>%
  dplyr::distinct(c_label, spp_id, .keep_all = T) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(multiples = ifelse(n() > 1, TRUE, FALSE)) %>%
  dplyr::filter(multiples ==T) %>%
  dplyr::select(spp_id, substrate, c_label)
