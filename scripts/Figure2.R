library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(imager)
library(jpeg)

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
                     "No Worm" = "#222222",
                     "multiple" = "#FFFFFF") #2

substrate_palette <- c("Leaf litter" = "#E68FAC",
                       "Fruit/nut/veg" = "#0067A5",
                       "Flower" = "#DCD300",
                       "Fungus" = "#604E97",
                       "Compost" = "#F6A600",
                       "Other" = "#B3446C")

allPalette <- c("#222222", "#b3b3b3","#848482", "#008856")
names(allPalette) = c("No nematode", "Tracks only", "PCR -", "PCR +")

#######################################
# Figure 2A                           #
#######################################
#Figure 2A define rh postitive c-labels
rh_positive_c_labels <- cso %>%
  dplyr::filter(pcr_rhpositive == 1) %>%
  dplyr::distinct(c_label) 

#Figure 2A df for all collections broken into braod collection categories (i.e., No worm, tracks, pcr-, pcr+)
worms <- df %>%
  dplyr::filter(worms_on_sample != "?") %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass"), 
                                                              "Other",substrate))))) %>%
  dplyr::mutate(worms_on_sample = ifelse(worms_on_sample == "No", "No nematode",
                                         ifelse(worms_on_sample == "Tracks", "Tracks only",
                                                ifelse(worms_on_sample == "Yes", "PCR -", worms_on_sample)))) %>%
  dplyr::mutate(worms_on_sample = ifelse(c_label %in% rh_positive_c_labels$c_label, "PCR +", worms_on_sample)) %>%
  dplyr::distinct(c_label, .keep_all = TRUE) %>%
  dplyr::group_by(worms_on_sample, fixed_substrate) %>%
  dplyr::mutate(worm_per_substrate = n()) %>% 
  dplyr::distinct(worms_on_sample, fixed_substrate, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_substrate) %>%
  dplyr::mutate(total_substrates = sum(worm_per_substrate), perc_worm_sub = worm_per_substrate / total_substrates * 100) %>%
  dplyr::arrange(total_substrates) %>%
  dplyr::select(fixed_substrate, worm_per_substrate, total_substrates, perc_worm_sub, worms_on_sample) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(worms_on_sample = factor(worms_on_sample, levels = c("No nematode", "Tracks only", "PCR -", "PCR +"))) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = rev(c("Leaf litter","Fruit/nut/veg", "Flower", "Fungus",
                                                                         "Isopod",  "Millipede", "Slug", "Other"))))

# plot for all collections
Fig2A <- ggplot(worms) +
  geom_bar(stat = "identity", aes(x = factor(fixed_substrate), y = perc_worm_sub, fill = worms_on_sample), colour = "black") + 
  scale_fill_manual(values=c(allPalette))+
  coord_flip() + 
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  labs(fill = "", x = "", y = "Percentage of all collections") +
  geom_text(aes(x=fixed_substrate, y=113, label=paste0("n=", total_substrates)), 
            position = position_dodge(width=1), size = 2.5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 120))

#######################################
# Figure 2B                           #
#######################################
newdf <- cso %>%
  dplyr::filter(pcr_rhpositive == 1) %>%
  dplyr::mutate(species_family = ifelse(spp_id %in% c("Oscheius sp.",
                                                      "Panagrolaimus sp.",
                                                      "Chabertia ovina", 
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
                                        "Other PCR +", spp_id)) %>%
  dplyr::select(c_label, spp_id, species_family, substrate) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass"), 
                                                              "Other",substrate))))) %>%
  dplyr::distinct(c_label, species_family, .keep_all=T) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(plot_class = ifelse(n() > 1, "multiple", species_family)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::group_by(plot_class, fixed_substrate) %>%
  dplyr::mutate(worm_per_substrate = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_substrate) %>%
  dplyr::mutate(total_substrates = n(), perc_worm_sub = worm_per_substrate / total_substrates * 100) %>%
  dplyr::arrange(total_substrates) %>%
  #dplyr::select(fixed_substrate, plot_class, total_substrates, perc_worm_sub, worm_per_substrate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(plot_class, fixed_substrate, .keep_all = T) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = rev(names(substrate_palette)))) %>%
  dplyr::mutate(plot_class = factor(plot_class, levels = c("multiple", "Other PCR +", "C. briggsae", "C. tropicalis", "C. sp. 53", "C. elegans")))

# Fig2B plot for rhabditida positive collections
Fig2B <- ggplot(data = newdf) +
  geom_bar(stat = "identity", aes(x = factor(fixed_substrate), y = perc_worm_sub, fill = plot_class), colour = "black") +
  scale_fill_manual(values=c(species_palette)) +
  coord_flip() + 
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  labs(fill = "", x = "", y = "Percentage of PCR-positive collections") +
  geom_text(aes(x=fixed_substrate, y=109, label=paste0("n=",total_substrates)), 
            position = position_dodge(width=1), size = 2.5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 120))

# Identify multiples for plotting
multiples <- cso %>%
  dplyr::filter(pcr_rhpositive == 1) %>%
  dplyr::mutate(species_family = ifelse(spp_id %in% c("Oscheius sp.",
                                                      "Panagrolaimus sp.",
                                                      "Chabertia ovina", 
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
                                        "Other PCR +", spp_id)) %>%
  dplyr::select(c_label, spp_id, species_family, substrate) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass"), 
                                                              "Other",substrate))))) %>%
  dplyr::distinct(c_label, species_family, .keep_all=T) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(plot_class = ifelse(n() > 1, "multiple", species_family)) %>%
  dplyr::filter(plot_class == "multiple") %>%
  dplyr::arrange(species_family) %>%
  dplyr::select(fixed_substrate, species_family, c_label, plot_class) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(label = paste0(species_family, collapse = ", ")) %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::arrange(fixed_substrate, label) %>%
  dplyr::select(fixed_substrate, label, c_label, plot_class)
  
##################################
# Plot 2A and 2B together        #
##################################
# Put fig 2A and fig 2B together
panel_substrate <- cowplot::plot_grid(Fig2A, Fig2B, labels = c("A","B"), ncol=1, align = "v")
ggsave('plots/Fig2AB_substrate_panel.pdf', height = 6.59, width = 5)

##################################
# C Image gallery                #
##################################
# create function to load and rasterize images from a set of urls in the cso or df
rplots <- function(x) {
  ggplot() + annotation_raster(raster = imager::load.image(x), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
}

# Collection plots with intersting images. 
mdf <- cso %>%
  dplyr::filter(pcr_rhpositive == 1) %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::arrange(spp_id, substrate) %>%
  dplyr::mutate(., plot_labels = paste0(c_label,"\n",substrate,"\n",spp_id)) %>%
  dplyr::select(c_label, substrate, spp_id, photo_url_thumb, plot_labels )

spp <- rep(c("C. elegans", "C. briggsae", "C. tropicalis", "C. sp. 53"), each = 3)
type <- rep(c("Leaf litter", "Flower", "Fruit/nut/veg"), 4)
c_label <- c("C-3133", "C-0736", "C-2830", NA,"C-0257", "C-0072",
             "C-1083", "C-2909", "C-2906", "C-2735", "C-0846", "C-2845")


# Make df with dummy grey box image for plotting blanks
sub_df <- as.data.frame(cbind(spp, type, c_label)) %>%
  dplyr::left_join(., mdf) %>%
  dplyr::mutate(photo_url_thumb = ifelse(is.na(photo_url_thumb),"https://drive.google.com/uc?export=download&id=10D7zXBStCPcdHHk_UE2r6QBZYqHiBNL6", photo_url_thumb))

pl <- lapply(sub_df$photo_url_thumb, rplots)
margin = theme(plot.margin = unit(c(.5,.5,.5,.5), "mm"))
gl <- gridExtra::grid.arrange(grobs = lapply(pl, "+", margin), ncol = 3)

z <- textGrob("")
a <- textGrob("C. elegans")
b <- textGrob("C. briggsae")
c <- textGrob("C. tropicalis")
d<- textGrob("C. sp. 53")

sa <- textGrob("Leaf litter")
sb <- textGrob("Flower")
sc <- textGrob("Fruit/nut/veg")

yaxis <- grid.arrange(z,a,d,c,b, ncol = 1)
xaxis <- grid.arrange(sa,sb,sc, nrow = 1)
lay <- rbind(c(1,2,2,2),
             c(1,3,3,3),
             c(1,3,3,3),
             c(1,3,3,3),
             c(1,3,3,3))

full <- grid.arrange(yaxis, xaxis, gl, layout_matrix = lay, ncol = 4, nrow =5)

plot_save <- cowplot::plot_grid(full)
plot_save
# the height and width ratio accounts for image dimesions and column numbers to avoid transforming image dimensions. Height = 1.648 times Width
ggsave('plots/20190219_Reduced_Substrate_images_test.pdf', height = 8.24, width = 5)


#######################################
# Figure 2B  with out multiples       #
#######################################
newdf_no_multiples <- cso %>%
  dplyr::filter(pcr_rhpositive == 1) %>%
  dplyr::mutate(species_family = ifelse(spp_id %in% c("Oscheius sp.",
                                                      "Panagrolaimus sp.",
                                                      "Chabertia ovina", 
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
                                        "Other PCR +", spp_id)) %>%
  dplyr::select(c_label, spp_id, species_family, substrate) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit/nut/veg",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass"), 
                                                              "Other",substrate))))) %>%
  dplyr::mutate(species_family = factor(species_family, levels = c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae","Other PCR +"))) %>%
  dplyr::arrange(species_family) %>%
  dplyr::distinct(c_label, .keep_all=T) %>%
  dplyr::group_by(species_family, fixed_substrate) %>%
  dplyr::mutate(worm_per_substrate = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_substrate) %>%
  dplyr::mutate(total_substrates = n(),
                perc_worm_sub = worm_per_substrate / total_substrates * 100) %>%
  dplyr::arrange(total_substrates) %>%
  #dplyr::select(fixed_substrate, plot_class, total_substrates, perc_worm_sub, worm_per_substrate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(species_family, fixed_substrate, .keep_all = T) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = rev(names(substrate_palette)))) %>%
  dplyr::mutate(species_family = factor(species_family, levels = rev(c("C. elegans", "C. sp. 53", "C. tropicalis", "C. briggsae","Other PCR +"))))
  
  

# Fig2B plot for rhabditida positive collections
Fig2B_no_multiples <- ggplot(data = newdf_no_multiples) +
  geom_bar(stat = "identity", aes(x = factor(fixed_substrate), y = perc_worm_sub, fill = species_family), colour = "black") +
  scale_fill_manual(values=c(species_palette)) +
  coord_flip() + 
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  labs(fill = "", x = "", y = "Percentage of PCR-positive collections") +
  geom_text(aes(x=fixed_substrate, y=109, label=paste0("n=",total_substrates)), 
            position = position_dodge(width=1), size = 2.5) +
  guides(fill = guide_legend(reverse = T)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 120))


panel_substrate_no_multiples <- cowplot::plot_grid(Fig2A, Fig2B_no_multiples, labels = c("A","B"), ncol=1, align = "v")
ggsave('plots/Fig2AB_substrate_panel_no_multiples.pdf', height = 6.59, width = 5)
