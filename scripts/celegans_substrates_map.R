library(tidyverse)
library(ggmap)
library(memoise)
library(lubridate)
library(cowplot)
library(pals)
library(grid)
library(gridExtra)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# identify collection locations and type strain
ce <- cso %>%
  dplyr::filter(spp_id == "C. elegans") %>%
  dplyr::mutate(num_of_ce_isolates = n(),
                num_of_ce_collections = length(unique(c_label))) %>%
  dplyr::select(strain, num_of_ce_isolates, num_of_ce_collections, isotype, c_label, s_label, strain, island, substrate, latitude, longitude, spp_id, photo_url_thumb)

# Collection plots with intersting images. 
ce_sub_images <- ce %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::arrange(substrate) %>%
  dplyr::mutate(plot_labels = paste0(c_label,"\n",substrate,"\n",spp_id)) %>%
  dplyr::select(c_label, substrate, spp_id, photo_url_thumb, plot_labels)

# create function to load and rasterize images from a set of urls in the cso or df
rplots <- function(x) {
  ggplot() + annotation_raster(raster = imager::load.image(x), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
}

# generate a list of plots using rplots function
pl <- lapply(ce_sub_images$photo_url_thumb, rplots)

# plot all plots in the pl using grid.arrange
gl <- gridExtra::grid.arrange(grobs=pl)

# Overlay grid with labels working but not best way to do it. font size should be very small to plot all
ids <- function(id){
  ggplot() + ggtitle(id) + theme(plot.title = element_text(size = 8, face = "bold", color = "white"))
}

id_plots <- lapply(ce_sub_images$plot_labels, ids)
id_print <- gridExtra::grid.arrange(grobs=id_plots)

# ploting images and labels. this takes a lot of time, be patient. maybe 10-15 minutes
ggplot()
grid.draw(gl)
grid.draw(id_print)

#ggsave('plots/coiwi_substrate_images.pdf', height = 8.24, width = 5)
