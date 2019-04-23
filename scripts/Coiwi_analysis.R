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
coiwi <- cso %>%
  dplyr::filter(spp_id == "C. sp. 53") %>%
  dplyr::mutate(num_of_coiwi_isolates = n(),
                num_of_coiwi_collections = length(unique(c_label))) %>%
  dplyr::select(strain, num_of_coiwi_isolates, num_of_coiwi_collections, c_label, s_label, strain, island, substrate, latitude, longitude, spp_id, photo_url_thumb) %>%
  dplyr::mutate(strain = ifelse(s_label == "S-05114", "ECA1100", strain)) %>%
  dplyr::mutate(spp_id = ifelse(spp_id == "C. sp. 53", "C. oiwi", spp_id))
                  


# plot map of C. oiwi other map source
mapImageData2 <- get_map(center = c(lon = mean(coiwi$longitude), lat = mean(coiwi$latitude)),
                        maptype = "terrain-background",
                        source = "stamen",      
                        zoom = 7,
                        scale = "auto")
ggmap(mapImageData2, extent = "device") + # removes axes, etc.
  geom_point(aes(x = longitude,
                 y = latitude),
             data = coiwi,
             fill = "#875692",
             size = 2,
             shape = 21,
             colour = "black") +
  #coord_fixed(1) + 
  #coord_map(projection = "mercator",
            #xlim = c(-153.7671, -160.7605),
            #ylim = c(17.3977, 24.3923)) + 
  geom_text(data = coiwi, aes(x = longitude, y = latitude, label = strain), 
            size = 3, vjust = 0, hjust = -0.1)

ggsave('plots/coiwi_map_2.png', width = 7.5, height = 5)
ggsave('plots/coiwi_map_2.pdf', width = 7.5, height = 7.5)

# show collection images for C. oiwi

##################################
# C Image gallery C. oiwi        #
##################################
# create function to load and rasterize images from a set of urls in the cso or df
rplots <- function(x) {
  ggplot() + annotation_raster(raster = imager::load.image(x), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
}

# Collection plots with intersting images. 
coiwi_sub_images <- coiwi %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::arrange(substrate) %>%
  dplyr::mutate(plot_labels = paste0(c_label,"\n",substrate,"\n",spp_id)) %>%
  dplyr::select(c_label, substrate, spp_id, photo_url_thumb, plot_labels)

# generate a list of plots using rplots function
pl <- lapply(coiwi_sub_images$photo_url_thumb, rplots)

# plot all plots in the pl using grid.arrange
gl <- gridExtra::grid.arrange(grobs=pl)

# Overlay grid with labels working but not best way to do it. font size should be very small to plot all
ids <- function(id){
  ggplot() + ggtitle(id) + theme(plot.title = element_text(size = 8, face = "bold", color = "white"))
}

id_plots <- lapply(coiwi_sub_images$plot_labels, ids)
id_print <- gridExtra::grid.arrange(grobs=id_plots)

# ploting images and labels. this takes a lot of time, be patient. maybe 10-15 minutes
ggplot()
grid.draw(gl)
grid.draw(id_print)

#ggsave('plots/coiwi_substrate_images.pdf', height = 8.24, width = 5)

# # plot map of C. oiwi using google maps
# need to update ggmap
#devtools::install_github("dkahle/ggmap")
#######################register_google(key = "AIzaSyCYANDExSNdNMQAhwJetJXZCCxwVZxrkf8")######################

# mapImageData <- get_googlemap(center = c(lon = mean(coiwi$longitude), lat = mean(coiwi$latitude)),
#                               zoom = 7,
#                               scale = 2,
#                               maptype = "terrain",
#                               format = "png8")
# ggmap(mapImageData, extent = "device") + # removes axes, etc.
#   geom_point(aes(x = longitude,
#                  y = latitude),
#              data = coiwi,
#              fill = "#875692",
#              size = 2,
#              shape = 21,
#              colour = "black") +
#   coord_fixed(xlim = c(-154, -161), ylim = c(18.5, 22.5)) +
#   geom_text(data = coiwi, aes(x = longitude, y = latitude, label = strain), 
#             size = 3, vjust = 0, hjust = -0.1)
# 
# ggsave('plots/coiwi_map.png', width = 7.5, height = 5)