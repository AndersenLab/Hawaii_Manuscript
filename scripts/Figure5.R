# load necessary packages
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

#====================#
# Gridsect Functions #
#====================#
isotype_palette <- c("ECA760" = "#654522" ,
                     "ECA778" =  "#8DB600",
                     "ECA768" = "#882D17" , 
                     "ECA812" = "#DCD300" ,
                     "ECA730" = "#B3446C" , 
                     "ECA712" =  "#F6A600", 
                     "ECA777" = "#604E97" , 
                     "ECA807" = "#F99379")

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

substrate_shapes <- c("Leaf litter" = 21,
                      "Fruit/nut/vegetable"= 24,
                      "Rotting flower"= 22,
                      "Fungus"= 23,
                      "Isopod"= 25)

adjust_x <- function(x, n, rn) {
  
  if (n == 1) {
    return(x)
  } else if ((n == 2 && rn == 1) || (n == 3 && rn == 2) || (n == 4 && rn == 1)  || (n == 4 && rn == 3) || (n %in% c(6, 7) && rn %in% c(3,5)) ) {
    return(x - 0.09)
  } else if (n == 2 && rn == 2  || (n == 3 && rn == 3) || (n == 4 && rn == 2) || (n == 4 && rn == 4) || (n %in% c(6, 7) && rn %in% c(2,4)) ) {
    return(x + 0.09)
  } else if (n == 3 && rn >= 2) {
    return(x + 0.09)
  } else if (n == 5 && rn == 2) {
    return(x - 0.09)
  } else if (n == 5 && rn == 3) {
    return(x + 0.09)
  } else if (n == 5 && rn == 4) {
    return (x - 0.09)
  } else if (n == 5 && rn == 5) {
    return (x + 0.09)
  } else if (n == 7 && rn == 7) {
    return (x + 0.25)
  }
  return(x)
}

adjust_y <- function(y, n, rn) {
  if (n == 1) {
    return(y)
  } else if ( (n == 3 && rn == 1 ) || (n == 4 && rn <= 2) || (n == 5 && rn %in% c(4, 5)) || (n %in% c(6, 7) && rn %in% c(4,5))) {
    return(y - 0.09)
  } else if ( (n == 3 && rn >= 2 ) || (n == 4 && rn >= 3) || (n == 5 && rn %in% c(2, 3))  || (n %in% c(6, 7) && rn %in% c(2,3)) ) {
    return (y + 0.09)
  } else if ((n == 5 && rn == 1) || (n == 6 && rn == 1)) {
    return (y + 0.25)
  } else if ( (n %in% c(6, 7) && rn == 6)) {
    return (y - 0.25)
  }
  return(y)
}

set_size <- function(n) {
  if (n == 1) {
    return(5)
  } else {
    return(2.4)
  }
}


plot_gridsect <- function(gn, cso = cso) {
  angles = list(A = 1,
                B = 2,
                C = 3,
                D = 4,
                E = 5,
                F = 6)
  
  mm <- cso %>%
    dplyr::filter(grid_num == gn, !is.na(grid_num)) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(x = gridsect_radius * (sin( (angles[[gridsect_direction]] - 1) * (pi/3) )),
                  y = gridsect_radius * (cos( (angles[[gridsect_direction]] - 1) * (pi/3) )),
                  label = paste0(gridsect_direction, gridsect_radius)) %>%
    dplyr::group_by(gridsect_radius, gridsect_direction) %>%
    dplyr::mutate(n = n(),
                  rn = row_number()) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(x = adjust_x(x, n, rn),
                  y = adjust_y(y, n, rn),
                  size_c = set_size(n)) %>%
    dplyr::ungroup()
  
  ggplot(mm, aes(x = x, y = y, label = label)) +
    annotate("path",
             x=0+3.5*cos(seq(0,2*pi,length.out=100)),
             y=0+3.5*sin(seq(0,2*pi,length.out=100))) +
    theme_void() +
    theme(legend.position = "none") +
    scale_size_area(guide = "none")
}

#=========#
# Species # 
#=========#

sp = cso %>%
  dplyr::filter(!is.na(grid_num)) %>%
  #dplyr::mutate(spp_id = ifelse(is.na(spp_id), "Not genotyped", spp_id)) %>%
  dplyr::mutate(plot_type = ifelse(worms_on_sample %in% c("No","?"), "No Worm",
                                   ifelse(is.na(spp_id), "Not genotyped",
                                          ifelse(pcr_rhpositive == 0, "PCR -",
                                                 ifelse(spp_id %in% c("Chabertia ovina",
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
                                                                      NA), "Other PCR +", spp_id)))))

gridsect_list <- lapply(c(1:15, 17:21), function(x){ 
  if (x == 21) {
    plot_gridsect(5, sp) + 
      geom_text(size= 8)
  } else {
    plot_gridsect(x, sp) +
      geom_point(aes(fill = plot_type,  size = size_c, shape = substrate), stroke = 0.3) +  # scale_shape_manual(values=c(3, 16, 17))+
      scale_fill_manual("Species", values = c(species_palette), drop = FALSE) +
      scale_shape_manual("Substrate", values = c(substrate_shapes), drop = FALSE)
  }
})
legend <- cowplot::get_legend(gridsect_list[[3]] + theme(legend.position = "bottom",
                                                          legend.text = element_text(size = 10),
                                                          legend.title = element_text(size = 12)))

plot_grid(plot_grid(plotlist=gridsect_list, hjust = -1, label_size = 20,
                    nrow = 4, ncol=5, rel_heights = c(5,5,5,5), labels=c(1:15, 17:20)), legend,
          rel_heights = c(20, 1), nrow=2)

cowplot::ggsave("plots/gridsect_species_19_shapes.pdf", width = 28.56, height = 24)
cowplot::ggsave("plots/gridsect_species_19_shapes.png", width = 28.56, height = 24)

# Gridsects 1 and 3 by species
fig_5a <- cowplot::plot_grid(gridsect_list[[1]],
          gridsect_list[[3]],
          gridsect_list [[20]],
          hjust = -0.09,
          labels = c("A", "B", "C"),
          cols = 3)

cowplot::ggsave(filename = "plots/Figure5a_substrates.pdf", plot = fig_5a, height = 3*57.5, width = 3*170, units = "mm")

# Gridsects 1 and 3 by isotypes
isotype = cso %>%
  dplyr::mutate(isotype = ifelse(is.na(isotype), "No isotype", isotype))

gridsect_list <- lapply(c(1:15, 17:21), function(x){ 
  if (x == 21) {
    plot_gridsect(5, isotype) + 
      geom_text(size= 8)
  } else {
    plot_gridsect(x, isotype) +
      geom_point(aes(fill = isotype,  size = size_c,  shape = substrate), stroke = 0.3) +
      scale_fill_manual("Isotypes", values = c(isotype_palette), drop = FALSE) +
      scale_shape_manual("Substrate", values = c(substrate_shapes), drop = FALSE)
  }
})
legend <- cowplot::get_legend(gridsect_list[[11]] + theme(legend.position = "bottom",
                                                          legend.text = element_text(size = 10),
                                                          legend.title = element_text(size = 12)))

fig_5b <- cowplot::plot_grid(gridsect_list[[1]],
                             gridsect_list[[3]],
                             gridsect_list [[20]],
                             hjust = -0.09,
                             labels = c("A", "B", "C"),
                             cols = 3)

cowplot::ggsave(filename = "plots/Figure5b_substrates.pdf", plot = fig_5b, height = 3*57.5, width = 3*170, units = "mm")

# make nice legends
legend_isotypes <- get_legend(ggplot(isotype %>% dplyr::filter(isotype %in% names(isotype_palette))) +
                       geom_bar(aes(x = strain, fill = isotype)) +
                       scale_fill_manual(values=c(isotype_palette)) +
                       labs(fill = "isotypes"))
legend_species <- get_legend(ggplot(sp %>% dplyr::filter(plot_type %in% names(species_palette)) %>%
                      dplyr::mutate(plot_type = factor(plot_type, levels = names(species_palette)))) +
                      geom_bar(aes(x = strain, fill = plot_type)) +
                      scale_fill_manual(values=c(species_palette)) +
                      labs(fill = "species"))

legend_substrates <- get_legend(ggplot(sp %>% dplyr::mutate(substrate = factor(substrate, levels = names(substrate_shapes)))) +
                      geom_point(aes(x = substrate_temperature, y = substrate_temperature, shape = substrate)) +
                      scale_shape_manual("substrates", values = c(substrate_shapes), drop = FALSE) +
                      labs(shape = "isotypes"))

plot_5ab_legends <- plot_grid(legend_isotypes,legend_species,legend_substrates)
ggsave(filename = "plots/Figure5ab_legends.pdf", plot = plot_5ab_legends, height = 170, width = 170, units = "mm")

# number of isotypes on a c_plate (2 collections are excluded because no isotypes were recovered from them do to extinction before sequencing)
num_isotypes_per_sample <- cso %>%
  dplyr::filter(!is.na(s_label)) %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(c_label, isotype, .keep_all = TRUE) %>%
  dplyr::select(isotype, spp_id, c_label, s_label, substrate, latitude, longitude) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(distinct_isotypes = n()) %>%
  dplyr::distinct(c_label, .keep_all=TRUE)

######################
# Describe gridsects and isotypes #
######################
# Find mean number of biologically distinct categories on a gridsect
grid_cat_numbers <- sp %>%
  dplyr::filter(grid_num != 16) %>%
  dplyr::filter(plot_type %in% c("PCR -", "C. briggsae", "Oscheius sp.", "C. elegans")) %>%
  dplyr::group_by(grid_num) %>%
  dplyr::mutate(num_species = length(unique(plot_type))) %>%
  dplyr::ungroup() %>%
  dplyr::add_row(grid_num = c(4,14), num_species = c(1,1)) %>% # correct for fact that we actually isolated a nematode from these gridsects we just never genotyped it
  dplyr::distinct(grid_num, .keep_all = T) %>%
  dplyr::mutate(mean_cat_num = mean(num_species))

# number of islands we did gridsects on
island_grids <- cso %>%
  dplyr::distinct(grid_num, island) %>%
  dplyr::filter(!is.na(grid_num))

# gridsect temp averages
grid_temps <- sp %>%
  dplyr::filter(grid_num != 16 & !is.na(ambient_temperature)) %>%
  dplyr::distinct(grid_num, c_label, .keep_all = T) %>%
  dplyr::group_by(grid_num) %>%
  dplyr::mutate(avg_temp = mean(substrate_temperature),
                min_temp = min(substrate_temperature),
                max_temp = max(substrate_temperature),
                avg_temp_amb = mean(ambient_temperature)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(grid_num, .keep_all = T) %>%
  dplyr::select(grid_num, avg_temp_amb, avg_temp, min_temp, max_temp)

# fraction of samples with multiple isotypes
iso_frac <- cso %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(c_label, isotype) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(iso_num = n()) %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(iso_num) %>%
  dplyr::mutate(iso_num_total = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac_iso = iso_num_total/n())




