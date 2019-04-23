library(janitor)
library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

#=========================#
# Generate Summary Tables #
#=========================#

df %>%
  dplyr::select(isotype,
                strain,
                c_label,
                s_label,
                worms_on_sample,
                approximate_number_of_worms,
                males_observed,
                latitude,
                longitude,
                geometry,
                substrate,
                substrate_other,
                substrate_notes,
                landscape,
                sky_view,
                gridsect,
                gridsect_direction,
                gridsect_radius,
                substrate_temperature,
                substrate_moisture,
                substrate_moisture_issue,
                datetime,
                date) %>% readr::write_tsv("table/collection_summary.tsv")


#===================#
# Species Collected #
#===================#

# Species by island - copied from 20180122-map-overview.R
table_cso <- cso %>%
  dplyr::mutate(plot_type = ifelse(worms_on_sample %in% c("No","?"), "No Worm",
                                   ifelse(is.na(spp_id), "Not genotyped",
                                          ifelse(pcr_rhpositive == 0, "PCR -",
                                                 ifelse(spp_id == "Unknown", "Not genotyped",
                                                        spp_id))))) %>%
  dplyr::distinct(c_label, plot_type, .keep_all = T) %>% 
  dplyr::mutate(plot_type = forcats::as_factor(plot_type),
                plot_type = forcats::fct_relevel(plot_type,
                                                 "C. elegans",
                                                 "C. sp. 53",
                                                 "C. tropicalis",
                                                 "C. kamaaina",
                                                 "C. briggsae",
                                                 "Panagrolaimus sp.",
                                                 "Oscheius sp.",
                                                 "Teratorhabditis sp.",
                                                 "Rhabditis terricola",
                                                 "Choriorhabditis sp.",
                                                 "Mesorhabditis sp.",
                                                 "Chabertia ovina",
                                                 "Heterhabditis zealandica",
                                                 "PCR -",
                                                 "Not genotyped",
                                                 "No Worm")) %>%
  dplyr::arrange(desc(plot_type)) %>%
  dplyr::rename(collection_category = plot_type) %>%
  dplyr::group_by(collection_category, island) %>%
  dplyr::summarize(worm_isolates=as.integer(n())) %>%
  tidyr::spread(island, worm_isolates, fill=0) %>%
  janitor::adorn_totals('row') %>%
  janitor::adorn_totals('col') %>% readr::write_tsv("table/collection_category_identified_by_island.tsv")




#==================================#
# Number of worm isolates / sample #
#==================================#

cso %>%
  dplyr::group_by(c_label) %>%
  dplyr::filter(!is.na(s_label)) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::summarize(m = mean(n))

# Species by sample

cso %>% 
  dplyr::group_by(c_label) %>%
  dplyr::filter(!is.na(spp_id)) %>%
  dplyr::distinct(c_label, spp_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n >= 2) %>%
  tidyr::nest(spp_id, .key = "species") %>%
  dplyr::mutate(has_unknown = purrr::map_lgl(species,  ~ "Unknown" %in% .x$spp_id)) %>%
  dplyr::mutate(n_species = purrr::map_int(species, ~ length(.x$spp_id))) %>% 
  dplyr::select(has_unknown) %>% table()


#=======================#
# Isotypes by substrate #
#=======================#

cso %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(isotype, substrate) %>%
  dplyr::group_by(substrate) %>%
  dplyr::summarize(n=n())

#=======================#
# Isotypes by landscape #
#=======================#

cso %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(isotype, landscape) %>%
  dplyr::group_by(landscape) %>%
  dplyr::summarize(n=n())

