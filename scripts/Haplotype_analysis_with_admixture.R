library(tidyverse)
library(ggmap)
library(memoise)
library(lubridate)
library(cowplot)
library(pals)
library(grid)
library(gridExtra)
library(DT)
library(FSA

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))
        
# load data
load('data/fulcrum/df.Rda')
                
# load Q files
base.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/Manuscripts/Hawaii_Manuscript/"
setwd(glue::glue("{base.dir}Data/ADMIXTURE/249_Hawaii/BEST_K/"))
qfile_name <- grep(pattern = paste0("Q$"), value = T, x = list.files())
qfile <- pophelper::readQ(files = qfile_name)

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

#####################################################################
# Confirm most common haplotype in admix pop B is the global sweep  #
#####################################################################
#load haplotypes
haplotypes <- read.table("~/Dropbox/AndersenLab/Hawaii_manuscript/data/HAPLOTYPE/haplotype.tsv", sep = "\t", header = TRUE)

