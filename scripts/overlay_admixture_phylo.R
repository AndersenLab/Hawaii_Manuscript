#!/usr/bin/env Rscript
library(tidyverse)
library(ggtree)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load in wi strain info for lat long coordinates for map
df <- data.table::fread("data/WI_strain_list.csv")
strain_colors <- dplyr::filter(df, state== "Hawaii")

tree_275 <- ape::read.tree(glue::glue("data/tree/249-hawaii_genome.raxml.bestTree"))

kpops <- 3
K <- as.numeric(strsplit(grep(".Q", list.files("data/ADMIXTURE/BEST_K/"), value = T)[kpops], split = "\\.")[[1]][4])

# load Q files
qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files("data/ADMIXTURE/BEST_K/"))
qfile <- pophelper::readQ(files = paste0("data/ADMIXTURE/BEST_K/",qfile_name))[[1]]
# add pop names
colnames(qfile) <- LETTERS[1:K]
sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples.txt", header = F) %>% dplyr::pull(V1))
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
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

# make long and determin order of plotting
long_admix_pops <- qfile %>%
  dplyr::mutate(samples = sample_names) %>%
  tidyr::gather(cluster, frac_cluster, -samples) %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(max_frac = max(frac_cluster)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))%>%
  dplyr::filter(max_frac == frac_cluster)

branch_strains <- list(A = as.character(dplyr::filter(long_admix_pops, cluster == "A", frac_cluster > 0.9) %>% dplyr::pull(samples)),
                       B = as.character(dplyr::filter(long_admix_pops, cluster == "B", frac_cluster > 0.9) %>% dplyr::pull(samples)),
                       C = as.character(dplyr::filter(long_admix_pops, cluster == "C", frac_cluster > 0.9) %>% dplyr::pull(samples)),
                       D = as.character(dplyr::filter(long_admix_pops, cluster == "D", frac_cluster > 0.9) %>% dplyr::pull(samples)),
                       E = as.character(dplyr::filter(long_admix_pops, cluster == "E", frac_cluster > 0.9) %>% dplyr::pull(samples)),
                       F = as.character(dplyr::filter(long_admix_pops, cluster == "F", frac_cluster > 0.9) %>% dplyr::pull(samples)),
                       G = as.character(dplyr::filter(long_admix_pops, frac_cluster < 0.9) %>% dplyr::pull(samples)))

tree_pt_h <- ggtree::groupOTU(tree_275, branch_strains)



tree_pt<- ggtree(tree_275,
                 branch.length="rate", layout="equal_angle")+xlim(NA,0.3)


ggtree(tree_pt_h,
       branch.length="rate", 
       # layout="equal_angle",
       aes(color=group)) + 
  scale_color_manual(values=ancestry.colours) + 
  theme(legend.position="right")+
  geom_tiplab(align = T)+
  # geom_hilight(1, "steelblue") +
  theme_tree2() 

