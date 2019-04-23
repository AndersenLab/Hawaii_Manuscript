#!/usr/bin/env Rscript
library(tidyverse)
library(ggtree)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#define palette
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "Admixed" = "grey")

# load in wi strain info for lat long coordinates for map
df <- data.table::fread("data/WI_strain_list.csv")
strain_colors <- dplyr::filter(df, state== "Hawaii")

tree_275 <- ape::read.tree(glue::glue("data/tree/249-hawaii_genome.raxml.bestTree"))

tree_pt<- ggtree(tree_275,
                 branch.length="rate",layout="equal_angle")+xlim(NA,0.3)

hawaii_strains <- df%>%
  dplyr::mutate(strain_loc = ifelse(strain %in% strain_colors$isotype, "Hawaii", "Not")) %>%
  dplyr::filter(release %in% c("20160408","20170531") | strain %in% strain_colors$isotype) %>%
  dplyr::filter(reference_strain == 1) %>%
  dplyr::select(isotype, lat = latitude, long = longitude, strain_loc) %>%
  dplyr::mutate(color = isotype ) %>%
  as.data.frame()


colored_tree <- tree_pt %<+% hawaii_strains + 
  geom_tiplab(aes(color = strain_loc)) +
  scale_color_manual(values = c("#D7263D", "gray40")) + theme_tree2()

branch_strains <- list(CONNECT = dplyr::filter(hawaii_strains, strain_loc == "Hawaii") %>% dplyr::pull(isotype))

tree_pt_h <- ggtree::groupOTU(tree_275, branch_strains)

# define colors
highlight_color <- "#D7263D"
background_color <- "#000F08"

#plot tree with Hawaiian isotypes colored red
hi_tree <- ggtree(tree_pt_h,
       branch.length="rate", 
       layout="equal_angle",
       aes(color=group)) + 
  scale_color_manual(values=c(background_color, highlight_color), 
                     name = "Presence of TALT", 
                     breaks=c("0", "TALT"),
                     labels=c("FALSE", "TRUE")) + 
  theme(legend.position="right")+
  # geom_hilight(1, "steelblue") +
  theme_tree2() 
hi_tree

########### updated tree based with admix branches
#loading admixture data
admix_pops <- readr::read_tsv("~/Dropbox/AndersenLab/Hawaii_manuscript/data/ADMIXTURE/admix_assignments_K=6.tsv")

# assign groups
admix <- list(A = dplyr::filter(admix_pops, pop_assignment == "A") %>% dplyr::pull(isotype),
                B = dplyr::filter(admix_pops, pop_assignment == "B") %>% dplyr::pull(isotype),
                C = dplyr::filter(admix_pops, pop_assignment == "C") %>% dplyr::pull(isotype),
                D = dplyr::filter(admix_pops, pop_assignment == "D") %>% dplyr::pull(isotype),
                E = dplyr::filter(admix_pops, pop_assignment == "E") %>% dplyr::pull(isotype),
                F = dplyr::filter(admix_pops, pop_assignment == "F") %>% dplyr::pull(isotype))
  
# place groups in tree  
admix_groups <- ggtree::groupOTU(tree_275, admix)

# plot tree with admixture groups on branches
admix_tree <- ggtree(admix_groups,
       branch.length="rate", 
       layout="equal_angle",
       aes(color=group)) +
  #geom_tiplab(aes(subset=(label %in% c("ECA812","ECA701")))) +
  #geom_tiplab() +
  scale_color_manual(values=c(ancestry.colours)) +
  theme(legend.position="right") +
  # geom_hilight(1, "steelblue") +
  theme_tree2() 
admix_tree       

# plotting hawaii and admixture trees above one another to compare 

admix_and_hi_trees <- cowplot::plot_grid(hi_tree, admix_tree, ncol = 1)


ggsave(admix_and_hi_trees, filename = "plots/SuppFigX_whole_pop_tree_admix_colors.pdf", height = 7.5, width = 7.5)
ggsave(admix_and_hi_trees, filename = "plots/SuppFigX_whole_pop_tree_admix_colors.png", height = 7.5, width = 7.5, dpi = 300) 

# regular tree so we can see
# plot tree with admixture groups on branches
admix_tree <- ggtree(admix_groups,
                     branch.length="rate",
                     aes(color=group)) +
  #geom_tiplab(aes(subset=(label %in% c("ECA812","ECA701")))) +
  geom_tiplab() +
  scale_color_manual(values=c(ancestry.colours)) +
  theme(legend.position="right") +
  # geom_hilight(1, "steelblue") +
  theme_tree2() 
admix_tree       

ggsave(admix_tree, filename = "plots/SuppFigX_whole_pop_tree_admix_colors_CLEAR.pdf", height = 35, width = 20)
