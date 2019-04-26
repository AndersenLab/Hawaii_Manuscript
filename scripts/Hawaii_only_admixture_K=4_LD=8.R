library(tidyverse)

ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# vcf name says 330, but it is a lie
system("bcftools query -l data/ANNOTATE_VCF/ce_norm_HI_ONLY.vcf.gz > data/ANNOTATE_VCF/samples_HI_ONLY.txt")

# get sample strain name information from HI_ONLY VCF
sample_names <- data.table::fread("data/ANNOTATE_VCF/samples_HI_ONLY.txt", header = F) %>% dplyr::pull(V1)


############################## K=4 HI_ONLY  
# load Q file for K = 4
qfile <- pophelper::readQ(files = "data/ADMIXTURE_LD8_HI_ONLY/BEST_K/LD_0.8_MAF_0.025.4.Q")

# label Q file rownames and colnames
rownames(qfile[[1]]) <- sample_names
K <- 4
colnames(qfile[[1]]) <- LETTERS[1:K]

################################################################
# find admixture assignment for HAwaiian strains               #
################################################################
# assign Hawaii isotypes
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
  tidyr::spread(pop, frac_pop) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(alt_pop_name = case_when(
                                        pop_assignment == "A" ~ "whatA",
                                        pop_assignment == "B" ~ "whatB",
                                        pop_assignment == "C" ~ "whatC",
                                        pop_assignment == "D" ~ "whatD"))

########################
# Plot admix
########################
# make long and determin order of plotting
long_admix_pops <- admix_pops %>%
  dplyr::arrange(pop_assignment, desc(max_pop_frac)) %>%
  dplyr::mutate(isotype = factor(isotype, levels = unique(isotype)))

# establish plot order of strains based on anc pop and max fraction
plot_order <- long_admix_pops %>%
  tidyr::gather(cluster, frac_cluster, -isotype, -max_pop_frac, -alt_pop_name, - pop_assignment) %>%
  dplyr::arrange(isotype)

 admix_plot <-  ggplot(plot_order) +
  geom_bar(stat = "identity", 
           aes(x = isotype, 
               y = frac_cluster, 
               fill = cluster)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1),    
        axis.text.y=element_blank(),
        axis.title.y = element_text(angle = 0, vjust = .5),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

  ggsave(admix_plot, filename = "plots/ADMIXTURE_LD8_HI_ONLY_K=4.pdf", height = 4, width = 12)
  ggsave(admix_plot, filename = "plots/ADMIXTURE_LD8_HI_ONLY_K=4.png", height = 4, width = 12,dpi=300)  
  
############################## K=3 HI_ONLY  
# load Q file for K = 3
  qfile3 <- pophelper::readQ(files = "data/ADMIXTURE_LD8_HI_ONLY/BEST_K/LD_0.8_MAF_0.025.3.Q")
  
# label Q file rownames and colnames
  rownames(qfile3[[1]]) <- sample_names
  K <- 3
  colnames(qfile3[[1]]) <- LETTERS[1:K]
  
  # assign Hawaii isotypes
  admix_pops3 <- qfile3[[1]] %>%
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
    tidyr::spread(pop, frac_pop) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(alt_pop_name = case_when(
      pop_assignment == "A" ~ "whatA",
      pop_assignment == "B" ~ "whatB",
      pop_assignment == "C" ~ "whatC"))
  
  ########################
  # Plot admix3
  ########################
  # make long and determin order of plotting
  long_admix_pops3 <- admix_pops3 %>%
    dplyr::arrange(pop_assignment, desc(max_pop_frac)) %>%
    dplyr::mutate(isotype = factor(isotype, levels = unique(isotype)))
  
  # establish plot order of strains based on anc pop and max fraction
  plot_order3 <- long_admix_pops3 %>%
    tidyr::gather(cluster, frac_cluster, -isotype, -max_pop_frac, -alt_pop_name, - pop_assignment) %>%
    dplyr::arrange(isotype)
  
  admix_plot3 <-  ggplot(plot_order3) +
    geom_bar(stat = "identity", 
             aes(x = isotype, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 90, hjust = 1),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(angle = 0, vjust = .5),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))
  
  ggsave(admix_plot3, filename = "plots/ADMIXTURE_LD8_HI_ONLY_K=3.pdf", height = 4, width = 12)
  ggsave(admix_plot3, filename = "plots/ADMIXTURE_LD8_HI_ONLY_K=3.png", height = 4, width = 12,dpi=300)  
  