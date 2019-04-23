#!/usr/bin/env Rscript

library(pophelper)
library(tidyverse)
library(ggthemes)

# hawaii strains
strain_islands <- c("XZ1514" = "#E69F00", "XZ1516" = "#E69F00","XZ1513" = "#E69F00","ECA372" = "#E69F00","ECA701" = "#E69F00","XZ1515" = "#E69F00",
                    "CB4856" = "#56B4E9",
                    "ECA369" = "#009E73","ECA738" = "#009E73",
                    "QX1792" = "#0072B2", "QX1794" = "#0072B2", "QX1793" = "#0072B2", "QX1791" = "#0072B2", "ECA740" = "#0072B2", "ECA741" = "#0072B2", "ECA363" = "#0072B2", "ECA743" = "#0072B2", "ECA742" = "#0072B2",
                    "ECA760" = "#CC79A7","ECA768" = "#CC79A7","ECA777" = "#CC79A7","ECA706" = "#CC79A7","ECA705" = "#CC79A7","ECA703" = "#CC79A7","ECA807" = "#CC79A7","ECA778" = "#CC79A7",
                    "ECA812" = "#CC79A7","ECA710" = "#CC79A7","ECA744" = "#CC79A7","ECA745" = "#CC79A7","ECA732" = "#CC79A7","ECA733" = "#CC79A7","ECA746" = "#CC79A7","DL238" = "#CC79A7",
                    "ECA347" = "#CC79A7","ECA730" = "#CC79A7","ECA724" = "#CC79A7","ECA722" = "#CC79A7","ECA189" = "#CC79A7","ECA191" = "#CC79A7","ECA723" = "#CC79A7","ECA712" = "#CC79A7",
                    "ECA396" = "#CC79A7")

# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")


# generate sample list
# vcf name says 330, but it is a lie
system("bcftools query -l data/ANNOTATE_VCF/Ce330_annotated.vcf.gz > data/ANNOTATE_VCF/samples.txt")

# strain names
# get sample information
sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples.txt", header = F) %>% dplyr::pull(V1))

# plot panel A, K by CV summary
k_summary <- data.table::fread("data/ADMIXTURE_LD8/CV_Summary/admix_replicates_CV.tsv",header = T) 

ksum_plot <- ggplot(k_summary)+
  aes(x = factor(K), y = CV)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = .1)+
  theme_bw()+
  labs(x = "K")


# generate K summary plot - Supplemental figure XX
admix_plots <- list()
for(kpops in 1:length(grep(".Q", list.files("data/ADMIXTURE_LD8/BEST_K/"), value = T))){
  K <- as.numeric(strsplit(grep(".Q", list.files("data/ADMIXTURE_LD8/BEST_K/"), value = T)[kpops], split = "\\.")[[1]][4])
  
  # load Q files
  qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files("data/ADMIXTURE_LD8/BEST_K/"))
  qfile <- pophelper::readQ(files = paste0("data/ADMIXTURE_LD8/BEST_K/",qfile_name))[[1]]
  # add pop names
  colnames(qfile) <- LETTERS[1:K]
  
  qfile <- qfile %>%
    dplyr::mutate(samples = sample_names)
  
  write.table(qfile, file = glue::glue("data/ADMIXTURE_LD8/BEST_K/K{K}_Processed_Ancestry.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
  
  # make long and determin order of plotting
  long_admix_pops <- qfile %>%
    dplyr::mutate(samples = sample_names) %>%
    tidyr::gather(cluster, frac_cluster, -samples) %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(max_frac = max(frac_cluster)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  # establish plot order of strains based on anc pop and max fraction
  plot_order <- long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  admix_plots[[kpops]] <- long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
    theme_bw() +
    theme(axis.text.x=element_blank(),    
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
  
  # extract representative strains from each ancesteral population for generating neighbor-net
  if(!exists("representative_K_strains")){
    representative_K_strains <- dplyr::filter(plot_order, frac_cluster > 0.999, !samples %in% names(strain_islands)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(sample_n = 1:n()) %>%
      dplyr::top_n(3, sample_n) %>%
      dplyr::mutate(K_size = K)
  } else {
    representative_K_strains <- dplyr::filter(plot_order, frac_cluster > 0.999, !samples %in% names(strain_islands)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(sample_n = 1:n()) %>%
      dplyr::top_n(3, sample_n) %>%
      dplyr::mutate(K_size = K) %>%
      dplyr::bind_rows(representative_K_strains, .)
  }
  
}

# make panel B
admixture_plots <- cowplot::plot_grid(admix_plots[[1]],
                   admix_plots[[2]],
                   admix_plots[[3]],
                   admix_plots[[4]],
                   admix_plots[[5]],
                   ncol = 1)

# make final figure
ksummary_plot <- cowplot::plot_grid(ksum_plot,
                   admixture_plots,
                   ncol = 2,
                   labels = c("A", "B"),
                   rel_widths = c(0.5, 1))

# get big legend
admix_legend <- cowplot::plot_grid(cowplot::get_legend(admix_plots[[5]]))

# save
ggsave(ksummary_plot, filename = "plots/SuppFigureXX_KSummary_LD8.pdf", height = 8, width = 12)
ggsave(ksummary_plot, filename = "plots/SuppFigureXX_KSummary_LD8.png", height = 8, width = 12,dpi=300)
# ggsave(final_plot, filename = "plots/Figure6_LD8.png", height = 10, width = 12, dpi = 300)
# ggsave(admix_legend, filename = "plots/admix_legend_LD8.pdf", height = 8, width = 12)
# generate input for running SplitsTree


strain_vector <- paste(sort(names(strain_islands)), sep = ",", collapse = ",")
system(glue::glue("bcftools view -s {strain_vector} data/ANNOTATE_VCF/Ce330_annotated.vcf.gz -Oz -o data/ANNOTATE_VCF/Hawaii.vcf.gz"))
# Generate nexus file, .nexus file is used for SplitsTree
system(glue::glue("python scripts/vcf2phylip.py -i data/ANNOTATE_VCF/Hawaii.vcf.gz -m 43 --fasta --nexus --nexus-binary"))

for(apops in unique(representative_K_strains$K_size)){
  # extract strain names corresponding to different K sizes
  koi <- dplyr::filter(representative_K_strains, K_size == apops)
  strain_vector <- paste(sort(c(as.character(koi$samples), names(strain_islands))), sep = ",", collapse = ",")
  # extract number of strain variable
  n_strains <- length(sort(c(as.character(koi$samples), names(strain_islands))))
  # Subset VCF to contain strains of interest - Hawaiian and representative strains
  system(glue::glue("bcftools view -s {strain_vector} data/ANNOTATE_VCF/Ce330_annotated.vcf.gz -Oz -o data/ANNOTATE_VCF/Hawaii_K{apops}.vcf.gz"))
  # Generate nexus file, .nexus file is used for SplitsTree
  system(glue::glue("python scripts/vcf2phylip.py -i data/ANNOTATE_VCF/Hawaii_K{apops}.vcf.gz -m {n_strains} --fasta --nexus --nexus-binary"))
}

# load .nexus files into SplitsTree and output plot svg into plots folder

write.table(representative_K_strains, 
            file = "data/ADMIXTURE/BEST_K/Hawaiia_plus_RepStrains_by_K.tsv",
            quote = F,
            col.names = T, 
            row.names = F, 
            sep = "\t")










