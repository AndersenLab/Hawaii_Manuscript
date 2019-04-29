library(ggtern)


# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# set color palette
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

# assign Hawaii isotypes
hi_only_samples <- read.csv(file = "data/fulcrum/hawaii_isotypes.csv") 

# load admix info for full 276_set
admix <- data.table::fread("data/ADMIXTURE_LD8/BEST_K/K7_Processed_Ancestry.tsv") %>%
  dplyr::rename(isotype = samples) %>%
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
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE"))


#generate ternery plots
plotFG <- admix %>%
   dplyr::filter(!pop_assignment %in% c("A", "B", "C", "E")) %>%
   dplyr::select(isotype, pop_assignment, D, F, G)

A <-  ggtern(data=plotFG,aes(F,D,G,fill=pop_assignment)) +
    scale_fill_manual(values = ancestry.colours) +
    theme_classic() +
    #theme_rgbw() +
    geom_point(size = 3, shape = 21) + 
    tern_limit(T = 1.1, L = 1.1, R = 1.1) +
    labs(x="F",y="D",z="G",title="")
A

ggsave(paste("plots/Ternary_FDG_K=7_LD8_HI.png"), width = 3.75, height = 3.75)
ggsave(paste("plots/Ternary_FDG__K=7_LD8_HI.pdf"), width = 3.75, height = 3.75, useDingbats=FALSE)

  
plotF <- admix %>%
# dplyr::group_by(isotype) %>%
# dplyr::mutate(world = sum(A, B, D, E)) %>%
# dplyr::ungroup() %>%
dplyr::filter(!pop_assignment %in% c("A", "B", "E", "G")) %>%
dplyr::select(isotype, pop_assignment, D, F, C)

B <- ggtern(data=plotF,aes(F,D,C,fill=pop_assignment)) +
  scale_fill_manual(values = ancestry.colours) +
  theme_classic() +
  #theme_rgbw() +
  geom_point(size = 3, shape = 21) + 
  tern_limit(T = 1.1, L = 1.1, R = 1.1) +
  labs(x="F",y="D",z="C",title="")
B

ggsave(paste("plots/Ternary_FDC_K=7_LD8_HI.png"), width = 3.75, height = 3.75)
ggsave(paste("plots/Ternary_FDC__K=7_LD8_HI.pdf"), width = 3.75, height = 3.75, useDingbats=FALSE)

plotG <- admix %>%
  # dplyr::group_by(isotype) %>%
  # dplyr::mutate(world = sum(A, B, D, E)) %>%
  # dplyr::ungroup() %>%
  dplyr::filter(!pop_assignment %in% c("A", "B", "E", "F")) %>%
  dplyr::select(isotype, pop_assignment, D, G, C)

C <- ggtern(data=plotG,aes(G,D,C,fill=pop_assignment)) +
  scale_fill_manual(values = ancestry.colours) +
  theme_classic() +
  #theme_rgbw() +
  geom_point(size = 3, shape = 21) + 
  tern_limit(T = 1.1, L = 1.1, R = 1.1) +
  labs(x="G",y="D",z="C",title="")
C

ggsave(paste("plots/Ternary_GDC_K=7_LD8_HI.png"), width = 3.75, height = 3.75)
ggsave(paste("plots/Ternary_GDC__K=7_LD8_HI.pdf"), width = 3.75, height = 3.75, useDingbats=FALSE)

