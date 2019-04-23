library(tidyverse)
load("data/HAPLOTYPE/haplotype_plot_df.Rda")

get_group_number = function(){
  i = 0
  function(){
    i <<- i+1
    i
  }
}
group_number = get_group_number()

strain_islands <- c("XZ1514" = "#E69F00", "XZ1516" = "#E69F00","XZ1513" = "#E69F00","ECA372" = "#E69F00","ECA701" = "#E69F00","XZ1515" = "#E69F00",
                    "CB4856" = "#56B4E9","ECA928" = "#56B4E9","ECA923" = "#56B4E9",
                    "ECA369" = "#009E73","ECA738" = "#009E73",
                    "QX1792" = "#0072B2", "QX1794" = "#0072B2", "QX1793" = "#0072B2", "QX1791" = "#0072B2", "ECA740" = "#0072B2", "ECA741" = "#0072B2", "ECA363" = "#0072B2", "ECA743" = "#0072B2", "ECA742" = "#0072B2",
                    "ECA760" = "#CC79A7","ECA768" = "#CC79A7","ECA777" = "#CC79A7","ECA706" = "#CC79A7","ECA705" = "#CC79A7","ECA703" = "#CC79A7","ECA807" = "#CC79A7","ECA778" = "#CC79A7",
                    "ECA812" = "#CC79A7","ECA710" = "#CC79A7","ECA744" = "#CC79A7","ECA745" = "#CC79A7","ECA732" = "#CC79A7","ECA733" = "#CC79A7","ECA746" = "#CC79A7","DL238" = "#CC79A7",
                    "ECA347" = "#CC79A7","ECA730" = "#CC79A7","ECA724" = "#CC79A7","ECA722" = "#CC79A7","ECA189" = "#CC79A7","ECA191" = "#CC79A7","ECA723" = "#CC79A7","ECA712" = "#CC79A7",
                    "ECA396" = "#CC79A7")

mcolor_grp <- plot_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
mcolor <- mcolor_grp$color
names(mcolor) <- mcolor_grp$haplotype

strain_labels <- plot_df %>%
  dplyr::select(isotype, plotpoint)


pr_plot_df <- plot_df %>%
  dplyr::filter(isotype %in% names(strain_islands)) %>%
  dplyr::group_by(chromosome, isotype) %>%
  dplyr::mutate(total_swept_haplotype = sum(isotype_swept_haplotype_length, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(total_swept_haplotype), plotpoint) %>%
  dplyr::mutate(f_isotype = factor(isotype, levels = unique(isotype))) %>%
  dplyr::group_by(f_isotype) %>%
  dplyr::mutate(plot_order = group_number()) %>%
  dplyr::ungroup() 

pr_plot_df%>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plot_order - 0.5, ymax = plot_order + 0.5,
             fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = pr_plot_df$plot_order,
                     labels = pr_plot_df$isotype,
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme(legend.position="none",
        strip.background.x = element_blank())

ggsave("plots/haplotype.pdf", height = 8, width = 12)





# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  Generate Population Summary File for TREEMIX  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# define K
K = 7

# define outgroup
outgroup_strain <- "XZ1516"

sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples.txt", header = F) %>% dplyr::pull(V1))

# load P files
pfile_name <- grep(pattern = glue::glue("{K}\\.P$"), value = T, x = list.files("data/ADMIXTURE_LD8/BEST_K/"))
pfile <- pophelper::readQ(files = paste0("data/ADMIXTURE_LD8/BEST_K/",pfile_name))[[1]]

# label P file rownames and colnames
colnames(pfile) <- LETTERS[1:K]


treemix_input <- apply(pfile, MARGIN = c(1,2), function(x){
  f <- round(x*length(sample_names),digits = 0)
  paste(length(sample_names)-f, f, sep = ",")
})

# find outgroup population
# load Q files
qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files("data/ADMIXTURE_LD8/BEST_K/"))
qfile <- pophelper::readQ(files = paste0("data/ADMIXTURE_LD8/BEST_K/",qfile_name))[[1]]
# add pop names
colnames(qfile) <- LETTERS[1:K]
qfile$strain <- sample_names

outgroup_population <- qfile%>%
  dplyr::filter(strain == outgroup_strain)%>%
  tidyr::gather(Population, Frequency, -strain)%>%
  dplyr::filter(Frequency == max(Frequency))


write.table(x = treemix_input,
            file = glue::glue("data/ADMIXTURE_LD8/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
            quote = F,
            col.names = T,
            row.names = F)


treemix_input_name <- strsplit(glue::glue("data/ADMIXTURE_LD8/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
                               split = "/")[[1]][3]

# setwd(paste0(base.dir, "data/"))

system(glue::glue("gzip data/ADMIXTURE_LD8/{treemix_input_name}"))

outgroup_pop <- outgroup_population$Population[1]

output_name <- strsplit(strsplit(glue::glue("data/ADMIXTURE_LD8/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
                                 split = "/")[[1]][3], split = "\\.txt")[[1]][1]

system(glue::glue("treemix -i data/ADMIXTURE_LD8/{treemix_input_name}.gz -o data/ADMIXTURE_LD8/{output_name} -root {outgroup_pop} -k 500"))

source("scripts/PLOT_TREEMIX.R")
# 
# treemix.dir <- paste0(base.dir, "Data/Treemix/")
# 
# dir.create(paste0(plot.dir, "/Treemix"), showWarnings = F)
# dir.create(paste0(plot.dir, "/Treemix/Residual_Plots"), showWarnings = F)
# 
# 
# setwd(treemix.dir)

system(paste0(paste0("printf \"", paste(LETTERS[1:K], collapse = '\n')), '\n\" > data/ADMIXTURE_LD8/poporder'))

treemix_stem <- glue::glue("data/ADMIXTURE_LD8/{output_name}")

pdf(glue::glue("plots/K{K}_TREEMIX_phylogeny_NO_MIGRATION.pdf"))
plot_tree(stem = treemix_stem)
dev.off()

pdf(glue::glue("plots/K{K}_TREEMIX_phylogeny_NO_MIGRATION_Residuals.pdf"))
plot_resid(stem = treemix_stem, "data/ADMIXTURE_LD8/poporder")
dev.off()

# treemix with migration

# run treemix with migration -
setwd(treemix.dir)
dir.create("migrations")

for(m in 1:4){
  system(glue::glue("treemix -i data/ADMIXTURE_LD8/{treemix_input_name}.gz -o data/ADMIXTURE_LD8/{output_name}_{m} -root {outgroup_pop} -k 500 -m {m} -se"))
}

# # # # #  PLOT TREEMIX WITH MIGRATIONS

for(migration in gsub("\\.llik","",grep(glue::glue("K-{K}"), grep("llik",list.files("data/ADMIXTURE_LD8/"),value = T), value = T))){
  
  system(paste0(paste0("printf \"", paste(LETTERS[1:K], collapse = '\n')), '\n\" > data/ADMIXTURE_LD8/poporder'))
  
  pdf(glue::glue("plots/TREEMIX_phylogeny_{migration}_migrations.pdf"))
  plot_tree(stem =  glue::glue("data/ADMIXTURE_LD8/{migration}"))
  dev.off()
  
  pdf(glue::glue("plots/TREEMIX_phylogeny_{migration}_migrations_RESIDUALS.pdf"))
  plot_resid(stem = glue::glue("data/ADMIXTURE_LD8/{migration}"), "data/ADMIXTURE_LD8/poporder")
  dev.off()
  # 
}

