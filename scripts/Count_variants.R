library(tidyverse)
library(ggthemes)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# strain names for Hawaii Manuscript 276set
# get sample information for 276 strain set 
sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples.txt", header = F) %>% dplyr::pull(V1))
strain_vector <- paste(sort(sample_names), sep = ",", collapse = ",")

# get just hawaii 43 strain subset
hi_sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples_HI_ONLY.txt", header = F) %>% dplyr::pull(V1))
hi_strain_vector <- paste(sort(hi_sample_names), sep = ",", collapse = ",")

# get just hawaii 26 strain subset (new collections)
hi_new_strain_vector <- c("ECA701","ECA703","ECA705","ECA706","ECA710",
                          "ECA712","ECA722","ECA723","ECA724","ECA730",
                          "ECA732","ECA733","ECA738","ECA740","ECA741",
                          "ECA742","ECA743","ECA744","ECA745","ECA746",
                          "ECA760","ECA768","ECA777","ECA778","ECA807",
                          "ECA812")

# Filter VCFs to variant sets for counting
# Filter 330 set VCF to include just 276 strains for Hawaii Manuscript and exclude sites where all strains are REF
system(glue::glue("bcftools view -s {strain_vector} data/ANNOTATE_VCF/WI.20180527.hard-filter.vcf.gz | bcftools view -c 1 -Oz -o data/ANNOTATE_VCF/276_set_hard-filter.vcf.gz"))

# Filter 276 set VCF to include just 43 HI strains and exclude sites where all strains are REF
system(glue::glue("bcftools view -s {hi_strain_vector} data/ANNOTATE_VCF/276_set_hard-filter.vcf.gz | bcftools view -c 1 -Oz -o data/ANNOTATE_VCF/43_set_HI_only_hard-filter.vcf.gz"))

# Filter 276 set VCF to include just 26 new HI strains and exclude sites where all strains are REF
system(glue::glue("bcftools view -s {hi_new_strain_vector} data/ANNOTATE_VCF/WI.20180527.hard-filter.vcf.gz | bcftools view -c 1 -Oz -o data/ANNOTATE_VCF/26_new_HI_hard-filter.vcf.gz"))

# Count variants in VCFs
# count variants in 276 strain VCF (2,805,683 SNV and small INDELS)
system(glue::glue("bcftools view -H data/ANNOTATE_VCF/276_set_hard-filter.vcf.gz | wc -l")) <- var_count_276

# count variants in HI only 43 strain VCF (2,336,143 SNV and small INDELS)
system(glue::glue("bcftools view -H data/ANNOTATE_VCF/43_set_HI_only_hard-filter.vcf.gz | wc -l")) <- var_count_43

#count variants in new HI 26 strain VCF (526,492 SNV and small INDELS)
system(glue::glue("bcftools view -H data/ANNOTATE_VCF/26_new_HI_hard-filter.vcf.gz | wc -l"))

