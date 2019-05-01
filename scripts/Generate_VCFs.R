
library(tidyverse)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# ARGS
# 1 - VCF path
VCF_PATH="data/ANNOTATE_VCF/Ce330_annotated.vcf.gz"
# 2 - LD to prune
ld=0.8
# 3 - MAF
maf=0.004
# 3 - output path
OUTPUT_PATH="data/ANNOTATE_VCF/"

# generate markers to be kept
system(glue::glue("plink --vcf {VCF_PATH} --snps-only --biallelic-only --maf {maf} --set-missing-var-ids @:# --indep-pairwise 50 10 {ld} --allow-extra-chr --out {OUTPUT_PATH}LD_{ld}"))
# generate plink LD pruned VCF
system(glue::glue("plink --vcf {VCF_PATH} --snps-only --biallelic-only --maf {maf} --set-missing-var-ids @:# --extract {OUTPUT_PATH}LD_{ld}.prune.in --recode vcf --allow-extra-chr --out {OUTPUT_PATH}LD_{ld}"))
# compress vcf

# fix sample names and compress
system(glue::glue("bcftools query -l {VCF_PATH} > {OUTPUT_PATH}LD_{ld}_samples.txt"))
system(glue::glue("bcftools reheader -s {OUTPUT_PATH}LD_{ld}_samples.txt {OUTPUT_PATH}LD_{ld}.vcf | bcftools view -Oz -o {OUTPUT_PATH}LD_{ld}.vcf.gz"))

# rm uncompressed vcf
system(glue::glue("rm {OUTPUT_PATH}LD_{ld}.vcf"))
