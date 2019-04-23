library(tidyverse)
library(lubridate)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# looking at numbers sterile
df1 <- cso %>%
  dplyr::filter(state_of_plate_evanston %in% c("sterile", "Sterile")) %>%
  dplyr::filter(!is.na(spp_id))
