library(tidyverse)
library(lubridate)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
load('data/fulcrum/df.Rda')

# Look at which samples were putatively isolated in Evanston
evanston <- df%>%
 dplyr::filter(po_in_evanston == "TRUE") %>%
dplyr::distinct(c_label)

# Calculate time diffs between collection and isolation 
times <- df %>%
  dplyr::select(c_label, po_date, po_time, datetime, po_in_evanston) %>%
  dplyr::filter(!is.na(po_date)) %>%
  tidyr::unite(col = isolation_datetime, po_date, po_time, sep = " ", remove = F) %>%
  dplyr::mutate(collection_datetime = lubridate::ymd_hms(datetime, tz = "HST"),
                isolation_datetime = lubridate::ymd_hms(isolation_datetime, tz = "HST")) %>%
  dplyr::select(c_label, collection_datetime, isolation_datetime, po_in_evanston) %>%
  dplyr::mutate(time_diff_h = as.numeric(difftime(isolation_datetime, collection_datetime, units = "hours"))) %>%
  dplyr::group_by(po_in_evanston) %>%
  dplyr::mutate(mean_time_diff_h = mean(time_diff_h),
                std_time_diff_h = sd(time_diff_h))

# Plot time histograms
time_diffs <- ggplot(times) +
  aes(x = time_diff_h, fill = po_in_evanston) +
  geom_histogram() +
  theme(axis.title = element_text(size = 9, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9, color = "black"),
        legend.title = element_text(size = 9, color = "black")) +
  labs(fill = "Isolated at NU", x = "Time from collection to isolation (h)", y = "count")

time_diffs
ggsave('plots/time_diffs.pdf', height = 5, width = 5)
ggsave('plots/time_diffs.png', height = 5, width = 5)
