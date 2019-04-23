map_file <- data.table::fread("data/ADMIXTURE/BEST_K/LD_0.1_MAF_0.004.map")
p_file <- data.table::fread("data/ADMIXTURE/BEST_K/LD_0.1_MAF_0.004.6.P", header = F, col.names = LETTERS[1:6]) %>%
  dplyr::mutate(CHROM = map_file$V1, POS = map_file$V4) %>%
  dplyr::mutate(CHROM = gsub("23","X",CHROM))

ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

p_df <- p_file %>%
  dplyr::select(B, f = F, CHROM, POS) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sum_freq= sum(B,f),
                low_freq = ifelse(sum(B,f) < 0.01, "BOTH_LOW","BOTH_HIGH"),
                delta_freq = abs(B-f)) %>%
  tidyr::gather(population, frequency, -CHROM, -POS,-low_freq,-delta_freq,-sum_freq) %>%
  dplyr::filter(low_freq!="BOTH_LOW", sum_freq > 0.5)

p_df%>%
  ggplot() +
  aes(x = POS, y = delta_freq, color = sum_freq) +
  geom_point()+
  facet_grid(low_freq~CHROM, scales = "free_x", space = "free")+
  scale_color_viridis_c()
                
p_df%>%
  dplyr::filter(population == "f") %>%
  ggplot() +
  aes(x = POS, y = frequency, color = sum_freq) +
  geom_point()+
  facet_grid(low_freq~CHROM, scales = "free_x", space = "free")+
  scale_color_viridis_c()



test_df <- p_file %>%
  dplyr::select(B, f = F, CHROM, POS) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sum_freq= sum(B,f),
                low_freq = ifelse(sum(B,f) < 0.01, "BOTH_LOW","BOTH_HIGH"),
                delta_freq = abs(B-f)) %>%
  tidyr::gather(population, frequency, -CHROM, -POS,-low_freq,-delta_freq,-sum_freq) %>%
  dplyr::filter(low_freq!="BOTH_LOW", sum_freq > 1, delta_freq < 0.25) %>%
  tidyr::spread(population, frequency)


test_df %>%
  ggplot()+
  aes(x = POS) +
  geom_histogram() +
  facet_grid(.~CHROM, space = "free")
