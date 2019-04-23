# Setup
library(cowplot)
library(ggnetwork)
library(intergraph)
library(igraph)

# load fulcrum cso df
load('data/fulcrum/df.Rda')


cso_species <- cso %>%
  dplyr::mutate(spp_id = ifelse(spp_id == "Heterhabditis zealandica", "Heterohabditis zealandica", spp_id)) %>%
  dplyr::filter(!is.na(spp_id), spp_id != "no match", spp_id != "Unknown", spp_id != "unknown")%>%
  dplyr::select(spp_id, c_label)%>%
  dplyr::group_by(c_label)%>%
  dplyr::summarise(species_ct = length(unique(spp_id)))%>%
  dplyr::filter(species_ct > 1)

species_found_together <- cso%>%
  dplyr::filter(c_label %in% cso_species$c_label)%>%
  dplyr::select(c_label, spp_id)%>%
  dplyr::filter(!is.na(spp_id), spp_id != "Unknown")%>%
  dplyr::group_by(c_label)

species_found_together$spp_id <- gsub(" ", 
                                      "\n",
                                      species_found_together$spp_id)

# function to build input for network from relational observations
# df contains a grouping column, relationship column, observation number column
# ~ ~ #  grouping column - is used to group by and look for relationships
# ~ ~ #  relationship column - end up being the nodes of the network (can have multiple relationship columns)
# repeat_observations - count observations of the same individual
generate_relational_network_df <- function(df, 
                                           group_col, 
                                           relationship_col, 
                                           repeat_observations) {
  
  # process input variables 
  define_range <- df %>%
    dplyr::group_by(!!group_col) %>%
    dplyr::mutate(observation = 1:n())%>%
    dplyr::ungroup() %>%
    dplyr::mutate(dimensions = length(unique((!!relationship_col)))) %>%
    dplyr::select( .,
                   groupings = !!group_col, 
                   rel = !!relationship_col, 
                   observation, 
                   dimensions)
  
  if(!repeat_observations){
    define_range <- define_range %>%
      dplyr::distinct(groupings, rel, .keep_all =T)
  }
  
  # initialize matrix to fill with relationships
  empty_net <- matrix(nrow=length(unique(define_range$rel)),
                      ncol=length(unique(define_range$rel)),
                      dimnames = list(unique(define_range$rel),
                                      unique(define_range$rel)))
  
  # iterate through unique observations in relationship variable 
  for(i in unique(define_range$rel)){
    
    subet_groups <- dplyr::filter(define_range, rel == i)%>%
      dplyr::distinct(groupings, rel,.keep_all=T)%>%
      dplyr::mutate(id = paste(groupings, rel, observation, sep = "_"))
    
    found_together <- dplyr::filter(define_range, groupings %in% subet_groups$groupings)%>%
      dplyr::mutate(id = paste(groupings, rel, observation, sep = "_"))%>%
      dplyr::filter(!(id %in% subet_groups$id))
    
    for(j in found_together$rel){
      if(is.na(empty_net[j,i])){
        empty_net[j,i] <- 1
      } 
      else {
        empty_net[j,i] <- empty_net[j,i] + 1
      }
    }
  }
  
  empty_net[is.na(empty_net)] <-0
  
  full_net <- empty_net
  
  return(full_net)
}

function_net <- generate_relational_network_df(species_found_together, 
                                               quo(c_label), 
                                               quo(spp_id), 
                                               repeat_observations = F)

net=graph.adjacency(function_net,
                    mode="undirected",
                    weighted=TRUE,
                    diag=FALSE)


ggplot(ggnetwork::ggnetwork(net), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodes(color = "black", size = 8) +
  geom_nodelabel_repel(aes( label = vertex.names ),
                       fontface = "bold.italic", box.padding = unit(1, "lines"))+
  theme_blank()+
  geom_edgetext(aes(label = weight), color = "grey25") 

ggsave('plots/CoHabitat_Network.pdf', width = 5, height = 7)


c_found_tog <- species_found_together %>%
  dplyr::filter(grepl("C\\.", spp_id) | grepl("Oscheius", spp_id))

function_net <- generate_relational_network_df(c_found_tog, 
                                               quo(c_label), 
                                               quo(spp_id), 
                                               repeat_observations = F)

net=graph.adjacency(function_net,
                    mode="undirected",
                    weighted=TRUE,
                    diag=FALSE)


ggplot(ggnetwork::ggnetwork(net), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodes(color = "black", size = 8) +
  geom_nodelabel_repel(aes( label = vertex.names ),
                       fontface = "bold.italic", box.padding = unit(1, "lines"))+
  theme_blank()+
  geom_edgetext(aes(label = weight), color = "grey25") 

ggsave('plots/Figure5C.pdf', width = 5, height = 7)
