library(tidyverse)
library(igraph)

#read data for  analysis
case_data <- read_csv(file = "data/case_data.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

quarantine_list <-  case_data %>%
  dplyr::select(case.no, quarantine) %>%
  mutate(id = as.character(case.no)) %>%
  dplyr::select(id, quarantine)

##BAR CLUSTER GRAPH
#Bar Cluster Import
bar_cluster <- read_csv(file = "data/bar_cluster.csv")

#reframe bar cluster edges to line-listed data and colour by risk
nodes <- enframe(c(bar_cluster$from, bar_cluster$to)) %>%
  dplyr::select(value) %>%
  distinct() %>%
  transmute(id = value,
            label = value,
            to = value) %>%
  left_join(.,bar_cluster, by = 'to') %>%
  dplyr::select(id, cluster.generation, cluster.risk) %>%
  mutate(color = case_when(id == 'NA' ~ '#b43b15',
                           cluster.risk == 'social' ~ '#fbbf28',
                           cluster.risk == 'family' ~ '#2e7d80',
                           cluster.risk == 'work' ~ '#f56890',
                           cluster.risk == 'travel' ~ '#70a494')) %>%
  left_join(., quarantine_list, by = "id") %>%
  mutate(shape = case_when(quarantine == "Y" ~ 'csquare',
                           TRUE ~ 'circle')) %>%
  distinct()

#generate edges for plotting
network <- bar_cluster %>%
  mutate(source = from, target = to) %>%
  dplyr::select(source, target) %>%
  as.data.frame() %>%
  graph_from_data_frame()

#recode network parameters extracting data from line-listed nodes
V(network)$color <- nodes %>% pull(color)
V(network)$size <- 7
V(network)$frame.color <- "black"
V(network)$shape <- nodes %>% pull(shape)

gp <- layout_with_graphopt(network, charge=0.000001, spring.length = 0)

#Plot bar cluster
plot(network, 
     edge.arrow.size=.05,
     vertex.label=NA,
     layout=gp)

##WEDDING CLUSTER GRAPH
wedding_cluster <- read_csv(file = "data/wedding_cluster.csv")

#reframe bar cluster edges to line-listed data and colour by risk
nodes <- enframe(c(wedding_cluster$from, wedding_cluster$to)) %>%
  dplyr::select(value) %>%
  distinct() %>%
  transmute(id = value,
            label = value,
            to = value) %>%
  left_join(.,wedding_cluster, by = 'to') %>%
  dplyr::select(id, cluster.generation, cluster.risk) %>%
  mutate(color = case_when(id == 'NA' ~ '#b43b15',
                           cluster.risk == 'social' ~ '#fbbf28',
                           cluster.risk == 'family' ~ '#2e7d80',
                           cluster.risk == 'work' ~ '#f56890',
                           cluster.risk == 'travel' ~ '#70a494')) %>%
  left_join(., quarantine_list, by = "id") %>%
  mutate(shape = case_when(quarantine == "Y" ~ 'csquare',
                           TRUE ~ 'circle')) %>%
  distinct()

#generate edges for plotting
network <- wedding_cluster %>%
  mutate(source = from, target = to) %>%
  dplyr::select(source, target) %>%
  as.data.frame() %>%
  graph_from_data_frame()

network <- graph_from_data_frame(links)

#recode network parameters extracting data from line-listed nodes
V(network)$color <- nodes %>% pull(color)
V(network)$size <- 10
V(network)$frame.color <- "black"
V(network)$shape <- nodes %>% pull(shape)

#Plot wedding cluster
plot(network, 
     edge.arrow.size=.05,
     vertex.label=NA, 
     layout=layout_with_fr)

##TEMPLE CLUSTER GRAPH
temple_cluster <- read_csv(file = "data/temple_cluster.csv")

#reframe temple cluster edges to line-listed data and colour by risk
nodes <- enframe(c(temple_cluster$from, temple_cluster$to)) %>%
  dplyr::select(value) %>%
  distinct() %>%
  transmute(id = value,
            label = value,
            to = value) %>%
  left_join(.,temple_cluster, by = 'to') %>%
  dplyr::select(id, cluster.generation, cluster.risk) %>%
  mutate(color = case_when(id == 'NA' ~ '#b43b15',
                           cluster.risk == 'social' ~ '#fbbf28',
                           cluster.risk == 'family' ~ '#2e7d80',
                           cluster.risk == 'work' ~ '#f56890',
                           cluster.risk == 'travel' ~ '#70a494')) %>%
  left_join(., quarantine_list, by = "id") %>%
  mutate(shape = case_when(quarantine == "Y" ~ 'csquare',
                           TRUE ~ 'circle')) %>%
  distinct()

#generate edges for plotting
network <- temple_cluster %>%
  mutate(source = from, target = to) %>%
  dplyr::select(source, target) %>%
  as.data.frame() %>%
  graph_from_data_frame()

#recode network parameters extracting data from line-listed nodes
V(network)$color <- nodes %>% pull(color)
V(network)$size <- 15
V(network)$frame.color <- "black"
V(network)$shape <- nodes %>% pull(shape)

#Plot temple cluster
plot(network, 
     edge.arrow.size=.1,
     vertex.label=NA, 
     layout=layout_with_dh)

##OTHER CLUSTER GRAPH

#exclude cases from three previous clusters (cluster.id 46, 77, 80) from the transmission_pairs dataset which includes all resolved pairs
other_clusters <- transmission_pairs %>%
  filter(cluster.id != 46,
         cluster.id != 77,
         cluster.id != 80) %>%
  dplyr::select(infector.case, infectee.case, cluster.risk, cluster.generation, pair.type)

#generate list of imported source only pairs
infectee <- transmission_pairs %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector <-  transmission_pairs %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infector.case')

infectors_only <- infector %>%
  left_join(., infectee, by = 'value') %>%
  filter(is.na(key.y)) %>%
  dplyr::select(value) %>%
  distinct()

source_imported <- other_clusters %>%
  dplyr::select(infector.case, pair.type) %>%
  filter(infector.case %in% infectors_only$value) %>%
  mutate(id = infector.case) %>%
  dplyr::select(-infector.case) %>%
  filter(pair.type == "imported") %>%
  dplyr::select(id, pair.type)

#reframe all other cluster edges to line-listed data and colour by risk and source imported
nodes <- enframe(c(other_clusters$infector.case, other_clusters$infectee.case)) %>%
  dplyr::select(value) %>%
  distinct() %>%
  transmute(id = value,
            label = value,
            infectee.case = value) %>%
  left_join(., other_clusters, by = "infectee.case") %>%
  dplyr::select(id, label, cluster.risk) %>%
  mutate(color = case_when(id == 'NA' ~ '#b43b15',
                           cluster.risk == 'social' ~ '#fbbf28',
                           cluster.risk == 'family' ~ '#2e7d80',
                           cluster.risk == 'work' ~ '#f56890',
                           cluster.risk == 'travel' ~ '#70a494',
                           TRUE ~ '#b43b15')) %>%
  left_join(., source_imported, by = "id") %>%
  mutate(color = case_when(pair.type == 'imported' ~ 'black',
                           TRUE ~ color)) %>%
  left_join(., quarantine_list, by = "id") %>%
  mutate(shape = case_when(quarantine == "Y" ~ 'csquare',
                           TRUE ~ 'circle')) %>%
  dplyr::select(-cluster.risk, -pair.type) %>%
  distinct()
  
  #generate edges for plotting
  network <- other_clusters %>%
    mutate(source = infector.case, target = infectee.case) %>%
    dplyr::select(source, target) %>%
    as.data.frame() %>%
    graph_from_data_frame()
  
  #recode network parameters extracting data from line-listed nodes
  V(network)$color <- nodes %>% pull(color)
  V(network)$size <- 5
  V(network)$frame.color <- "black"
  V(network)$shape <- nodes %>% pull(shape)
  
  #Plot temple cluster
  plot(network, 
       edge.arrow.size=.05,
       vertex.label=NA, 
       layout=layout_with_fr)
  
  

  
