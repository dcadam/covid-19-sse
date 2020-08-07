library(tidyverse)
library(igraph)

case_data <- read_csv(file = "data/case_data.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

quarantine_list <-  case_data %>%
  dplyr::select(case.no, quarantine) %>%
  mutate(id = as.character(case.no)) %>%
  dplyr::select(id, quarantine)

temple_cluster <- read_csv(file = "data/temple_cluster.csv")
wedding_cluster <- read_csv(file = "data/wedding_cluster.csv")
bar_cluster <- read_csv(file = "data/bar_cluster.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

complete_pairs <- transmission_pairs %>%
  filter(cluster.id != 46,
         cluster.id != 77,
         cluster.id != 80) %>%
  mutate(from = infector.case, to = infectee.case) %>%
  dplyr::select(from, to, cluster.risk, cluster.generation) %>%
  full_join(., temple_cluster) %>%
  full_join(., wedding_cluster) %>%
  full_join(., bar_cluster)


###ENCODE CASE POSITIONS WITHIN THE ENTIRE NETWORK AS SOURCE, INTERMEDIATE, OR TERMINAL
infectee <- complete_pairs %>%
  transmute(infector.case = from,
            infectee.case = to) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector <-  complete_pairs %>%
  transmute(infector.case = from,
            infectee.case = to) %>%
  gather() %>%
  filter(key == 'infector.case')

middle_cases <- infector %>%
  left_join(., infectee, by = 'value') %>%
  filter(key.y != 'NA') %>%
  dplyr::select(value) %>%
  distinct() %>%
  transmute(case.no = value, 
            position = "middle")

terminal_cases <- infectee %>% 
  dplyr::select(value) %>%
  filter(!value %in% middle_cases$position) %>%
  distinct() %>%
  transmute(case.no = value, 
            position = "terminal")

sources <- infector %>% 
  dplyr::select(value) %>%
  filter(!value %in% middle_cases$case.no) %>%
  distinct() %>%
  transmute(case.no = value, 
            position = "source")

#find those whole are terminal but are linked to the bar directly and exclude them due to their ambigious position in the chain
bars_direct <- complete_pairs %>%
  transmute(infector.case = from,
            infectee.case = to) %>%
  filter(infector.case == "BARS") %>%
  dplyr::select(infectee.case)

terminal_cases_xbar <- terminal_cases %>%
  filter(!case.no %in% bars_direct$infectee.case) %>%
  distinct() %>%
  transmute(case.no = case.no, 
            position = "terminal")

#combine source, middle, and terminal adding a position coulmn, then terminal Y/N, then quarantine column 
position_list <- sources %>%
  full_join(middle_cases) %>%
  full_join(terminal_cases_xbar)

#excluding cluster sources cases, cross table and count
position_list %>%
  filter(position != "source") %>%
  mutate(terminal = case_when(position == "terminal" ~ "Y",
                              TRUE ~ "N")) %>%
  mutate(id = case.no) %>%
  dplyr::select(-case.no) %>%
  left_join(quarantine_list, by = "id") %>%
  mutate(quarantine = case_when(is.na(quarantine) ~ "N", 
                                TRUE ~ quarantine)) %>%
  group_by(terminal, quarantine) %>%
  count()

#odds ratio
(45/144)/(1/46)




  
