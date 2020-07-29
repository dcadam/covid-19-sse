library(tidyverse)
library(lubridate)
library(ggsci)
library(forcats)

#import data
case_data <- read_csv(file = "data/case_data.csv")

case_data <- case_data %>%
  mutate(onset.date = dmy(onset.date),
         confirm.date = dmy(confirm.date),
         epi.date = dmy(epi.date))

#Descriptive Cluster analysis
case_data %>%
  group_by(cluster.category) %>%
  summarise(n = n(), pct = n/1038)

#find overseas cluster values
case_data %>%
  filter(cluster.category == "Cluster of imported cases") %>%
  group_by(cluster.id) %>%
  n_groups()

case_data %>%
  filter(cluster.category == "Cluster of imported cases") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  median()

case_data %>%
  filter(cluster.category == "Cluster of imported cases") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  range()


#find local cluster values
case_data %>%
  filter(cluster.category == "Cluster initiated by local case") %>%
  group_by(cluster.id) %>%
  n_groups()

case_data %>%
  filter(cluster.category == "Cluster initiated by local case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  median()

case_data %>%
  filter(cluster.category == "Cluster initiated by local case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  range()

#find imported cluster values
case_data %>%
  filter(cluster.category == "Cluster initiated by imported case") %>%
  group_by(cluster.id) %>%
  n_groups()

case_data %>%
  filter(cluster.category == "Cluster initiated by imported case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  median()

case_data %>%
  filter(cluster.category == "Cluster initiated by imported case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  range()

#Plot Figure 1
ggplot(data = case_data) +
  geom_bar(aes(x = epi.date, fill = cluster.category), colour = 'black', size = 0.05, binwidth = 1, alpha = 0.9) +
  labs(x="Illness Onset Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'right') +
  scale_x_date(date_breaks = "2 days", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_y_continuous("Daily SARS-CoV-2 Cases in Hong Kong (N)", expand = c(0,0), limits = c(0,50)) +
  scale_fill_brewer(type = "div", palette = 4)




